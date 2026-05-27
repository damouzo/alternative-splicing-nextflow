#!/usr/bin/env python3
"""
Generate synthetic test data for alternative-splicing-nextflow validation.

Simulates a Skipped Exon (SE) event between two conditions (3 replicates each):

  ctrl  (n=3)  high cassette-exon inclusion  PSI ≈ 0.80
  treat (n=3)  cassette-exon skipping        PSI ≈ 0.20

Gene model — single gene on chromosome "1" (Ensembl style, no chr prefix):

  Exon 1 (constitutive)  : 1000 – 1100  (101 bp)
  Exon 2 (cassette SE)   : 1200 – 1300  (101 bp)
  Exon 3 (constitutive)  : 1500 – 1600  (101 bp)

  TX1 / ENST00000000001  inclusion isoform  E1-E2-E3  (303 bp)
  TX2 / ENST00000000002  skipping  isoform  E1-E3     (202 bp)

Additionally, 10 background genes (ENSG00000000002–11) are written to the GTF
and quant.sf files to give IsoformSwitchAnalyzeR's sva::num.sv() enough features.
Each background gene has 2 isoforms with equal expression across conditions.
Background genes have no junction reads in the BAM (rMATS ignores them).

Junction reads (100 bp PE, anchors 40+60 or 60+40 bp):

  E1→E2  CIGAR 40M99N60M   pos 1060  (0-based)
  E2→E3  CIGAR 60M199N40M  pos 1240  (0-based)
  E1→E3  CIGAR 40M399N60M  pos 1060  (0-based)  ← skipping

PSI formula  (rMATS JC mode):
  IncLevel = IJC_avg / (IJC_avg + 2 × SJC)

  ctrl  IJC_avg=(80+80)/2=80  SJC=10  PSI=80/100=0.80
  treat IJC_avg=(20+20)/2=20  SJC=40  PSI=20/100=0.20

Outputs under <outdir>/ (default: test/):

  ref/test.genome.fa      single-chromosome FASTA (22000 bp)
  ref/test.gtf            22-transcript gene model (TESTGENE + 10 background)
  samples/<s>/<s>.bam     coordinate-sorted BAM + .bai
  samples/<s>/quant.sf    Salmon-format quantification (for ISAR)
  samplesheet.csv
  comparisons.csv
  params.yaml             ready-to-use params file (run_majiq disabled)

Requires: pysam >= 0.20
  conda install -c bioconda pysam
  pip install pysam

Run:
  python bin/generate_demo_data.py --outdir test

Then launch the full pipeline (rMATS + ISAR, MAJIQ disabled):
  nextflow run main.nf -profile docker,test -params-file test/params.yaml
  nextflow run main.nf -profile apptainer,test -params-file test/params.yaml
"""

import argparse
import csv
import random
import sys
import textwrap
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.exit(
        "ERROR: pysam is required to generate BAM files.\n"
        "  conda install -c bioconda pysam\n"
        "  pip install pysam"
    )

# ── Gene / genome constants ───────────────────────────────────────────────────

CHROM      = "1"
CHROM_LEN  = 22000   # extended to accommodate background genes
READ_LEN   = 100

# Exon boundaries — 1-based, inclusive (standard GTF / SAM coordinates)
E1_START, E1_END = 1000, 1100
E2_START, E2_END = 1200, 1300
E3_START, E3_END = 1500, 1600

TX1_LEN    = (E1_END - E1_START + 1) + (E2_END - E2_START + 1) + (E3_END - E3_START + 1)  # 303
TX2_LEN    = (E1_END - E1_START + 1) + (E3_END - E3_START + 1)                             # 202
# Effective length approximation: drop 50 bp to account for fragment ends
TX1_EFFLEN = TX1_LEN - 50   # 253
TX2_EFFLEN = TX2_LEN - 50   # 152

GENE_ID   = "ENSG00000000001"
TX1_ID    = "ENST00000000001"   # inclusion isoform
TX2_ID    = "ENST00000000002"   # skipping  isoform
GENE_NAME = "TESTGENE"

# Background genes — constitutive expression, no SE event, no junction reads.
# Each has 2 isoforms: a "long" (2-exon) and "short" (1-exon) variant.
# Needed so sva::num.sv() has ≥ 3 features (ndf = min(22,6) - 2 = 4 > 0).
N_BG_GENES = 10
BG_GENE_OFFSET = 2000    # first background gene starts at position 2000 on chr1
BG_GENE_STRIDE = 2000    # 2 kb between gene starts


def _bg_gene(i: int) -> dict:
    """Return coordinates and IDs for background gene i (1-indexed)."""
    base = BG_GENE_OFFSET + (i - 1) * BG_GENE_STRIDE
    # Two exons separated by a 200 bp intron
    ex_a_s, ex_a_e = base + 100, base + 299    # 200 bp
    ex_b_s, ex_b_e = base + 500, base + 699    # 200 bp
    gene_id  = f"ENSG{i + 1:09d}"
    tx_l_id  = f"ENST{(i * 2 + 1):09d}"        # long  (ex_a + ex_b)
    tx_s_id  = f"ENST{(i * 2 + 2):09d}"        # short (ex_a only)
    gene_name = f"BGENE{i:02d}"
    return dict(
        base=base, ex_a_s=ex_a_s, ex_a_e=ex_a_e,
        ex_b_s=ex_b_s, ex_b_e=ex_b_e,
        gene_id=gene_id, gene_name=gene_name,
        tx_l_id=tx_l_id, tx_s_id=tx_s_id,
    )


BG_GENES = [_bg_gene(i) for i in range(1, N_BG_GENES + 1)]

# Sample definitions: (name, condition, replicate)
SAMPLES = [
    ("ctrl_rep1",  "ctrl",  1),
    ("ctrl_rep2",  "ctrl",  2),
    ("ctrl_rep3",  "ctrl",  3),
    ("treat_rep1", "treat", 1),
    ("treat_rep2", "treat", 2),
    ("treat_rep3", "treat", 3),
]

# Junction read counts per replicate: (n_E1→E2, n_E2→E3, n_E1→E3_skip)
# ctrl  PSI: IJC=(80+80)/2=80  SJC=10  → 80/(80+20)=0.80
# treat PSI: IJC=(20+20)/2=20  SJC=40  → 20/(20+80)=0.20
JUNCTION_COUNTS = {
    "ctrl":  [(80, 80, 10), (78, 82,  9), (82, 79, 11)],
    "treat": [(20, 20, 40), (18, 22, 38), (22, 18, 42)],
}

# Salmon NumReads per replicate: (tx1_inclusion, tx2_skipping)
SALMON_COUNTS = {
    "ctrl":  [(800, 200), (820, 180), (790, 210)],
    "treat": [(200, 800), (180, 820), (210, 790)],
}


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--outdir", default="test",
        help="Root output directory (default: test/)",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for genome sequence (default: 42)",
    )
    return p.parse_args()


# ── Genome FASTA ──────────────────────────────────────────────────────────────

def make_genome_seq(seed: int) -> str:
    """Deterministic pseudo-random genome sequence."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=CHROM_LEN))


def write_fasta(path: Path, chrom: str, seq: str) -> None:
    with open(path, "w") as fh:
        fh.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    print(f"  wrote {path}")


# ── GTF annotation ────────────────────────────────────────────────────────────

def _gtf_attr(gene_id: str, **kwargs) -> str:
    parts = [f'gene_id "{gene_id}"']
    parts += [f'{k} "{v}"' for k, v in kwargs.items()]
    return "; ".join(parts) + ";"


def write_gtf(path: Path) -> None:
    # Each row: chrom, source, feature, start, end, score, strand, frame, attributes
    rows = [
        # Gene
        (CHROM, "test", "gene", E1_START, E3_END, ".", "+", ".",
         _gtf_attr(GENE_ID, gene_name=GENE_NAME)),
        # TX1 — inclusion isoform (E1-E2-E3)
        (CHROM, "test", "transcript", E1_START, E3_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX1_ID,
                   gene_name=GENE_NAME, transcript_name=f"{GENE_NAME}-201")),
        (CHROM, "test", "exon", E1_START, E1_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX1_ID, gene_name=GENE_NAME, exon_number="1")),
        (CHROM, "test", "exon", E2_START, E2_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX1_ID, gene_name=GENE_NAME, exon_number="2")),
        (CHROM, "test", "exon", E3_START, E3_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX1_ID, gene_name=GENE_NAME, exon_number="3")),
        # TX2 — skipping isoform (E1-E3)
        (CHROM, "test", "transcript", E1_START, E3_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX2_ID,
                   gene_name=GENE_NAME, transcript_name=f"{GENE_NAME}-202")),
        (CHROM, "test", "exon", E1_START, E1_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX2_ID, gene_name=GENE_NAME, exon_number="1")),
        (CHROM, "test", "exon", E3_START, E3_END, ".", "+", ".",
         _gtf_attr(GENE_ID, transcript_id=TX2_ID, gene_name=GENE_NAME, exon_number="2")),
    ]

    # Background genes — 2 isoforms each (long = 2 exons, short = 1 exon)
    for bg in BG_GENES:
        gid, gname = bg["gene_id"], bg["gene_name"]
        tx_l, tx_s = bg["tx_l_id"], bg["tx_s_id"]
        a_s, a_e   = bg["ex_a_s"], bg["ex_a_e"]
        b_s, b_e   = bg["ex_b_s"], bg["ex_b_e"]
        rows += [
            (CHROM, "test", "gene",       a_s, b_e, ".", "+", ".", _gtf_attr(gid, gene_name=gname)),
            (CHROM, "test", "transcript", a_s, b_e, ".", "+", ".",
             _gtf_attr(gid, transcript_id=tx_l, gene_name=gname, transcript_name=f"{gname}-201")),
            (CHROM, "test", "exon",       a_s, a_e, ".", "+", ".",
             _gtf_attr(gid, transcript_id=tx_l, gene_name=gname, exon_number="1")),
            (CHROM, "test", "exon",       b_s, b_e, ".", "+", ".",
             _gtf_attr(gid, transcript_id=tx_l, gene_name=gname, exon_number="2")),
            (CHROM, "test", "transcript", a_s, a_e, ".", "+", ".",
             _gtf_attr(gid, transcript_id=tx_s, gene_name=gname, transcript_name=f"{gname}-202")),
            (CHROM, "test", "exon",       a_s, a_e, ".", "+", ".",
             _gtf_attr(gid, transcript_id=tx_s, gene_name=gname, exon_number="1")),
        ]

    with open(path, "w") as fh:
        fh.write("##gff-version 2\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}  ({N_BG_GENES} background genes + TESTGENE)")


# ── BAM helpers ───────────────────────────────────────────────────────────────

def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _make_bam_header(sample_name: str) -> dict:
    return {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
        "PG": [{"ID": "generate_demo_data", "PN": "generate_demo_data", "VN": "1.0"}],
        "RG": [{"ID": "1", "SM": sample_name, "LB": "lib1", "PL": "ILLUMINA"}],
    }


def _make_read(
    header,
    name: str,
    flag: int,
    ref_id: int,
    pos: int,           # 0-indexed
    cigar: str,
    mate_pos: int,      # 0-indexed
    seq: str,
) -> "pysam.AlignedSegment":
    a = pysam.AlignedSegment(header)
    a.query_name           = name
    a.flag                 = flag
    a.reference_id         = ref_id
    a.reference_start      = pos
    a.mapping_quality      = 255
    a.cigarstring          = cigar
    a.next_reference_id    = ref_id
    a.next_reference_start = mate_pos
    a.template_length      = 0       # omit TLEN — not used by rMATS/MAJIQ
    a.query_sequence       = seq
    a.query_qualities      = pysam.qualitystring_to_array("I" * len(seq))
    a.set_tag("NH", 1)   # unique mapper
    a.set_tag("NM", 0)   # no mismatches
    return a


# ── BAM writer ────────────────────────────────────────────────────────────────

def write_bam(
    bam_path: Path,
    genome: str,
    n_e1e2: int,
    n_e2e3: int,
    n_e1e3: int,
    sample_name: str,
) -> None:
    """
    Write a coordinate-sorted BAM with junction reads encoding a SE event.

    Read layout (all 100 bp):
      read1 (FLAG 99)  — junction-spanning forward read
      read2 (FLAG 147) — exon-body reverse mate in E3

    Junction CIGAR strings (0-indexed positions):
      E1→E2  40M 99N 60M  start=1060  seq=genome[1060:1100]+genome[1199:1259]
      E2→E3  60M199N 40M  start=1240  seq=genome[1240:1300]+genome[1499:1539]
      E1→E3  40M399N 60M  start=1060  seq=genome[1060:1100]+genome[1499:1559]
    """
    # 0-indexed exon boundaries
    E1_END0   = E1_END          # exclusive end  = 1-based end = 1100
    E2_START0 = E2_START - 1   # 0-indexed start = 1199
    E2_END0   = E2_END          # exclusive end  = 1300
    E3_START0 = E3_START - 1   # 0-indexed start = 1499

    # Read1 start positions (0-indexed)
    POS_E1E2 = E1_END0 - 40         # 1060
    POS_E2E3 = E2_END0 - 60         # 1240
    POS_E1E3 = E1_END0 - 40         # 1060
    POS_MATE = E3_START0             # 1499

    # Intron lengths (number of bases skipped by N in CIGAR)
    N_INTRON1  = E2_START0 - E1_END0          # 1199 - 1100 = 99
    N_INTRON2  = E3_START0 - E2_END0          # 1499 - 1300 = 199
    N_INTRON12 = E3_START0 - E1_END0          # 1499 - 1100 = 399

    CIGAR_E1E2 = f"40M{N_INTRON1}N60M"        # 40M99N60M
    CIGAR_E2E3 = f"60M{N_INTRON2}N40M"        # 60M199N40M
    CIGAR_E1E3 = f"40M{N_INTRON12}N60M"       # 40M399N60M
    CIGAR_MATE = f"{READ_LEN}M"               # 100M

    # Read sequences — extracted from genome at the aligned positions
    seq_e1e2 = genome[POS_E1E2 : E1_END0]  + genome[E2_START0 : E2_START0 + 60]
    seq_e2e3 = genome[POS_E2E3 : E2_END0]  + genome[E3_START0 : E3_START0 + 40]
    seq_e1e3 = genome[POS_E1E3 : E1_END0]  + genome[E3_START0 : E3_START0 + 60]
    seq_mate = _revcomp(genome[POS_MATE : POS_MATE + READ_LEN])

    assert len(seq_e1e2) == READ_LEN, f"E1→E2 seq length {len(seq_e1e2)} ≠ {READ_LEN}"
    assert len(seq_e2e3) == READ_LEN, f"E2→E3 seq length {len(seq_e2e3)} ≠ {READ_LEN}"
    assert len(seq_e1e3) == READ_LEN, f"E1→E3 seq length {len(seq_e1e3)} ≠ {READ_LEN}"
    assert len(seq_mate) == READ_LEN, f"mate seq length {len(seq_mate)} ≠ {READ_LEN}"

    header_dict = _make_bam_header(sample_name)
    ref_id = 0

    # FLAG 99  = 0x1 (paired) | 0x2 (proper) | 0x20 (mate_rev) | 0x40 (read1)
    # FLAG 147 = 0x1 (paired) | 0x2 (proper) | 0x10 (rev)      | 0x80 (read2)
    FLAG_R1 = 0x1 | 0x2 | 0x20 | 0x40   # 99
    FLAG_R2 = 0x1 | 0x2 | 0x10 | 0x80   # 147

    reads = []
    read_num = 0

    def add_pair(r1_pos: int, r1_cigar: str, r1_seq: str) -> None:
        nonlocal read_num
        read_num += 1
        name = f"{sample_name}.{read_num:06d}"
        header = pysam.AlignmentHeader.from_dict(header_dict)
        reads.append(_make_read(header, name, FLAG_R1, ref_id, r1_pos, r1_cigar, POS_MATE, r1_seq))
        reads.append(_make_read(header, name, FLAG_R2, ref_id, POS_MATE, CIGAR_MATE, r1_pos, seq_mate))

    for _ in range(n_e1e2):
        add_pair(POS_E1E2, CIGAR_E1E2, seq_e1e2)
    for _ in range(n_e2e3):
        add_pair(POS_E2E3, CIGAR_E2E3, seq_e2e3)
    for _ in range(n_e1e3):
        add_pair(POS_E1E3, CIGAR_E1E3, seq_e1e3)

    # Sort by (ref, pos) before writing — pysam.index() requires coordinate order
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))

    # Write sorted BAM directly (header already declares SO:coordinate after sort)
    sorted_header_dict = {**header_dict, "HD": {"VN": "1.6", "SO": "coordinate"}}
    tmp_bam = bam_path.with_suffix(".unsorted.bam")
    sorted_header = pysam.AlignmentHeader.from_dict(sorted_header_dict)
    with pysam.AlignmentFile(str(tmp_bam), "wb", header=sorted_header) as out:
        for r in reads:
            # Re-attach to the sorted header so the reference name resolves correctly
            out.write(r)

    pysam.sort("-o", str(bam_path), str(tmp_bam))
    pysam.index(str(bam_path))
    tmp_bam.unlink()

    total_reads = read_num * 2
    print(f"  wrote {bam_path}  ({total_reads} reads, PSI target: "
          f"{'0.80' if n_e1e3 < n_e1e2 else '0.20'})")


# ── Salmon quant.sf ───────────────────────────────────────────────────────────

def _tpm(num_reads: int, eff_len: int, total_rate: float) -> float:
    """TPM = (reads/eff_len) / total_rate * 1e6"""
    return (num_reads / eff_len) / total_rate * 1e6


def write_quant_sf(path: Path, n_tx1: int, n_tx2: int) -> None:
    """
    Write a minimal Salmon quant.sf for TESTGENE isoforms + background genes.

    TPM is computed consistently from NumReads and EffectiveLength so that
    IsoformSwitchAnalyzeR imports the correct isoform fractions.
    Background genes each have 150 reads for the long isoform and 50 for the
    short — constant across conditions so they don't confound the test.
    """
    # Collect all (tx_id, length, eff_len, n_reads) entries
    entries = [
        (TX1_ID, TX1_LEN, TX1_EFFLEN, n_tx1),
        (TX2_ID, TX2_LEN, TX2_EFFLEN, n_tx2),
    ]
    for bg in BG_GENES:
        tx_l, tx_s = bg["tx_l_id"], bg["tx_s_id"]
        a_s, a_e   = bg["ex_a_s"], bg["ex_a_e"]
        b_s, b_e   = bg["ex_b_s"], bg["ex_b_e"]
        len_long   = (a_e - a_s + 1) + (b_e - b_s + 1)   # 400 bp
        len_short  = (a_e - a_s + 1)                       # 200 bp
        eff_long   = max(len_long  - 50, 1)                # 350 bp
        eff_short  = max(len_short - 50, 1)                # 150 bp
        entries.append((tx_l, len_long,  eff_long,  150))
        entries.append((tx_s, len_short, eff_short,  50))

    # Compute TPM jointly over all transcripts
    rates = [n / eff for _, _, eff, n in entries]
    total_rate = sum(rates)

    with open(path, "w") as fh:
        fh.write("Name\tLength\tEffectiveLength\tNumReads\tTPM\n")
        for (tx_id, length, eff_len, n_reads), rate in zip(entries, rates):
            tpm = rate / total_rate * 1e6
            fh.write(f"{tx_id}\t{length}\t{eff_len:.1f}\t{n_reads:.3f}\t{tpm:.6f}\n")


# ── Metadata files ────────────────────────────────────────────────────────────

def write_samplesheet(path: Path, outdir: Path) -> None:
    fieldnames = ["sample", "condition", "replicate", "bam", "bai", "salmon_dir"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for sample, condition, replicate in SAMPLES:
            sample_dir = outdir / "samples" / sample
            writer.writerow({
                "sample":     sample,
                "condition":  condition,
                "replicate":  replicate,
                "bam":        str(sample_dir / f"{sample}.bam"),
                "bai":        str(sample_dir / f"{sample}.bam.bai"),
                "salmon_dir": str(sample_dir),
            })
    print(f"  wrote {path}")


def write_comparisons(path: Path) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["group1", "group2"])
        writer.writeheader()
        writer.writerow({"group1": "ctrl", "group2": "treat"})
    print(f"  wrote {path}")


def write_params_yaml(path: Path, outdir: Path) -> None:
    abs_out = outdir.resolve()
    content = textwrap.dedent(f"""\
        # Synthetic SE test data — ctrl vs treat (3 replicates each).
        # MAJIQ is disabled: requires academic licence + user-built container.
        #
        # Launch:
        #   nextflow run main.nf -profile docker,test   -params-file {path}
        #   nextflow run main.nf -profile apptainer,test -params-file {path}
        #
        # To also run MAJIQ, add:
        #   --run_majiq true --majiq_sif <path/to/majiq.sif> --majiq_license <path/to/license>

        input:        {abs_out}/samplesheet.csv
        comparisons:  {abs_out}/comparisons.csv
        outdir:       {abs_out}/results

        gtf:          {abs_out}/ref/test.gtf
        genome_fasta: {abs_out}/ref/test.genome.fa
        use_gffread:  true

        strandedness: unstranded
        read_length:  {READ_LEN}

        run_rmats: true
        run_majiq: false
        run_isar:  true
    """)
    with open(path, "w") as fh:
        fh.write(content)
    print(f"  wrote {path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir).resolve()

    print(f"\nGenerating synthetic SE test data → {outdir}\n")

    ref_dir     = outdir / "ref"
    samples_dir = outdir / "samples"
    ref_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Genome FASTA ──
    print("[1/5] genome FASTA")
    genome = make_genome_seq(args.seed)
    write_fasta(ref_dir / "test.genome.fa", CHROM, genome)

    # ── 2. GTF ──
    print("[2/5] GTF annotation")
    write_gtf(ref_dir / "test.gtf")

    # ── 3. BAM + quant.sf per sample ──
    print("[3/5] BAM files + Salmon quant.sf")
    all_jct    = JUNCTION_COUNTS["ctrl"]  + JUNCTION_COUNTS["treat"]
    all_salmon = SALMON_COUNTS["ctrl"]    + SALMON_COUNTS["treat"]

    for (sample, condition, _replicate), jct, sal in zip(SAMPLES, all_jct, all_salmon):
        sample_dir = samples_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        n_e1e2, n_e2e3, n_e1e3 = jct
        write_bam(
            bam_path    = sample_dir / f"{sample}.bam",
            genome      = genome,
            n_e1e2      = n_e1e2,
            n_e2e3      = n_e2e3,
            n_e1e3      = n_e1e3,
            sample_name = sample,
        )

        n_tx1, n_tx2 = sal
        quant_path = sample_dir / "quant.sf"
        write_quant_sf(quant_path, n_tx1, n_tx2)
        n_bg = N_BG_GENES * 2
        print(f"  wrote {quant_path}  (tx1={n_tx1}, tx2={n_tx2}, +{n_bg} bg isoforms)")

    # ── 4. Samplesheet + comparisons ──
    print("[4/5] samplesheet and comparisons")
    write_samplesheet(outdir / "samplesheet.csv", outdir)
    write_comparisons(outdir / "comparisons.csv")

    # ── 5. params.yaml ──
    print("[5/5] params.yaml")
    write_params_yaml(outdir / "params.yaml", outdir)

    print(f"\nDone. Total samples: {len(SAMPLES)}")
    print(f"  ctrl  PSI target: 0.80  (IJC≈80, SJC≈10)")
    print(f"  treat PSI target: 0.20  (IJC≈20, SJC≈40)")
    print(f"\nLaunch:")
    print(f"  nextflow run main.nf -profile docker,test \\")
    print(f"    -params-file {outdir}/params.yaml\n")


if __name__ == "__main__":
    main()
