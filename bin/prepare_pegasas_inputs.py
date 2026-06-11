#!/usr/bin/env python3
"""
Prepare PEGASAS input matrices from pipeline outputs.

Sample ids for the PSI matrix come from the rMATS BAM list (b1.txt / b2.txt
order) and are passed in via --g1-ids / --g2-ids. They are validated against
the per-event column counts of SE.MATS.JC.txt to make sure the assignment is
unambiguous — if the counts do not match, the script aborts with a clear
error instead of silently producing wrong correlations.

Outputs:
  gene_exp_bySample.tsv  — rows=samples, cols=genes (TPM, samples ordered by group)
  PSI_bySample.tsv       — rows=SE events, cols=samples (IncLevel1/IncLevel2 merged)
  group_info.tsv         — sample_id<TAB>group (two lines per sample: sample<TAB>group)
  group_order.txt        — comma-separated group order for PEGASAS heatmap
"""

import argparse
import csv
import sys
import os


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("salmon_tpm",     help="salmon.merged.gene_tpm.tsv (genes × samples)")
    p.add_argument("rmats_se",       help="SE.MATS.JC.txt from rMATS output")
    p.add_argument("group_info_in",  help="TSV: sample_id<TAB>group (one row per sample)")
    p.add_argument("--g1-ids",       dest="g1_ids", default="",
                   help="Comma-separated sample_ids in the order rMATS POST wrote b1.txt")
    p.add_argument("--g2-ids",       dest="g2_ids", default="",
                   help="Comma-separated sample_ids in the order rMATS POST wrote b2.txt")
    p.add_argument("--out-dir",      default=".", dest="out_dir",
                   help="Output directory [.]")
    p.add_argument("--min-samples",  type=int, default=3, dest="min_samples",
                   help="Min samples with valid PSI per event [3]")
    return p.parse_args()


def parse_id_list(raw: str) -> list:
    return [s.strip() for s in raw.split(",") if s.strip()]


def load_group_info(fin: str) -> dict:
    """Return {sample_id: group}."""
    mapping = {}
    with open(fin) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            mapping[parts[0]] = parts[1]
    return mapping


def transpose_tpm(fin: str, sample_order: list, out_path: str) -> None:
    """
    Convert salmon.merged.gene_tpm.tsv (genes × samples) to
    PEGASAS format (samples × genes, first col = sample_id).
    """
    with open(fin) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Build gene list and expression dict: {gene_name: {sample: tpm}}
        genes = []
        expr  = {}
        for row in reader:
            gene = row.get("gene_name") or row.get("gene_id")
            if not gene:
                continue
            genes.append(gene)
            expr[gene] = {}
            for k, v in row.items():
                if k not in ("gene_id", "gene_name"):
                    try:
                        expr[gene][k] = float(v)
                    except (ValueError, TypeError):
                        expr[gene][k] = 0.0

    # Identify samples present in both TPM file and group info
    all_tpm_samples = set()
    if genes:
        all_tpm_samples = set(expr[genes[0]].keys())
    valid_samples = [s for s in sample_order if s in all_tpm_samples]

    if not valid_samples:
        sys.exit("[ERROR] No samples matched between TPM file and group info")

    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SampleID"] + genes)
        for sample in valid_samples:
            row = [sample] + [expr[g].get(sample, 0.0) for g in genes]
            writer.writerow(row)

    print(f"[INFO] Gene expression matrix: {len(valid_samples)} samples × {len(genes)} genes → {out_path}")


def build_psi_matrix(fin: str, g1_ids: list, g2_ids: list, out_path: str,
                     min_samples: int) -> None:
    """
    Parse SE.MATS.JC.txt and output PSI matrix (events × samples).

    Sample ids are taken from --g1-ids / --g2-ids (the order rMATS POST used
    for b1.txt / b2.txt) and validated against the per-event IncLevel1 /
    IncLevel2 column counts.
    """
    with open(fin) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    if not rows:
        sys.exit("[ERROR] SE.MATS.JC.txt is empty")

    first = rows[0]
    n1 = len(first["IncLevel1"].split(","))
    n2 = len(first["IncLevel2"].split(","))

    if len(g1_ids) != n1:
        sys.exit(
            f"[ERROR] SE.MATS.JC.txt has {n1} samples in IncLevel1 but "
            f"--g1-ids declares {len(g1_ids)} ({g1_ids}). The BAM list "
            f"order used by rMATS POST is the source of truth; this mismatch "
            f"would silently misassign PSI to the wrong sample."
        )
    if len(g2_ids) != n2:
        sys.exit(
            f"[ERROR] SE.MATS.JC.txt has {n2} samples in IncLevel2 but "
            f"--g2-ids declares {len(g2_ids)} ({g2_ids}). The BAM list "
            f"order used by rMATS POST is the source of truth; this mismatch "
            f"would silently misassign PSI to the wrong sample."
        )

    all_samples = g1_ids + g2_ids

    header_cols = ["AC", "GeneName", "chr", "strand",
                   "exonStart", "exonEnd", "upstreamEE", "downstreamES"]

    events_written = 0
    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header_cols + all_samples)
        for row in rows:
            psi_vals1 = row["IncLevel1"].split(",")
            psi_vals2 = row["IncLevel2"].split(",")
            psi_all   = psi_vals1 + psi_vals2

            # Defensive: should never trigger after the validation above, but
            # different rows could in theory have variable N if rMATS ever
            # emitted them.
            if len(psi_all) != len(all_samples):
                sys.exit(
                    f"[ERROR] Row {row.get('ID', '?')} has "
                    f"{len(psi_all)} PSI values but expected "
                    f"{len(all_samples)} (g1={n1}, g2={n2}). Aborting to "
                    f"avoid silent misassignment."
                )

            valid = sum(1 for v in psi_all if v not in ("", "NA", "na"))
            if valid < min_samples:
                continue

            psi_clean = ["NA" if v in ("", "NA", "na") else v for v in psi_all]

            event_id = "{}_{}_{}_{}_{}_{}_{}_{}".format(
                row.get("ID", ""),
                row.get("GeneID", ""),
                row.get("chr", ""),
                row.get("strand", ""),
                row.get("exonStart_0base", row.get("exonStart", "")),
                row.get("exonEnd", ""),
                row.get("upstreamEE", ""),
                row.get("downstreamES", ""),
            )

            meta = [
                event_id,
                row.get("GeneID", ""),
                row.get("chr", ""),
                row.get("strand", ""),
                row.get("exonStart_0base", row.get("exonStart", "")),
                row.get("exonEnd", ""),
                row.get("upstreamEE", ""),
                row.get("downstreamES", ""),
            ]
            writer.writerow(meta + psi_clean)
            events_written += 1

    print(f"[INFO] PSI matrix: {events_written} events × {len(all_samples)} samples → {out_path}")


def write_group_info(sample_order: list, group_map: dict, out_path: str) -> None:
    with open(out_path, "w") as fh:
        for sample in sample_order:
            group = group_map.get(sample, "unknown")
            fh.write(f"{sample}\t{group}\n")
    print(f"[INFO] Group info written → {out_path}")


def write_group_order(group_map: dict, out_path: str) -> None:
    """Write unique groups in the order they first appear."""
    seen   = []
    unique = []
    for g in group_map.values():
        if g not in seen:
            seen.append(g)
            unique.append(g)
    with open(out_path, "w") as fh:
        fh.write(",".join(unique) + "\n")
    print(f"[INFO] Group order: {unique} → {out_path}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    group_map    = load_group_info(args.group_info_in)
    g1_ids       = parse_id_list(args.g1_ids)
    g2_ids       = parse_id_list(args.g2_ids)

    if not g1_ids or not g2_ids:
        sys.exit(
            "[ERROR] --g1-ids and --g2-ids are required. They are propagated "
            "from the rMATS BAM list order so that PSI values are assigned "
            "to the correct sample in the matrix."
        )

    # Every rMATS sample must be present in group_info — otherwise the
    # group_info used downstream would miss samples and the correlation would
    # silently skip them.
    missing = [s for s in (g1_ids + g2_ids) if s not in group_map]
    if missing:
        sys.exit(
            f"[ERROR] {len(missing)} sample(s) from rMATS are missing in "
            f"group_info: {missing}. Add them to the group_info TSV."
        )

    sample_order = list(group_map.keys())

    tpm_out = os.path.join(args.out_dir, "gene_exp_bySample.tsv")
    transpose_tpm(args.salmon_tpm, sample_order, tpm_out)

    psi_out = os.path.join(args.out_dir, "PSI_bySample.tsv")
    build_psi_matrix(args.rmats_se, g1_ids, g2_ids, psi_out, args.min_samples)

    grp_out = os.path.join(args.out_dir, "group_info.tsv")
    write_group_info(sample_order, group_map, grp_out)

    ord_out = os.path.join(args.out_dir, "group_order.txt")
    write_group_order(group_map, ord_out)


if __name__ == "__main__":
    main()
