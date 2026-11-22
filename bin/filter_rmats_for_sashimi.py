#!/usr/bin/env python3
"""
Filter rMATS output for sashimi plot generation.

Two modes:

* Default (top-N): selects the top N events per event type ranked by a combined
  priority score  priority = -log10(FDR) * |IncLevelDifference|.
  Only events passing the FDR and |dPSI| cutoffs are considered. Writes one
  filtered file per event type plus site_classes.tsv.

* Targeted (--gene-must FILE): selects every significant event whose geneSymbol
  is in the provided list. The list may be one symbol per line or a CSV with a
  gene_symbol/gene column; blank lines and '#' are ignored and matching is
  case-insensitive. No top-N cap. Writes one filtered file per event type plus
  targeted_classes.tsv. Genes in the list without any significant event in this
  contrast are recorded with note == "no_significant_events".
"""

import argparse
import json
import math
import os
import sys


EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]


def priority_score(fdr: float, dpsi: float) -> float:
    """Higher is more significant + larger effect."""
    fdr_clamped = max(fdr, 1e-300)
    return -math.log10(fdr_clamped) * abs(dpsi)


def read_gene_must(path: str) -> "set[str]":
    """Read a gene list from a plain one-symbol-per-line file or a CSV.

    Accepted layouts: one symbol per line, a single-column CSV, or a CSV with a
    'gene_symbol'/'gene' column (other columns ignored). Blank lines and '#' are
    skipped; matching is case-insensitive.
    """
    genes = set()
    with open(path) as fh:
        lines = [ln.strip() for ln in fh]
    lines = [ln for ln in lines if ln and not ln.startswith("#")]
    if not lines:
        return genes

    # Detect a CSV: any line with a comma or tab delimiter.
    csv_like = any("," in ln or "\t" in ln for ln in lines)
    if csv_like:
        import csv
        import io
        delimiter = "," if any("," in ln for ln in lines) else "\t"
        reader = csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter)
        rows = [r for r in reader if r]
        if not rows:
            return genes
        header = [c.strip().lower() for c in rows[0]]
        col = 0
        for name in ("gene_symbol", "gene", "symbol"):
            if name in header:
                col = header.index(name)
                break
        else:
            # No recognised header: treat the first row as data.
            for row in rows:
                if row and row[0].strip():
                    genes.add(row[0].strip().upper())
            return genes
        for row in rows[1:]:
            if len(row) > col and row[col].strip():
                genes.add(row[col].strip().upper())
        return genes

    # Single-column CSV: drop a recognised header row if present (keep a
    # single-gene list intact even if the symbol looks like a header).
    if len(lines) > 1 and lines[0].strip('"').lower() in ("gene_symbol", "gene", "symbol"):
        lines = lines[1:]
    for line in lines:
        genes.add(line.strip('"').upper())
    return genes


def load_novel_ids(rmats_dir: str, event_type: str) -> "set[str]":
    """De novo (unannotated) splice sites: match IDs against fromGTF.novelSpliceSite.

    rMATS 4.3 no longer encodes them as negative coordinates, so the event ID is
    matched against the novel-splice-site file instead.
    """
    novel_ids = set()
    nf = os.path.join(rmats_dir, f"fromGTF.novelSpliceSite.{event_type}.txt")
    if os.path.isfile(nf):
        with open(nf) as fh:
            fh.readline()  # header
            for ln in fh:
                novel_ids.add(ln.rstrip("\n").split("\t")[0])
    return novel_ids


def filter_events(rmats_dir: str, out_dir: str, top_n: int,
                  fdr_cutoff: float, dpsi_cutoff: float) -> None:
    os.makedirs(out_dir, exist_ok=True)
    selected = {}
    # (event_type, rank, event_id, site_class) for every plotted event. rank is
    # the 1-based row order in the .top.txt file, which is exactly the <rank>
    # prefix rmats2sashimiplot writes into the PDF filename, so downstream
    # consumers can recover each PDF's splice-site class.
    manifest_rows = []

    for event_type in EVENT_TYPES:
        jc_file = os.path.join(rmats_dir, f"{event_type}.MATS.JC.txt")
        if not os.path.isfile(jc_file):
            continue

        novel_ids = load_novel_ids(rmats_dir, event_type)

        with open(jc_file) as fh:
            header = fh.readline().rstrip("\n").split("\t")

        try:
            fdr_idx  = header.index("FDR")
            dpsi_idx = header.index("IncLevelDifference")
        except ValueError:
            print(f"[WARN] {event_type}: missing FDR or IncLevelDifference column — skipping",
                  file=sys.stderr)
            continue

        candidates = []
        n_novel = 0
        with open(jc_file) as fh:
            fh.readline()  # skip header
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(fdr_idx, dpsi_idx):
                    continue
                try:
                    fdr_val  = float(parts[fdr_idx])
                    dpsi_val = float(parts[dpsi_idx])
                except ValueError:
                    continue

                # novelSS events are kept: rMATS 4.3 no longer encodes them as
                # negative coordinates, so rmats2sashimiplot can draw them.
                if parts[0] in novel_ids:
                    n_novel += 1

                if fdr_val <= fdr_cutoff and abs(dpsi_val) >= dpsi_cutoff:
                    score = priority_score(fdr_val, dpsi_val)
                    candidates.append((score, parts[0], line))

        if n_novel > 0:
            print(f"[INFO] {event_type}: included {n_novel} novel-splice-site events",
                  file=sys.stderr)

        if not candidates:
            continue

        # Sort descending by priority score; take top N
        candidates.sort(key=lambda x: x[0], reverse=True)
        top = candidates[:top_n]

        out_file = os.path.join(out_dir, f"{event_type}.top.txt")
        with open(out_file, "w") as fh:
            # Write header as-is from original file
            with open(jc_file) as src:
                fh.write(src.readline())
            for rank, (_, event_id, line) in enumerate(top, start=1):
                fh.write(line)
                manifest_rows.append((
                    event_type,
                    rank,
                    event_id,
                    "de_novo" if event_id in novel_ids else "annotated",
                ))

        selected[event_type] = len(top)
        print(f"[INFO] {event_type}: selected {len(top)} / {len(candidates)} significant events")

    if not selected:
        print("[WARN] No significant events found for any event type. "
              "Check FDR/dPSI cutoffs.", file=sys.stderr)
        sys.exit(0)

    # Note that novel-splice-site events are kept and plotted: rmats2sashimiplot
    # can draw them from rMATS 4.3 coordinates, so annotation is a label, not a
    # filter. Every selected event is recorded below so the report and the PDF
    # index never mislabel a de novo plot as annotated.
    with open(os.path.join(out_dir, "site_classes.tsv"), "w") as fh:
        fh.write("event_type\trank\tevent_id\tsite_class\n")
        for event_type, rank, event_id, site_class in manifest_rows:
            fh.write(f"{event_type}\t{rank}\t{event_id}\t{site_class}\n")

    # Write a summary JSON for the report
    summary = {"top_n": top_n, "fdr_cutoff": fdr_cutoff, "dpsi_cutoff": dpsi_cutoff,
               "selected": selected}
    with open(os.path.join(out_dir, "filter_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


def select_targeted(rmats_dir: str, out_dir: str, gene_must: "set[str]",
                    fdr_cutoff: float, dpsi_cutoff: float) -> None:
    """Every significant event whose geneSymbol is in gene_must, per event type."""
    os.makedirs(out_dir, exist_ok=True)
    selected = {}
    manifest_rows = []
    matched_genes = set()

    for event_type in EVENT_TYPES:
        jc_file = os.path.join(rmats_dir, f"{event_type}.MATS.JC.txt")
        if not os.path.isfile(jc_file):
            continue

        novel_ids = load_novel_ids(rmats_dir, event_type)

        with open(jc_file) as fh:
            header = fh.readline().rstrip("\n").split("\t")

        try:
            fdr_idx  = header.index("FDR")
            dpsi_idx = header.index("IncLevelDifference")
            gene_idx = header.index("geneSymbol")
        except ValueError:
            print(f"[WARN] {event_type}: missing FDR/IncLevelDifference/geneSymbol "
                  f"column — skipping", file=sys.stderr)
            continue

        kept = []
        with open(jc_file) as fh:
            fh.readline()  # skip header
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(fdr_idx, dpsi_idx, gene_idx):
                    continue
                gene = parts[gene_idx].strip()
                if not gene or gene.upper() not in gene_must:
                    continue
                try:
                    fdr_val  = float(parts[fdr_idx])
                    dpsi_val = float(parts[dpsi_idx])
                except ValueError:
                    continue
                if fdr_val <= fdr_cutoff and abs(dpsi_val) >= dpsi_cutoff:
                    kept.append((parts[0], gene, fdr_val, dpsi_val, line))

        if not kept:
            continue

        out_file = os.path.join(out_dir, f"{event_type}.top.txt")
        with open(out_file, "w") as fh:
            with open(jc_file) as src:
                fh.write(src.readline())
            for rank, (event_id, gene, fdr_val, dpsi_val, line) in enumerate(kept, start=1):
                fh.write(line)
                site_class = "de_novo" if event_id in novel_ids else "annotated"
                manifest_rows.append((
                    event_type, rank, gene, event_id, dpsi_val, fdr_val, site_class))
                matched_genes.add(gene.upper())

        selected[event_type] = len(kept)
        print(f"[INFO] targeted {event_type}: selected {len(kept)} significant events")

    # Genes requested but with no significant event in this contrast.
    missing = sorted(gene_must - matched_genes)
    with open(os.path.join(out_dir, "targeted_classes.tsv"), "w") as fh:
        fh.write("event_type\trank\tgene\tevent_id\tdeltapsi\tfdr\tsite_class\tnote\n")
        for event_type, rank, gene, event_id, dpsi_val, fdr_val, site_class in manifest_rows:
            fh.write(f"{event_type}\t{rank}\t{gene}\t{event_id}\t"
                     f"{dpsi_val:.6g}\t{fdr_val:.6g}\t{site_class}\t\n")
        for gene in missing:
            fh.write(f"\t\t{gene}\t\t\t\t\tno_significant_events\n")

    summary = {
        "mode": "targeted",
        "fdr_cutoff": fdr_cutoff,
        "dpsi_cutoff": dpsi_cutoff,
        "requested_genes": len(gene_must),
        "genes_with_events": len(matched_genes),
        "genes_without_events": len(missing),
        "selected": selected,
    }
    with open(os.path.join(out_dir, "filter_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    if not selected:
        print("[WARN] No significant events found for any requested gene.",
              file=sys.stderr)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rmats_dir",   help="Directory with rMATS output files")
    p.add_argument("out_dir",     help="Output directory for filtered event files")
    p.add_argument("--top-n",     type=int,   default=10,   help="Events per event type [10]")
    p.add_argument("--fdr",       type=float, default=0.05, help="FDR cutoff [0.05]")
    p.add_argument("--dpsi",      type=float, default=0.1,  help="|dPSI| cutoff [0.1]")
    p.add_argument("--gene-must", default=None,
                   help="Targeted mode: path to a one-symbol-per-line gene list")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.gene_must:
        select_targeted(
            rmats_dir   = args.rmats_dir,
            out_dir     = args.out_dir,
            gene_must   = read_gene_must(args.gene_must),
            fdr_cutoff  = args.fdr,
            dpsi_cutoff = args.dpsi,
        )
    else:
        filter_events(
            rmats_dir   = args.rmats_dir,
            out_dir     = args.out_dir,
            top_n       = args.top_n,
            fdr_cutoff  = args.fdr,
            dpsi_cutoff = args.dpsi,
        )
