#!/usr/bin/env python3
"""
Finalise targeted sashimi plots:

  * rename each rmats2sashimiplot output  <rank>_<GENE>_<coords>.pdf
    to  A_<GENE>_<EVENT_TYPE>_<coords>.pdf  (type + coordinates keep two events
    of the same gene distinct, and the A_ prefix marks them as targeted),
  * write sashimi_out/targeted_classes.tsv mapping each final plot id back to
    its event metadata (gene, event_id, deltapsi, fdr, site_class).

Input: the filtered_events manifest written by filter_rmats_for_sashimi.py
(--gene-must mode) and the sashimi_out tree produced by rmats2sashimiplot.

Stdlib only.
"""

import argparse
import csv
import os
import re
import sys

EVENT_TYPES = ("SE", "A5SS", "A3SS", "MXE", "RI")
RANK_RE = re.compile(r"^(\d+)_(.*)$")


def load_manifest(path):
    """(event_type, rank) -> metadata dict."""
    meta = {}
    if not os.path.isfile(path):
        return meta
    with open(path, "r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                key = (row["event_type"], int(row["rank"]))
            except (KeyError, TypeError, ValueError):
                continue
            meta[key] = row
    return meta


def _read_tsv_rows(path):
    if not os.path.isfile(path):
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def plot_dir_for(root, event_type):
    """rmats2sashimiplot versions differ: <ETYPE>/Sashimi_plot/ or <ETYPE>/."""
    for candidate in (os.path.join(root, event_type, "Sashimi_plot"),
                      os.path.join(root, event_type)):
        if os.path.isdir(candidate) and any(
                n.endswith(".pdf") for n in os.listdir(candidate)):
            return candidate
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targeted-events", required=True,
                        help="dir with targeted_classes.tsv + <ETYPE>.top.txt")
    parser.add_argument("--sashimi-out", required=True,
                        help="sashimi_out tree produced by rmats2sashimiplot")
    args = parser.parse_args()

    meta = load_manifest(os.path.join(args.targeted_events, "targeted_classes.tsv"))
    root = args.sashimi_out
    rows = []
    n_renamed = 0

    for event_type in EVENT_TYPES:
        pdir = plot_dir_for(root, event_type)
        if not pdir:
            continue
        for name in sorted(os.listdir(pdir)):
            if not name.endswith(".pdf"):
                continue
            match = RANK_RE.match(name)
            if not match:
                continue
            rank = int(match.group(1))
            rest = match.group(2)  # <GENE>_<coords>
            record = meta.get((event_type, rank), {})
            gene = (record.get("gene") or "").strip()
            if not gene:
                gene = rest.split("_", 1)[0]
            coords = rest[len(gene) + 1:] if rest.startswith(gene + "_") else rest
            stem = "A_%s_%s_%s" % (gene, event_type, coords)
            final_name = stem + ".pdf"
            os.replace(os.path.join(pdir, name), os.path.join(pdir, final_name))
            n_renamed += 1
            rows.append({
                "event_type": event_type,
                "plot_id": stem,
                "gene": gene,
                "event_id": record.get("event_id", ""),
                "deltapsi": record.get("deltapsi", ""),
                "fdr": record.get("fdr", ""),
                "site_class": record.get("site_class", "annotated"),
                "note": "",
            })

    # Requested genes with no significant event are recorded too, so the manifest
    # documents them even though no PDF exists.
    for record in _read_tsv_rows(os.path.join(args.targeted_events, "targeted_classes.tsv")):
        if record.get("note") == "no_significant_events":
            rows.append({
                "event_type": "", "plot_id": "", "gene": (record.get("gene") or "").strip(),
                "event_id": "", "deltapsi": "", "fdr": "", "site_class": "",
                "note": "no_significant_events",
            })

    out_path = os.path.join(root, "targeted_classes.tsv")
    columns = ["event_type", "plot_id", "gene", "event_id",
               "deltapsi", "fdr", "site_class", "note"]
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print("[finalize_targeted_sashimi] renamed %d PDFs" % n_renamed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
