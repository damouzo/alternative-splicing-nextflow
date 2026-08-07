#!/usr/bin/env python3
"""
PEGASAS collation for a single contrast.

Gathers the per-contrast PEGASAS result into one directory and produces a
single long-format table of pathway KS scores for the report plots:

  all_pathways_scores.tsv — pathway, sample, group, KS_score, p_value, median_rank
                            (only the contrast samples, from the shared scores)

Also summarises, per pathway, how many splicing events were significantly
correlated (from the correlation step's *_high_cor_matrix.txt) so that a
cross-contrast upset/heatmap can be built:

  <comp_id>_sig_pathways.tsv — pathway, n_sig_events

The correlation output is copied alongside so the per-contrast report
receives a single pegasas directory.
"""

import argparse
import csv
import glob
import os
import shutil
import sys


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("shared_scores_dir", help="Shared PEGASAS pathway_out/ (contains <pathway>.scores.txt)")
    p.add_argument("correlation_out",    help="Per-contrast correlation_out/ directory")
    p.add_argument("--g1-ids",  dest="g1_ids", default="", help="Comma-separated contrast group1 sample ids")
    p.add_argument("--g2-ids",  dest="g2_ids", default="", help="Comma-separated contrast group2 sample ids")
    p.add_argument("--comp-id", required=True, help="Comparison id (used to name the sig summary)")
    p.add_argument("--out-dir", default=".", dest="out_dir", help="Output directory for all_pathways_scores.tsv + correlation_out [.]")
    p.add_argument("--sig-out", default=None, dest="sig_out",
                   help="Where to write <comp_id>_sig_pathways.tsv [out_dir]; in the pipeline it is set to the work root so it can be emitted separately from pegasas_out/")
    return p.parse_args()


def parse_id_list(raw: str) -> list:
    return [s.strip() for s in raw.split(",") if s.strip()]


def collect_scores(shared_dir: str, contrast_samples: set, out_path: str) -> list:
    """
    Read every <pathway>.scores.txt under shared_dir (sample<TAB>group<TAB>
    KS<TAB>p<TAB>median_rank, no header), subset to the contrast samples,
    and write a single long table with a pathway column.

    PEGASAS writes scores at shared_dir/<pathway>/<pathway>.scores.txt (one
    level of nesting), so glob recursively. shared_dir may be a symlink to the
    shared pathway run; glob/Python follow symlinks transparently.
    """
    rows = []
    found = []
    # Match shared_dir/<pathway>/<pathway>.scores.txt (and any deeper nesting).
    score_files = sorted(glob.glob(os.path.join(shared_dir, "*", "*.scores.txt")))
    for sfile in score_files:
        entry = os.path.basename(sfile)
        if not entry.endswith(".scores.txt"):
            continue
        pathway = entry[: -len(".scores.txt")]
        found.append(pathway)
        with open(sfile) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                sample, group, ks, pv, med = parts[0], parts[1], parts[2], parts[3], parts[4]
                if sample not in contrast_samples:
                    continue
                rows.append([pathway, sample, group, ks, pv, med])

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["pathway", "sample", "group", "KS_score", "p_value", "median_rank"])
        writer.writerows(rows)
    print(f"[INFO] all_pathways_scores.tsv: {len(rows)} rows across {len(found)} pathways → {out_path}")
    return found


def count_sig_events(correlation_dir: str, pathways: list, out_path: str, comp_id: str) -> None:
    """
    For each pathway, count rows in <pathway>_high_cor_matrix.txt (significant
    events only). Writes pathway<TAB>n_sig_events.
    """
    counts = []
    for pathway in pathways:
        # correlation_out/<pathway>/<pathway>_high_cor_matrix.txt
        high_file = os.path.join(correlation_dir, pathway, f"{pathway}_high_cor_matrix.txt")
        n = 0
        if os.path.isfile(high_file):
            with open(high_file) as fh:
                for i, line in enumerate(fh):
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    # Skip a header row (e.g. second field "pearson_r")
                    if i == 0 and len(parts) >= 2:
                        try:
                            float(parts[1])
                        except ValueError:
                            continue
                    if len(parts) >= 1:
                        n += 1
        counts.append([pathway, n])

    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["pathway", "n_sig_events"])
        writer.writerows(counts)
    print(f"[INFO] {comp_id}_sig_pathways.tsv: {sum(c for _, c in counts)} significant events "
          f"across {len(pathways)} pathways → {out_path}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    contrast_samples = set(parse_id_list(args.g1_ids) + parse_id_list(args.g2_ids))
    if not contrast_samples:
        sys.exit("[ERROR] --g1-ids and --g2-ids are required for collation")

    # all_pathways_scores.tsv (long, contrast subset)
    scores_out = os.path.join(args.out_dir, "all_pathways_scores.tsv")
    pathways = collect_scores(args.shared_scores_dir, contrast_samples, scores_out)

    # significant-event summary for the cross-contrast upset
    sig_dir = args.sig_out if args.sig_out else args.out_dir
    os.makedirs(sig_dir, exist_ok=True)
    sig_out = os.path.join(sig_dir, f"{args.comp_id}_sig_pathways.tsv")
    count_sig_events(args.correlation_out, pathways, sig_out, args.comp_id)

    # gather correlation output into the same directory
    dest_corr = os.path.join(args.out_dir, "correlation_out")
    if os.path.abspath(args.correlation_out) == os.path.abspath(dest_corr):
        pass  # already in place
    elif os.path.isdir(args.correlation_out):
        if os.path.exists(dest_corr):
            shutil.rmtree(dest_corr)
        shutil.copytree(args.correlation_out, dest_corr)
        print(f"[INFO] correlation_out gathered → {dest_corr}")


if __name__ == "__main__":
    main()
