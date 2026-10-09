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
import math
import os
import shutil
import sys


def bh_adjust(pvalues: list) -> list:
    """Benjamini-Hochberg adjusted p-values, order-preserving."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [0.0] * n
    previous = 1.0
    for position in range(n - 1, -1, -1):
        index = order[position]
        value = pvalues[index] * n / float(position + 1)
        previous = min(previous, value)
        adjusted[index] = min(previous, 1.0)
    return adjusted


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the regularized incomplete beta (Lentz)."""
    max_it, eps, fpmin = 200, 3e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_it + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_beta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) +
                a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(log_beta)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def pearson_two_sided_p(r: float, n_samples: int):
    """Two-sided p-value for a Pearson r using the t distribution (df = n-2).

    The permutation p-value is floored at 1/(N_PERMS+1), so BH over the full
    event x pathway matrix cannot reject anything at ~1M tests. The analytic
    p is continuous and gives BH real resolution; the permutation p is retained
    as an independent cross-check.
    """
    df = n_samples - 2
    if df <= 0 or r is None:
        return None
    if r >= 1.0:
        r = 1.0 - 1e-12
    if r <= -1.0:
        r = -1.0 + 1e-12
    t_stat = r * math.sqrt(df / (1.0 - r * r))
    x = df / (df + t_stat * t_stat)
    return _betai(df / 2.0, 0.5, x)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("shared_scores_dir", help="Shared PEGASAS pathway_out/ (contains <pathway>.scores.txt)")
    p.add_argument("correlation_out",    help="Per-contrast correlation_out/ directory")
    p.add_argument("--g1-ids",  dest="g1_ids", default="", help="Comma-separated contrast group1 sample ids")
    p.add_argument("--g2-ids",  dest="g2_ids", default="", help="Comma-separated contrast group2 sample ids")
    p.add_argument("--comp-id", required=True, help="Comparison id (used to name the sig summary)")
    p.add_argument("--fdr-cutoff", dest="fdr_cutoff", type=float, default=0.05,
                   help="BH FDR applied across the contrast event x pathway matrix [0.05]")
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
    has_background = False
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
                n_background = parts[5] if len(parts) > 5 else ""
                n_hits = parts[6] if len(parts) > 6 else ""
                if n_background != "" or n_hits != "":
                    has_background = True
                rows.append([pathway, sample, group, ks, pv, med, n_background, n_hits])

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        header = ["pathway", "sample", "group", "KS_score", "p_value", "median_rank"]
        if has_background:
            header += ["n_background", "n_hits"]
        writer.writerow(header)
        for row in rows:
            writer.writerow(row[:len(header)])
    print(f"[INFO] all_pathways_scores.tsv: {len(rows)} rows across {len(found)} pathways → {out_path}")
    return found


def count_sig_events(correlation_dir: str, pathways: list, out_path: str,
                     comp_id: str, fdr_cutoff: float, n_samples: int) -> None:
    """
    Count BH-significant event x pathway correlations per pathway.

    The analytic Pearson p-value (t distribution) is pooled across ALL
    pathway x event tests of the contrast and BH-adjusted once; n_sig_events is
    the count with padj < FDR. The permutation p is also BH-corrected and
    reported as n_sig_events_perm purely as a cross-check: its 1/N_PERMS floor
    makes it unable to reject at ~1M tests.

    Writes pathway, n_sig_events, n_sig_events_perm, n_tested, min_padj.
    """
    analytic = []        # (pathway, p_analytic)
    perm = []            # (pathway, p_permutation)
    n_tested = {}        # pathway -> rows with a usable analytic p
    for pathway in pathways:
        global_file = os.path.join(
            correlation_dir, pathway, f"{pathway}_global_cor_matrix.txt")
        n_tested[pathway] = 0
        if not os.path.isfile(global_file):
            continue
        with open(global_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    r = float(row.get("pearson_r"))
                except (TypeError, ValueError):
                    continue
                p_analytic = pearson_two_sided_p(r, n_samples)
                if p_analytic is None:
                    continue
                analytic.append((pathway, p_analytic))
                n_tested[pathway] += 1
                try:
                    perm.append((pathway, float(row.get("p_value"))))
                except (TypeError, ValueError):
                    pass

    adj_analytic = bh_adjust([p for _pathway, p in analytic])
    adj_perm = bh_adjust([p for _pathway, p in perm])

    counts = {pathway: 0 for pathway in pathways}
    counts_perm = {pathway: 0 for pathway in pathways}
    min_padj = {}
    for (pathway, _p), q in zip(analytic, adj_analytic):
        if pathway not in min_padj or q < min_padj[pathway]:
            min_padj[pathway] = q
        if q < fdr_cutoff:
            counts[pathway] += 1
    for (pathway, _p), q in zip(perm, adj_perm):
        if q < fdr_cutoff:
            counts_perm[pathway] += 1

    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["pathway", "n_sig_events", "n_sig_events_perm",
                         "n_tested", "min_padj"])
        for pathway in pathways:
            writer.writerow([
                pathway,
                counts.get(pathway, 0),
                counts_perm.get(pathway, 0),
                n_tested.get(pathway, 0),
                "" if pathway not in min_padj else "%.6g" % min_padj[pathway],
            ])
    print(f"[INFO] {comp_id}_sig_pathways.tsv: {sum(counts.values())} BH-significant events "
          f"(analytic, pooled {len(analytic)} tests, FDR < {fdr_cutoff}); "
          f"permutation cross-check: {sum(counts_perm.values())} → {out_path}")


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
    count_sig_events(args.correlation_out, pathways, sig_out, args.comp_id,
                     args.fdr_cutoff, len(contrast_samples))

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
