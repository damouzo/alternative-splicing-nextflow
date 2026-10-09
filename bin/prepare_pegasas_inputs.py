#!/usr/bin/env python3
"""
Prepare PEGASAS input matrices from pipeline outputs.

Sample ids for the PSI matrix come from the rMATS BAM list (b1.txt / b2.txt
order) and are passed in via --g1-ids / --g2-ids. They are validated against
the per-event column counts of SE.MATS.JC.txt to make sure the assignment is
unambiguous — if the counts do not match, the script aborts with a clear
error instead of silently producing wrong correlations. When the materialized
b1_samples.txt / b2_samples.txt written by rMATS POST are provided, their
order is also compared positionally against --g1-ids / --g2-ids.

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
    p.add_argument("--salmon-tpm", dest="salmon_tpm", default=None,
                   help="salmon.merged.gene_tpm.tsv (optional)")
    p.add_argument("rmats_se",       help="SE.MATS.JC.txt from rMATS output")
    p.add_argument("group_info_in",  help="TSV: sample_id<TAB>group (one row per sample)")
    p.add_argument("--g1-ids",       dest="g1_ids", default="",
                   help="Comma-separated sample_ids in the order rMATS POST wrote b1.txt")
    p.add_argument("--g2-ids",       dest="g2_ids", default="",
                   help="Comma-separated sample_ids in the order rMATS POST wrote b2.txt")
    p.add_argument("--g1-samples-file", dest="g1_samples_file", default=None,
                   help="b1_samples.txt written by rMATS POST: materialized "
                        "column->sample order for group 1 (validated positionally "
                        "against --g1-ids)")
    p.add_argument("--g2-samples-file", dest="g2_samples_file", default=None,
                   help="b2_samples.txt written by rMATS POST: materialized "
                        "column->sample order for group 2 (validated positionally "
                        "against --g2-ids)")
    p.add_argument("--out-dir",      default=".", dest="out_dir",
                   help="Output directory [.]")
    p.add_argument("--fdr-cutoff",   type=float, default=0.0, dest="fdr_cutoff",
                   help="Skip events with FDR above this threshold (0 = no filter). "
                        "Set to e.g. 0.05 to keep only significant events [0.0]")
    p.add_argument("--dpsi-cutoff",  type=float, default=0.0, dest="dpsi_cutoff",
                   help="Skip events with |IncLevelDifference| below this threshold "
                        "(0 = no filter). Set with --fdr-cutoff to mirror the "
                        "rMATS Significant flag (FDR < x AND |dPSI| >= y) [0.0]")
    p.add_argument("--min-samples",  type=int, default=3, dest="min_samples",
                   help="Min samples with valid PSI per event [3]")
    return p.parse_args()


def parse_id_list(raw: str) -> list:
    return [s.strip() for s in raw.split(",") if s.strip()]


def validate_order(ids: list, ids_file: str, group_label: str) -> None:
    """
    Identity-of-order guard: when rMATS POST materialized the column->sample
    mapping (bN_samples.txt, one id per line in bN.txt order), verify it
    matches the channel-provided ids positionally. A cardinality match alone
    would let a crossed order pass silently.
    """
    if not ids_file:
        return
    if not os.path.isfile(ids_file):
        print(f"[WARN] {group_label}: {ids_file} not found — "
              "falling back to cardinality-only validation")
        return
    with open(ids_file) as fh:
        reference = [ln.strip() for ln in fh if ln.strip()]
    if reference == ids:
        print(f"[INFO] {group_label}: sample order matches the materialized bN_samples.txt")
        return
    sys.exit(
        f"[ERROR] {group_label} sample order mismatch: --g1-ids/--g2-ids order "
        f"{ids} does not match the materialized order {reference}. PSI columns "
        f"would be assigned to the wrong sample."
    )


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

    Duplicated gene symbols (several gene_ids sharing one symbol) are collapsed
    by keeping the row with the highest mean TPM across samples.

    Keep the collapse rule in sync with
    prepare_pegasas_gene_matrix.py.transpose_tpm (the two scripts share no
    import path, so the logic is mirrored).
    """
    with open(fin) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        order = []
        chosen = {}   # symbol -> (mean_tpm, {sample: tpm})
        n_duplicates = 0
        for row in reader:
            gene = row.get("gene_name") or row.get("gene_id")
            if not gene:
                continue
            values = {}
            for k, v in row.items():
                if k not in ("gene_id", "gene_name"):
                    try:
                        values[k] = float(v)
                    except (ValueError, TypeError):
                        values[k] = 0.0
            mean_tpm = (sum(values.values()) / len(values)) if values else 0.0
            if gene not in chosen:
                order.append(gene)
                chosen[gene] = (mean_tpm, values)
            else:
                n_duplicates += 1
                if mean_tpm > chosen[gene][0]:
                    chosen[gene] = (mean_tpm, values)

    genes = order
    expr = {g: chosen[g][1] for g in genes}
    if n_duplicates:
        print(f"[INFO] collapsed {n_duplicates} duplicate gene symbol row(s) "
              "by highest mean TPM")

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


def build_psi_matrix(fin, g1_ids, g2_ids, out_path,
                     min_samples, fdr_cutoff, dpsi_cutoff):
    """
    Parse SE.MATS.JC.txt and output PSI matrix (events × samples).
    When fdr_cutoff > 0, events with FDR >= cutoff (or missing/NA FDR) are
    skipped to keep only statistically significant splicing events.
    When dpsi_cutoff > 0, events with |IncLevelDifference| < cutoff (or
    missing/NA) are additionally skipped, mirroring the rMATS Significant flag.

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

    total_events  = 0
    skipped_fdr   = 0
    skipped_dpsi  = 0
    events_written = 0
    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header_cols + all_samples)
        for row in rows:
            total_events += 1

            if fdr_cutoff > 0:
                fdr_str = row.get("FDR", "")
                try:
                    fdr = float(fdr_str) if fdr_str not in ("", "NA", "na", "NaN") else None
                except ValueError:
                    fdr = None
                if fdr is None or fdr >= fdr_cutoff:
                    skipped_fdr += 1
                    continue
            if dpsi_cutoff > 0:
                dpsi_str = row.get("IncLevelDifference", "")
                try:
                    dpsi = float(dpsi_str) if dpsi_str not in ("", "NA", "na", "NaN") else None
                except ValueError:
                    dpsi = None
                if dpsi is None or abs(dpsi) < dpsi_cutoff:
                    skipped_dpsi += 1
                    continue
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
                row.get("geneSymbol", row.get("GeneID", "")),
                row.get("chr", ""),
                row.get("strand", ""),
                row.get("exonStart_0base", row.get("exonStart", "")),
                row.get("exonEnd", ""),
                row.get("upstreamEE", ""),
                row.get("downstreamES", ""),
            ]
            writer.writerow(meta + psi_clean)
            events_written += 1

    if fdr_cutoff > 0 or dpsi_cutoff > 0:
        filters = []
        if fdr_cutoff > 0:
            filters.append(f"{skipped_fdr}/{total_events} by FDR >= {fdr_cutoff}")
        if dpsi_cutoff > 0:
            filters.append(f"{skipped_dpsi}/{total_events} by |dPSI| < {dpsi_cutoff}")
        print(f"[INFO] PSI matrix: {events_written} events × {len(all_samples)} samples "
              f"(removed: {'; '.join(filters)}) → {out_path}")
    else:
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

    validate_order(g1_ids, args.g1_samples_file, "Group1")
    validate_order(g2_ids, args.g2_samples_file, "Group2")

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
    contrast_samples = set(g1_ids + g2_ids)

    if args.salmon_tpm:
        tpm_out = os.path.join(args.out_dir, "gene_exp_bySample.tsv")
        transpose_tpm(args.salmon_tpm, sample_order, tpm_out)

    psi_out = os.path.join(args.out_dir, "PSI_bySample.tsv")
    build_psi_matrix(args.rmats_se, g1_ids, g2_ids, psi_out,
                     args.min_samples, args.fdr_cutoff, args.dpsi_cutoff)

    grp_out = os.path.join(args.out_dir, "group_info.tsv")
    write_group_info(sample_order, group_map, grp_out)

    # Restrict group_order to only the contrast samples so downstream steps
    # (prepareGeneMatrixOrdered.py) subset pathway scores to these samples.
    contrast_group_map = {s: g for s, g in group_map.items() if s in contrast_samples}
    ord_out = os.path.join(args.out_dir, "group_order.txt")
    write_group_order(contrast_group_map, ord_out)

    # Write contrast sample list so downstream steps can subset pathway scores.
    samples_out = os.path.join(args.out_dir, "contrast_samples.txt")
    with open(samples_out, "w") as fh:
        fh.write(",".join(g1_ids + g2_ids) + "\n")


if __name__ == "__main__":
    main()
