#!/usr/bin/env python3
"""
Build the shared PEGASAS gene expression matrix for the pathway (KS
enrichment) step.

The KS enrichment score of a sample is a deterministic function of that
sample's gene expression and the gene signature only — it does not depend
on the contrast or on the group. PEGASAS therefore computes pathway scores
ONCE for all samples and the per-contrast correlation step subsets them
(see prepare_pegasas_inputs.py). This matches the original PEGASAS design.

Outputs:
  gene_exp_bySample.tsv — rows=samples, cols=genes (TPM, samples = all in group_info)
  group_info.tsv        — sample_id<TAB>group (all samples, passed through)
"""

import argparse
import csv
import sys
import os


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("salmon_tpm",    help="salmon.merged.gene_tpm.tsv (genes × samples)")
    p.add_argument("group_info_in", help="TSV: sample_id<TAB>group (one row per sample)")
    p.add_argument("--out-dir",      default=".", dest="out_dir", help="Output directory [.]")
    return p.parse_args()


def load_group_info(fin: str) -> list:
    """Return ordered list of (sample_id, group), preserving file order."""
    pairs = []
    with open(fin) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            pairs.append((parts[0], parts[1]))
    return pairs


def transpose_tpm(fin: str, sample_order: list, out_path: str) -> None:
    """
    Convert salmon.merged.gene_tpm.tsv (genes × samples) to PEGASAS format
    (samples × genes, first col = sample_id).

    gene_name is not unique (several Ensembl gene_ids can share a symbol), so
    duplicated symbols are collapsed by keeping the row with the highest mean
    TPM across samples. That keeps the most-expressed gene_id per symbol and
    avoids duplicate matrix columns.

    Keep the collapse rule in sync with prepare_pegasas_inputs.py.transpose_tpm
    (the two scripts share no import path, so the logic is mirrored).
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

    all_tpm_samples = set()
    if genes:
        all_tpm_samples = set(expr[genes[0]].keys())
    valid_samples = [s for s in sample_order if s in all_tpm_samples]

    if not valid_samples:
        sys.exit("[ERROR] No samples matched between TPM file and group info")

    missing = [s for s in sample_order if s not in all_tpm_samples]
    if missing:
        print(f"[WARN] {len(missing)} group_info sample(s) absent from TPM: {missing}",
              file=sys.stderr)

    with open(out_path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SampleID"] + genes)
        for sample in valid_samples:
            row = [sample] + [expr[g].get(sample, 0.0) for g in genes]
            writer.writerow(row)

    print(f"[INFO] Shared gene expression matrix: {len(valid_samples)} samples × "
          f"{len(genes)} genes → {out_path}")


def write_group_info(pairs: list, out_path: str) -> None:
    with open(out_path, "w") as fh:
        for sample, group in pairs:
            fh.write(f"{sample}\t{group}\n")
    print(f"[INFO] Shared group info written ({len(pairs)} samples) → {out_path}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    pairs = load_group_info(args.group_info_in)
    if not pairs:
        sys.exit("[ERROR] group_info is empty")

    sample_order = [s for s, _ in pairs]

    tpm_out = os.path.join(args.out_dir, "gene_exp_bySample.tsv")
    transpose_tpm(args.salmon_tpm, sample_order, tpm_out)

    grp_out = os.path.join(args.out_dir, "group_info.tsv")
    write_group_info(pairs, grp_out)


if __name__ == "__main__":
    main()
