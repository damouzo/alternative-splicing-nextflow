#!/usr/bin/env python3
"""Convert GTF to GFF3, inferring gene records from transcript coordinates.

The input GTF may lack gene-level records; this script derives them by
grouping transcript coordinates per gene_id and emits a proper
gene -> transcript -> exon hierarchy required by MAJIQ v3.

The conversion preserves gene_name when present.
By default, missing gene_name values are tolerated.
Use --strict-gene-name to fail when any gene_id lacks gene_name.
"""

import argparse
import re
import sys
import os
from collections import OrderedDict
from urllib.parse import quote


def get_attr(attrs, key):
    m = re.search(r'(?:^|;)\s*' + key + r'\s+"([^"]+)"', attrs)
    return m.group(1) if m else None


def gff3_escape(value):
    return quote(str(value), safe='-_.~')


def main(gtf_file, gff3_file, strict_gene_name=False):
    genes = OrderedDict()   # gene_id -> coord dict
    records = []            # (feat, chrom, source, start, end, score, strand, tx_id, gene_id)
    gene_names = {}         # gene_id -> gene_name

    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feat, start, end, score, strand, _frame, attrs = parts
            gene_id = get_attr(attrs, 'gene_id')
            if not gene_id:
                continue
            gene_name = get_attr(attrs, 'gene_name')

            if gene_name:
                if gene_id in gene_names and gene_names[gene_id] != gene_name:
                    print(
                        f"ERROR: Inconsistent gene_name for gene_id '{gene_id}': "
                        f"'{gene_names[gene_id]}' vs '{gene_name}'",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                gene_names[gene_id] = gene_name

            if gene_id not in genes:
                genes[gene_id] = {
                    'chrom': chrom,
                    'source': source,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                }
            else:
                g = genes[gene_id]
                g['start'] = min(g['start'], int(start))
                g['end']   = max(g['end'],   int(end))

            if feat == 'transcript':
                transcript_id = get_attr(attrs, 'transcript_id') or gene_id + '_t'
                records.append(('transcript', chrom, source, int(start), int(end),
                                 score, strand, transcript_id, gene_id))

            elif feat == 'exon':
                transcript_id = get_attr(attrs, 'transcript_id') or gene_id + '_t'
                records.append(('exon', chrom, source, int(start), int(end),
                                 score, strand, transcript_id, gene_id))

    missing_gene_names = [gid for gid in genes if gid not in gene_names]
    if missing_gene_names and strict_gene_name:
        preview = ', '.join(missing_gene_names[:10])
        suffix = '' if len(missing_gene_names) <= 10 else f" ... (+{len(missing_gene_names) - 10} more)"
        print(
            "ERROR: Missing required gene_name for gene_id(s): "
            f"{preview}{suffix}",
            file=sys.stderr,
        )
        sys.exit(1)

    if missing_gene_names and not strict_gene_name:
        preview = ', '.join(missing_gene_names[:10])
        suffix = '' if len(missing_gene_names) <= 10 else f" ... (+{len(missing_gene_names) - 10} more)"
        print(
            "WARNING: Missing gene_name for gene_id(s); writing gene records without Name attribute: "
            f"{preview}{suffix}",
            file=sys.stderr,
        )

    written_genes = set()
    with open(gff3_file, 'w') as out:
        out.write('##gff-version 3\n')
        for rec in records:
            feat, chrom, source, start, end, score, strand, transcript_id, gene_id = rec
            gene_name = gene_names.get(gene_id)

            # Emit gene record once, immediately before its first transcript
            if gene_id not in written_genes:
                g = genes[gene_id]
                gene_attrs = ['ID=' + gff3_escape(gene_id)]
                if gene_name:
                    gene_attrs.append('Name=' + gff3_escape(gene_name))
                out.write('\t'.join([
                    g['chrom'], g['source'], 'gene',
                    str(g['start']), str(g['end']),
                    '.', g['strand'], '.',
                    ';'.join(gene_attrs),
                ]) + '\n')
                written_genes.add(gene_id)

            if feat == 'transcript':
                tx_attrs = [
                    'ID=' + gff3_escape(transcript_id),
                    'Parent=' + gff3_escape(gene_id),
                ]
                if gene_name:
                    tx_attrs.append('gene_name=' + gff3_escape(gene_name))
                out.write('\t'.join([
                    chrom, source, 'transcript',
                    str(start), str(end), score, strand, '.',
                    ';'.join(tx_attrs),
                ]) + '\n')
            elif feat == 'exon':
                out.write('\t'.join([
                    chrom, source, 'exon',
                    str(start), str(end), score, strand, '.',
                    'Parent=' + gff3_escape(transcript_id),
                ]) + '\n')

    n_with_gene_name = len(gene_names)
    n_total_genes = len(written_genes)
    n_without_gene_name = n_total_genes - n_with_gene_name
    print(
        f"Written {n_total_genes} genes to {gff3_file} "
        f"({n_with_gene_name} with gene_name, {n_without_gene_name} without gene_name)",
        flush=True,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert GTF to GFF3 for MAJIQ v3")
    parser.add_argument('input_gtf', help='Input GTF file')
    parser.add_argument('output_gff3', help='Output GFF3 file')
    parser.add_argument(
        '--strict-gene-name',
        action='store_true',
        help='Fail if any gene_id is missing gene_name',
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input_gtf):
        print(f"ERROR: GTF file not found: {args.input_gtf}", file=sys.stderr)
        sys.exit(1)
    main(args.input_gtf, args.output_gff3, strict_gene_name=args.strict_gene_name)
