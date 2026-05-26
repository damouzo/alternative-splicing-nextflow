#!/usr/bin/env python3
"""Convert GTF to GFF3, inferring gene records from transcript coordinates.

The input GTF may lack gene-level records; this script derives them by
grouping transcript coordinates per gene_id and emits a proper
gene -> transcript -> exon hierarchy required by MAJIQ v3.
"""

import re
import sys
import os
from collections import OrderedDict


def get_attr(attrs, key):
    m = re.search(r'(?:^|;)\s*' + key + r'\s+"([^"]+)"', attrs)
    return m.group(1) if m else None


def main(gtf_file, gff3_file):
    genes = OrderedDict()   # gene_id -> coord dict
    records = []            # (feat, chrom, source, start, end, score, strand, tx_id, gene_id)

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

            if feat == 'transcript':
                if gene_id not in genes:
                    genes[gene_id] = {
                        'chrom': chrom, 'source': source,
                        'start': int(start), 'end': int(end),
                        'strand': strand,
                    }
                else:
                    g = genes[gene_id]
                    g['start'] = min(g['start'], int(start))
                    g['end']   = max(g['end'],   int(end))
                transcript_id = get_attr(attrs, 'transcript_id') or gene_id + '_t'
                records.append(('transcript', chrom, source, int(start), int(end),
                                 score, strand, transcript_id, gene_id))

            elif feat == 'exon':
                transcript_id = get_attr(attrs, 'transcript_id') or gene_id + '_t'
                records.append(('exon', chrom, source, int(start), int(end),
                                 score, strand, transcript_id, gene_id))

    written_genes = set()
    with open(gff3_file, 'w') as out:
        out.write('##gff-version 3\n')
        for rec in records:
            feat, chrom, source, start, end, score, strand, transcript_id, gene_id = rec

            # Emit gene record once, immediately before its first transcript
            if gene_id not in written_genes:
                g = genes[gene_id]
                out.write('\t'.join([
                    g['chrom'], g['source'], 'gene',
                    str(g['start']), str(g['end']),
                    '.', g['strand'], '.',
                    'ID=' + gene_id,
                ]) + '\n')
                written_genes.add(gene_id)

            if feat == 'transcript':
                out.write('\t'.join([
                    chrom, source, 'transcript',
                    str(start), str(end), score, strand, '.',
                    'ID=' + transcript_id + ';Parent=' + gene_id,
                ]) + '\n')
            elif feat == 'exon':
                out.write('\t'.join([
                    chrom, source, 'exon',
                    str(start), str(end), score, strand, '.',
                    'Parent=' + transcript_id,
                ]) + '\n')

    print(f"Written {len(written_genes)} genes to {gff3_file}", flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.gtf> <output.gff3>", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(sys.argv[1]):
        print(f"ERROR: GTF file not found: {sys.argv[1]}", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
