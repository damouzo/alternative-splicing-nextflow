#!/usr/bin/env python3
"""
Build the navigable sashimi PDF index (shippable deliverable):

  deliverables/contrasts/<comparison_id>/plots/sashimi/sashimi_index.tsv

One row per final PDF. The final PDFs are published next to this index under
deliverables/contrasts/<comparison_id>/plots/sashimi/<EVENT_TYPE>/<file>.pdf,
and pdf_path is stored relative to the deliverables/ root (no absolute cluster
paths, no references to the audit-only raw/ layer), so the table stays valid
wherever the shipped deliverables/ folder is unpacked or renamed.

The rank and gene symbol are parsed from the plot filename written by
rmats2sashimiplot: <rank>_<GENE>_<rest>.pdf

Stdlib only.
"""

import argparse
import csv
import os
import re
import sys

EVENT_TYPES = ('SE', 'A5SS', 'A3SS', 'MXE', 'RI')
RANK_GENE_RE = re.compile(r'^(\d+)_([^_]+)_')


def parse_pdf_name(name):
    match = RANK_GENE_RE.match(name)
    if not match:
        return '', ''
    rank = int(match.group(1))
    gene = match.group(2)
    return rank, gene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-id', required=True)
    parser.add_argument('--sashimi-dir', required=True,
                        help='staged sashimi_out directory')
    parser.add_argument('--out-dir', default='.')
    args = parser.parse_args()

    staged_root = os.path.abspath(args.sashimi_dir)
    rows = []
    for event_type in EVENT_TYPES:
        plot_dir = os.path.join(staged_root, event_type, 'Sashimi_plot')
        if not os.path.isdir(plot_dir):
            continue
        for name in sorted(os.listdir(plot_dir)):
            if not name.endswith('.pdf'):
                continue
            rank, gene = parse_pdf_name(name)
            # Published layout, relative to the deliverables/ root:
            #   contrasts/<comparison_id>/plots/sashimi/<EVENT_TYPE>/<file>.pdf
            pdf_rel = '/'.join([
                'contrasts', args.comparison_id, 'plots', 'sashimi',
                event_type, name,
            ])
            rows.append({
                'comparison_id': args.comparison_id,
                'event_type': event_type,
                'rank': rank,
                'gene_symbol': gene,
                # rmats2sashimiplot cannot draw unannotated splice sites, so
                # filter_rmats_for_sashimi.py excludes de novo events upstream
                # and every indexed plot is annotated. Kept explicit so the
                # browser contract survives a future per-class split.
                'site_class': 'annotated',
                'plot_id': name[:-4],
                'pdf_path': pdf_rel,
            })

    out_path = os.path.join(args.out_dir, 'sashimi_index.tsv')
    columns = ['comparison_id', 'event_type', 'rank', 'gene_symbol',
               'site_class', 'plot_id', 'pdf_path']
    with open(out_path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print('[build_sashimi_index] %d PDFs indexed for %s' % (len(rows), args.comparison_id))
    return 0


if __name__ == '__main__':
    sys.exit(main())
