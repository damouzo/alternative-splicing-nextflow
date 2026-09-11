#!/usr/bin/env python3
"""
Build the navigable sashimi PDF index (CORE deliverable):

  plots/sashimi/sashimi_index.tsv

One row per final PDF under sashimi_out/<EVENT_TYPE>/Sashimi_plot/. Paths in
the index are relative to the run outdir (publish_root) so the file stays
valid after the pipeline finishes, regardless of where work/ lives.

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
    parser.add_argument('--publish-root', required=True,
                        help='absolute outdir the index paths are relative to')
    parser.add_argument('--publish-subdir', default='sashimi',
                        help='outdir-relative folder holding this comparison')
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
            pdf_path = os.path.join(plot_dir, name)
            rel_pdf = os.path.relpath(pdf_path, staged_root)
            # Published layout: <outdir>/<publish_subdir>/<comparison_id>/<dir name>/<rel>
            # (the staged dir keeps its original basename, e.g. sashimi_out)
            published = os.path.join(
                args.publish_subdir, args.comparison_id,
                os.path.basename(staged_root), rel_pdf).replace(os.sep, '/')
            rank, gene = parse_pdf_name(name)
            rows.append({
                'comparison_id': args.comparison_id,
                'event_type': event_type,
                'rank': rank,
                'gene_symbol': gene,
                'plot_id': name[:-4],
                'pdf_path': os.path.join(args.publish_root, published).replace(os.sep, '/'),
                'pdf_path_relative_to_outdir': published,
            })

    out_path = os.path.join(args.out_dir, 'sashimi_index.tsv')
    columns = ['comparison_id', 'event_type', 'rank', 'gene_symbol', 'plot_id',
               'pdf_path', 'pdf_path_relative_to_outdir']
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