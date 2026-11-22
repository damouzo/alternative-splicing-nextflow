#!/usr/bin/env python3
"""
Build the navigable sashimi PDF index (shippable deliverable):

  deliverables/contrasts/<comparison_id>/plots/rmats/sashimi/sashimi_index.tsv

One row per final PDF. The final PDFs are published next to this index under
deliverables/contrasts/<comparison_id>/plots/rmats/sashimi/<EVENT_TYPE>/<file>.pdf,
and pdf_path is stored relative to the deliverables/ root (no absolute cluster
paths, no references to the audit-only raw/ layer), so the table stays valid
wherever the shipped deliverables/ folder is unpacked or renamed.

Automatic (top-N) PDFs are named <rank>_<GENE>_<rest>.pdf; their splice-site
class comes from site_classes.tsv keyed by (event_type, rank). Targeted PDFs are
named A_<GENE>_<EVENT_TYPE>_<coords>.pdf and are resolved through
targeted_classes.tsv keyed by plot_id.

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
    """<rank>_<GENE>_<rest>.pdf -> (rank, gene) for automatic plots."""
    match = RANK_GENE_RE.match(name)
    if not match:
        return '', ''
    return int(match.group(1)), match.group(2)


def _has_pdfs(directory):
    return os.path.isdir(directory) and any(
        name.endswith('.pdf') for name in os.listdir(directory))


def _read_tsv(path):
    if not os.path.isfile(path):
        return []
    with open(path, 'r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def load_site_classes(sashimi_root):
    """(event_type, rank) -> 'annotated' | 'de_novo' from the SASHIMI_PLOTS manifest."""
    classes = {}
    for row in _read_tsv(os.path.join(sashimi_root, 'site_classes.tsv')):
        try:
            classes[(row['event_type'], int(row['rank']))] = row.get('site_class') or 'annotated'
        except (KeyError, ValueError, TypeError):
            continue
    return classes


def load_targeted(sashimi_root):
    """plot_id -> metadata row from the SASHIMI_TARGETED manifest."""
    return {row['plot_id']: row
            for row in _read_tsv(os.path.join(sashimi_root, 'targeted_classes.tsv'))
            if row.get('plot_id')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-id', required=True)
    parser.add_argument('--sashimi-dir', required=True,
                        help='staged sashimi_out directory')
    parser.add_argument('--out-dir', default='.')
    args = parser.parse_args()

    staged_root = os.path.abspath(args.sashimi_dir)
    site_classes = load_site_classes(staged_root)
    targeted = load_targeted(staged_root)
    rows = []
    for event_type in EVENT_TYPES:
        # rmats2sashimiplot versions differ: some write into <ETYPE>/Sashimi_plot/,
        # others straight into <ETYPE>/. Mirror RENDER_REPORT and SASHIMI_*.
        plot_dir = os.path.join(staged_root, event_type, 'Sashimi_plot')
        if not _has_pdfs(plot_dir):
            plot_dir = os.path.join(staged_root, event_type)
        if not _has_pdfs(plot_dir):
            continue
        for name in sorted(os.listdir(plot_dir)):
            if not name.endswith('.pdf'):
                continue
            plot_id = name[:-4]
            is_targeted = name.startswith('A_')
            if is_targeted:
                record = targeted.get(plot_id, {})
                gene = record.get('gene', '')
                rank = ''
                site_class = record.get('site_class', 'annotated')
                event_id = record.get('event_id', '')
                deltapsi = record.get('deltapsi', '')
                fdr = record.get('fdr', '')
                selection = 'targeted'
            else:
                rank, gene = parse_pdf_name(name)
                site_class = site_classes.get((event_type, rank), 'annotated')
                event_id = ''
                deltapsi = ''
                fdr = ''
                selection = 'auto'
            # Published layout, relative to the deliverables/ root:
            #   contrasts/<comparison_id>/plots/rmats/sashimi/<EVENT_TYPE>/<file>.pdf
            pdf_rel = '/'.join([
                'contrasts', args.comparison_id, 'plots', 'rmats', 'sashimi',
                event_type, name,
            ])
            rows.append({
                'comparison_id': args.comparison_id,
                'event_type': event_type,
                'rank': rank,
                'gene_symbol': gene,
                'site_class': site_class,
                'plot_id': plot_id,
                'selection': selection,
                'event_id': event_id,
                'deltapsi': deltapsi,
                'fdr': fdr,
                'pdf_path': pdf_rel,
            })

    out_path = os.path.join(args.out_dir, 'sashimi_index.tsv')
    columns = ['comparison_id', 'event_type', 'rank', 'gene_symbol',
               'site_class', 'plot_id', 'selection', 'event_id',
               'deltapsi', 'fdr', 'pdf_path']
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
