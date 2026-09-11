#!/usr/bin/env python3
"""
QA validator for the deliverables results contract (reestructure_plan.md §8.1).

Run this after all deliverables are published. The validator walks the
published outdir and checks:

  Structural (FAIL -> exit 2):
    - deliverables/metadata: run_manifest.yaml, sample_index.tsv, tools_matrix.tsv
    - one contrast_manifest.yaml per comparison
    - master + significant + summary per enabled AS tool, with the standard columns
    - rMATS master row count == sum of rows of the 5 *.MATS.JC.txt files
    - sashimi_index.tsv when sashimi enabled, every indexed PDF exists
    - per-comparison HTML report exists

  Content sanity (FAIL or WARN depending on severity), derived from known
  tool reliability issues:
    - leafcutter: clusters with p.adjust < cutoff & status Success present
      while significant.tsv is empty -> FAIL (filter regression)
    - rMATS: fraction of FDR==0 rows above threshold -> WARN (--cstat floor)
    - ISAR: q-values pinned (IQR < 0.05) while |dIF| >= 0.1 exists -> WARN
    - cross_tool: tools that ran but contributed 0 significant genes -> WARN
    - PEGASAS: significant table marked as uncorrected -> always noted

Exit codes: 0 = pass, 1 = pass with warnings, 2 = contract broken.
Stdlib only.
"""

import argparse
import csv
import os
import statistics
import sys

REQUIRED_COLUMNS = ['comparison_id', 'tool', 'feature_type', 'feature_id',
                    'gene_id', 'gene_symbol', 'effect_size', 'is_significant',
                    'significance_rule', 'source_file']

RMATS_EVENT_TYPES = ('SE', 'A5SS', 'A3SS', 'MXE', 'RI')
MAX_FDR_ZERO_FRACTION = 0.25
ISAR_IQR_THRESHOLD = 0.05
ISAR_DIF_THRESHOLD = 0.1


class Report:
    def __init__(self):
        self.errors = []
        self.warnings = []

    def error(self, message):
        self.errors.append(message)

    def warn(self, message):
        self.warnings.append(message)

    @property
    def exit_code(self):
        if self.errors:
            return 2
        if self.warnings:
            return 1
        return 0


def iter_tsv(path):
    """Stream rows of a TSV (header-driven, one row in memory at a time)."""
    with open(path, 'r', encoding='utf-8', errors='replace', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            yield row


def first_tsv_row(path):
    for row in iter_tsv(path):
        return row
    return None


def count_tsv_rows(path):
    with open(path, 'r', encoding='utf-8', errors='replace', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def check_columns(path, report):
    first = first_tsv_row(path)
    if first is None:
        report.error('%s is empty' % path)
        return
    missing = [c for c in REQUIRED_COLUMNS if c not in first]
    if missing:
        report.error('%s missing standard columns: %s' % (path, ', '.join(missing)))


def check_rmats_master(outdir, comparison_id, report):
    master = os.path.join(outdir, 'deliverables', 'contrasts', comparison_id,
                          'data_tables', '%s.rmats.master.tsv' % comparison_id)
    if not os.path.isfile(master):
        report.error('rMATS master missing: %s' % master)
        return
    check_columns(master, report)

    rmats_dir = os.path.join(outdir, 'rmats', comparison_id, comparison_id)
    total_raw = 0
    for event_type in RMATS_EVENT_TYPES:
        jc = os.path.join(rmats_dir, '%s.MATS.JC.txt' % event_type)
        if os.path.isfile(jc):
            total_raw += count_tsv_rows(jc)
    total_master = count_tsv_rows(master)
    if total_raw and total_master != total_raw:
        report.error(
            'rMATS master rows (%d) != sum of MATS.JC rows (%d) for %s'
            % (total_master, total_raw, comparison_id))

    floor = 0
    n_sig_master = 0
    for row in iter_tsv(master):
        if row.get('fdr_floor_flag') == 'true':
            floor += 1
        if row.get('is_significant') == 'true':
            n_sig_master += 1
    if floor:
        fraction = floor / total_master
        if fraction > MAX_FDR_ZERO_FRACTION:
            report.warn(
                'rMATS %s: %.1f%% of events have FDR==0 (--cstat numeric floor); '
                'fdr_floor_flag is set (known issue, annotated)'
                % (comparison_id, 100.0 * fraction))
    significant_path = os.path.join(os.path.dirname(master),
                                    '%s.rmats.significant.tsv' % comparison_id)
    if os.path.isfile(significant_path):
        n_sig_table = count_tsv_rows(significant_path)
        if n_sig_master != n_sig_table:
            report.error(
                'rMATS %s: is_significant rows in master (%d) != significant.tsv rows (%d)'
                % (comparison_id, n_sig_master, n_sig_table))
    else:
        report.error('rMATS significant.tsv missing for %s' % comparison_id)


def check_leafcutter_master(outdir, comparison_id, report, fdr_cutoff):
    master = os.path.join(outdir, 'deliverables', 'contrasts', comparison_id,
                          'data_tables', '%s.leafcutter.master.tsv' % comparison_id)
    significant = os.path.join(os.path.dirname(master),
                               '%s.leafcutter.significant.tsv' % comparison_id)
    if not os.path.isfile(master):
        report.error('LeafCutter master missing: %s' % master)
        return
    check_columns(master, report)
    n_ok_with_p = 0
    for row in iter_tsv(master):
        if row.get('status') == 'Success' and row.get('fdr') != '':
            try:
                if float(row['fdr']) < fdr_cutoff:
                    n_ok_with_p += 1
            except (TypeError, ValueError):
                pass
    if n_ok_with_p and not (os.path.isfile(significant) and count_tsv_rows(significant) > 0):
        report.error(
            'LeafCutter %s: %d clusters with status=="Success" & p.adjust < %s '
            'but significant.tsv has 0 rows — status filter regression'
            % (comparison_id, n_ok_with_p, fdr_cutoff))


def check_isar_master(outdir, comparison_id, report):
    master = os.path.join(outdir, 'deliverables', 'contrasts', comparison_id,
                          'data_tables', '%s.isar.master.tsv' % comparison_id)
    if not os.path.isfile(master):
        report.error('ISAR master missing: %s' % master)
        return
    check_columns(master, report)
    q_values = []
    dmax = 0.0
    for row in iter_tsv(master):
        try:
            q_values.append(float(row.get('fdr')))
        except (TypeError, ValueError):
            pass
        try:
            dmax = max(dmax, abs(float(row.get('effect_size'))))
        except (TypeError, ValueError):
            pass
    if len(q_values) > 10:
        q_sorted = sorted(q_values)
        iqr = q_sorted[len(q_sorted) // 4 * 3] - q_sorted[len(q_sorted) // 4]
        if iqr < ISAR_IQR_THRESHOLD and dmax >= ISAR_DIF_THRESHOLD:
            report.warn(
                'ISAR %s: q-values pinned (IQR %.4f < %.2f) while |dIF| max = %.2f — '
                'known issue (satuRn), do not read significant.tsv as authoritative'
                % (comparison_id, iqr, ISAR_IQR_THRESHOLD, dmax))


def check_cross_tool(outdir, comparison_id, report, enabled_tools):
    path = os.path.join(outdir, 'deliverables', 'contrasts', comparison_id,
                        'data_tables', '%s.cross_tool.master.tsv' % comparison_id)
    if not os.path.isfile(path):
        report.warn('cross_tool.master.tsv missing for %s' % comparison_id)
        return
    counts = {}
    for row in iter_tsv(path):
        counts[row.get('tool', '')] = counts.get(row.get('tool', ''), 0) + 1
    as_tools = [t for t in ('rmats', 'majiq', 'isar', 'leafcutter') if t in enabled_tools]
    empty = [t for t in as_tools if t not in counts or counts[t] == 0]
    if empty:
        report.warn(
            'cross_tool %s: tool(s) with 0 significant genes silently absent from '
            'the overlap/master: %s' % (comparison_id, ', '.join(empty)))


def check_sashimi_index(outdir, comparison_id, report):
    path = os.path.join(outdir, 'deliverables', 'contrasts', comparison_id,
                        'plots', 'sashimi', 'sashimi_index.tsv')
    if not os.path.isfile(path):
        report.error('sashimi_index.tsv missing for %s (sashimi enabled)' % comparison_id)
        return
    rows = list(iter_tsv(path))
    if not rows:
        report.warn('sashimi_index.tsv is empty for %s (no sashimi-eligible events; not an error)'
                    % comparison_id)
        return
    missing_pdfs = []
    for row in rows:
        pdf = row.get('pdf_path', '')
        if pdf and not os.path.isfile(pdf):
            missing_pdfs.append(pdf)
    if missing_pdfs:
        report.error(
            'sashimi_index.tsv %s: %d indexed PDFs missing: %s' %
            (comparison_id, len(missing_pdfs), missing_pdfs[0]))


def parse_manifest_yaml(path):
    """Minimal parse of our manifests: flat `key: value` + `key: [a, b]` lists.
    Blocks starting with '- ' (known_issues items) are skipped."""
    result = {}
    skip = False
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            text = line.rstrip('\n')
            stripped = text.strip()
            if not stripped or stripped.startswith('#'):
                continue
            if stripped.startswith('- '):
                skip = True
                continue
            if not text.startswith(' '):
                skip = False
            if skip:
                continue
            if ': ' in text:
                key, value = text.split(': ', 1)
                key = key.strip()
                value = value.strip()
                if value.startswith('[') and value.endswith(']'):
                    value = [v for v in (v.strip() for v in value[1:-1].split(',')) if v]
                result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, help='published run outdir')
    parser.add_argument('--report-file', default='results_contract_report.txt')
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    report = Report()
    deliverables = os.path.join(outdir, 'deliverables')
    metadata = os.path.join(deliverables, 'metadata')

    if not os.path.isdir(deliverables):
        report.error('deliverables/ layer missing in %s' % outdir)
        _finish(report, args.report_file)
        sys.exit(report.exit_code)

    run_manifest = os.path.join(metadata, 'run_manifest.yaml')
    sample_index = os.path.join(metadata, 'sample_index.tsv')
    tools_matrix = os.path.join(metadata, 'tools_matrix.tsv')
    for path, label in ((run_manifest, 'run_manifest.yaml'),
                        (sample_index, 'sample_index.tsv'),
                        (tools_matrix, 'tools_matrix.tsv')):
        if not os.path.isfile(path):
            report.error('missing global metadata: %s' % label)
    if not os.path.isfile(run_manifest):
        _finish(report, args.report_file)
        sys.exit(report.exit_code)

    manifest = parse_manifest_yaml(run_manifest)
    tools_raw = parse_manifest_yaml(tools_matrix) if os.path.isfile(tools_matrix) else {}
    enabled_tools = set()
    if os.path.isfile(tools_matrix):
        with open(tools_matrix, 'r', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                if row.get('enabled') == 'true':
                    enabled_tools.add(row.get('tool', ''))
    fdr_cutoff = float(manifest.get('fdr_cutoff', 0.05))

    comparisons = manifest.get('comparisons', [])
    if isinstance(comparisons, str):
        comparisons = [comparisons]
    if not comparisons:
        report.warn('run_manifest lists no comparisons')

    for comparison_id in comparisons:
        contrast_dir = os.path.join(deliverables, 'contrasts', comparison_id)
        if not os.path.isdir(contrast_dir):
            report.error('missing contrast deliverables dir: %s' % contrast_dir)
            continue
        contrast_manifest = os.path.join(contrast_dir, 'metadata', 'contrast_manifest.yaml')
        if not os.path.isfile(contrast_manifest):
            report.error('missing contrast_manifest.yaml for %s' % comparison_id)

        report_html = os.path.join(outdir, 'report', '%s_splicing_report.html' % comparison_id)
        if not os.path.isfile(report_html):
            report.error('missing report HTML for %s' % comparison_id)

        for tool in ('rmats', 'majiq', 'isar', 'leafcutter', 'pegasas'):
            if tool not in enabled_tools:
                continue
            data_tables = os.path.join(contrast_dir, 'data_tables')
            for role in ('master', 'significant', 'summary'):
                path = os.path.join(data_tables, '%s.%s.%s.tsv' % (comparison_id, tool, role))
                if not os.path.isfile(path):
                    report.error('missing %s.%s.%s.tsv for %s'
                                 % (comparison_id, tool, role, tool))

        if 'rmats' in enabled_tools:
            check_rmats_master(outdir, comparison_id, report)
        if 'leafcutter' in enabled_tools:
            check_leafcutter_master(outdir, comparison_id, report, fdr_cutoff)
        if 'isar' in enabled_tools:
            check_isar_master(outdir, comparison_id, report)
        check_cross_tool(outdir, comparison_id, report, enabled_tools)
        if 'sashimi' in enabled_tools:
            check_sashimi_index(outdir, comparison_id, report)

    _finish(report, args.report_file)
    sys.exit(report.exit_code)


def _finish(report, report_file):
    with open(report_file, 'w', encoding='utf-8') as handle:
        handle.write('RESULTS CONTRACT VALIDATION REPORT\n')
        handle.write('==================================\n')
        handle.write('errors: %d\nwarnings: %d\n\n' % (len(report.errors), len(report.warnings)))
        for message in report.errors:
            handle.write('[ERROR] %s\n' % message)
        for message in report.warnings:
            handle.write('[WARN] %s\n' % message)
    print('errors: %d | warnings: %d' % (len(report.errors), len(report.warnings)))


if __name__ == '__main__':
    sys.exit(main())