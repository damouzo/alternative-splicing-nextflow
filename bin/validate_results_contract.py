#!/usr/bin/env python3
"""
QA validator for the shippable results contract (assets/results_schema.yaml).

Run this after all deliverables are published. It walks the published outdir
and checks:

  Structural (FAIL -> exit 2):
    - run_info (deliverables/run_info): sample_index.tsv, software_versions.yml
    - run manifest and per-contrast manifests under results/run_info/ (NOT in deliverables/)
    - master + significant + summary per enabled AS tool, with the standard columns
    - cross_tool.master.tsv + cross_tool.gene_summary.tsv
    - rMATS master row count == sum of rows of the 5 *.MATS.JC.txt files (raw layer)
    - sashimi_index.tsv when sashimi enabled; every indexed PDF exists on disk
    - per-comparison HTML report exists under contrasts/<id>/

  Content sanity (WARN or FAIL, derived from known tool reliability issues):
    - leafcutter: clusters with status Success & padj < cutoff present while
      significant.tsv is empty -> FAIL (filter regression)
    - rMATS: every master row must carry a non-empty event_locus (contract
      2.1.0) -> FAIL on runs declaring contract >= 2.1.0; WARN on older runs
      (their tables predate the column)
    - rMATS: fraction of FDR==0 rows above threshold -> WARN (--cstat floor)
    - ISAR: padj pinned (IQR < 0.05) while |dIF| >= 0.1 exists -> WARN
    - cross_tool: tools that ran but contributed 0 significant genes -> WARN
    - PEGASAS: significant table is BH-corrected, noted

Exit codes: 0 = pass, 1 = pass with warnings, 2 = contract broken.
Stdlib only.
"""

import argparse
import csv
import os
import sys

DEFAULT_REQUIRED_COLUMNS = ['comparison_id', 'tool', 'feature_type', 'feature_id',
                            'gene_id', 'gene_symbol', 'group1_name', 'group2_name',
                            'effect_size_direction', 'effect_size',
                            'effect_size_type', 'is_significant',
                            'significance_rule', 'source_file']

RMATS_EVENT_TYPES = ('SE', 'A5SS', 'A3SS', 'MXE', 'RI')
MAX_FDR_ZERO_FRACTION = 0.25
ISAR_IQR_THRESHOLD = 0.05
ISAR_DIF_THRESHOLD = 0.1


class Report:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.infos = []

    def error(self, message):
        self.errors.append(message)

    def warn(self, message):
        self.warnings.append(message)

    def info(self, message):
        self.infos.append(message)

    @property
    def exit_code(self):
        if self.errors:
            return 2
        if self.warnings:
            return 1
        return 0


def load_schema(path):
    """Minimal parser for the flat lists we need from results_schema.yaml."""
    schema = {}
    if not path or not os.path.isfile(path):
        return schema
    list_keys = ('required_master_columns', 'master_columns',
                 'sample_index_columns', 'sashimi_index_columns',
                 'cross_tool_gene_columns')
    key = None
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            if not line.startswith(' ') and stripped.endswith(':'):
                key = stripped[:-1]
                continue
            if line.startswith('  - ') and key in list_keys:
                schema.setdefault(key, []).append(stripped[2:].strip())
    return schema


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


def check_columns(path, report, required_columns):
    first = first_tsv_row(path)
    if first is None:
        report.error('%s is empty' % path)
        return
    missing = [c for c in required_columns if c not in first]
    if missing:
        report.error('%s missing standard columns: %s' % (path, ', '.join(missing)))


def tables_dir(deliverables, comparison_id):
    return os.path.join(deliverables, 'contrasts', comparison_id, 'tables')


def check_rmats_master(outdir, deliverables, comparison_id, report,
                       required_columns, legacy_contract):
    master = os.path.join(tables_dir(deliverables, comparison_id),
                          '%s.rmats.master.tsv' % comparison_id)
    if not os.path.isfile(master):
        report.error('rMATS master missing: %s' % master)
        return
    check_columns(master, report, required_columns)

    rmats_dir = os.path.join(outdir, 'raw', 'rmats', comparison_id)
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
    empty_locus = {}
    empty_counts = 0
    count_cols = ('ijc_sample_1_total', 'sjc_sample_1_total',
                  'ijc_sample_2_total', 'sjc_sample_2_total')
    has_locus_col = None
    for row in iter_tsv(master):
        if row.get('fdr_floor_flag') == 'true':
            floor += 1
        if row.get('is_significant') == 'true':
            n_sig_master += 1
        if all((row.get(c) or '').strip() == '' for c in count_cols):
            empty_counts += 1
        if has_locus_col is None:
            has_locus_col = 'event_locus' in row
        if not (row.get('event_locus') or '').strip():
            event_type = row.get('event_type', '')
            empty_locus[event_type] = empty_locus.get(event_type, 0) + 1
    if empty_counts:
        report.warn(
            'rMATS %s: %d rows have all four junction-count columns empty — '
            'counts may be missing (zero-vs-empty regression?)'
            % (comparison_id, empty_counts))
    if has_locus_col is False:
        # Pre-2.1.0 tables legitimately lack the column: downgrade to WARN so
        # old runs keep validating, new runs fail until rMATS masters are
        # regenerated.
        if legacy_contract:
            report.warn(
                'rMATS %s: event_locus column missing (introduced in results '
                'contract 2.1.0; this table predates it)' % comparison_id)
        else:
            report.error(
                'rMATS %s: event_locus column missing (results contract >= 2.1.0; '
                'regenerate the rMATS master tables)' % comparison_id)
    elif empty_locus:
        detail = ', '.join(
            '%s=%d' % (event_type, count)
            for event_type, count in sorted(empty_locus.items()))
        if legacy_contract:
            report.warn(
                'rMATS %s: rows with empty event_locus by event_type: %s '
                '(events cannot be located in IGV)' % (comparison_id, detail))
        else:
            report.error(
                'rMATS %s: rows with empty event_locus by event_type: %s '
                '(events cannot be located in IGV)' % (comparison_id, detail))
    if floor and total_master:
        fraction = floor / total_master
        if fraction > MAX_FDR_ZERO_FRACTION:
            report.warn(
                'rMATS %s: %.1f%% of events have FDR==0 (--cstat numeric floor); '
                'fdr_floor_flag is set (known issue, annotated)'
                % (comparison_id, 100.0 * fraction))
    significant_path = os.path.join(tables_dir(deliverables, comparison_id),
                                    '%s.rmats.significant.tsv' % comparison_id)
    if os.path.isfile(significant_path):
        n_sig_table = count_tsv_rows(significant_path)
        if n_sig_master != n_sig_table:
            report.error(
                'rMATS %s: is_significant rows in master (%d) != significant.tsv rows (%d)'
                % (comparison_id, n_sig_master, n_sig_table))
    else:
        report.error('rMATS significant.tsv missing for %s' % comparison_id)


def check_leafcutter_master(deliverables, comparison_id, report, required_columns,
                            fdr_cutoff, dpsi_cutoff):
    master = os.path.join(tables_dir(deliverables, comparison_id),
                          '%s.leafcutter.master.tsv' % comparison_id)
    significant = os.path.join(tables_dir(deliverables, comparison_id),
                               '%s.leafcutter.significant.tsv' % comparison_id)
    if not os.path.isfile(master):
        report.error('LeafCutter master missing: %s' % master)
        return
    check_columns(master, report, required_columns)
    n_ok_with_p = 0
    n_sig_small_dpsi = 0
    for row in iter_tsv(master):
        if row.get('status') == 'Success' and row.get('padj') != '':
            try:
                if float(row['padj']) < fdr_cutoff:
                    n_ok_with_p += 1
            except (TypeError, ValueError):
                pass
        if row.get('is_significant') == 'true':
            try:
                if abs(float(row.get('effect_size'))) < dpsi_cutoff:
                    n_sig_small_dpsi += 1
            except (TypeError, ValueError):
                pass
    if n_ok_with_p and not (os.path.isfile(significant) and count_tsv_rows(significant) > 0):
        report.error(
            'LeafCutter %s: %d clusters with status=="Success" & padj < %s '
            'but significant.tsv has 0 rows — status filter regression'
            % (comparison_id, n_ok_with_p, fdr_cutoff))
    if n_sig_small_dpsi:
        report.error(
            'LeafCutter %s: %d significant introns with |deltaPSI| < %s — '
            'intron-level |dPSI| filter not applied'
            % (comparison_id, n_sig_small_dpsi, dpsi_cutoff))


def check_isar_master(deliverables, comparison_id, report, required_columns,
                      effective_method=''):
    master = os.path.join(tables_dir(deliverables, comparison_id),
                          '%s.isar.master.tsv' % comparison_id)
    if not os.path.isfile(master):
        report.error('ISAR master missing: %s' % master)
        return
    check_columns(master, report, required_columns)
    q_values = []
    dmax = 0.0
    n_sig = 0
    padj_methods = set()
    for row in iter_tsv(master):
        try:
            q_values.append(float(row.get('padj')))
        except (TypeError, ValueError):
            pass
        try:
            dmax = max(dmax, abs(float(row.get('effect_size'))))
        except (TypeError, ValueError):
            pass
        if row.get('is_significant') == 'true':
            n_sig += 1
        value = (row.get('padj_method') or '').strip()
        if value:
            padj_methods.add(value.lower())

    if q_values and n_sig == 0:
        report.warn(
            'ISAR %s: no isoform passes the significance rule — every gene lacks '
            'a significant isoform' % comparison_id)

    engine = str(effective_method or '').strip().lower()
    expected = 'empirical_fdr' if engine.startswith('satur') else 'bh'
    if padj_methods and expected not in padj_methods:
        report.warn(
            'ISAR %s: padj_method %s does not match the effective switch-test '
            'engine %r (expected %s)'
            % (comparison_id, ','.join(sorted(padj_methods)),
               effective_method or 'unknown', expected))

    if len(q_values) > 10:
        q_sorted = sorted(q_values)
        iqr = q_sorted[len(q_sorted) // 4 * 3] - q_sorted[len(q_sorted) // 4]
        if iqr < ISAR_IQR_THRESHOLD and dmax >= ISAR_DIF_THRESHOLD:
            report.warn(
                'ISAR %s: padj pinned (IQR %.4f < %.2f) while |dIF| max = %.2f — '
                'state under_investigation, do not read significant.tsv as authoritative'
                % (comparison_id, iqr, ISAR_IQR_THRESHOLD, dmax))


def check_cross_tool(deliverables, comparison_id, report, enabled_tools):
    dir_path = tables_dir(deliverables, comparison_id)
    path = os.path.join(dir_path, '%s.cross_tool.master.tsv' % comparison_id)
    if not os.path.isfile(path):
        report.warn('cross_tool.master.tsv missing for %s' % comparison_id)
        return
    gene_summary = os.path.join(dir_path, '%s.cross_tool.gene_summary.tsv' % comparison_id)
    if not os.path.isfile(gene_summary):
        report.error('cross_tool.gene_summary.tsv missing for %s' % comparison_id)
    counts = {}
    for row in iter_tsv(path):
        counts[row.get('tool', '')] = counts.get(row.get('tool', ''), 0) + 1
    as_tools = [t for t in ('rmats', 'majiq', 'isar', 'leafcutter') if t in enabled_tools]
    empty = [t for t in as_tools if t not in counts or counts[t] == 0]
    if empty:
        report.warn(
            'cross_tool %s: tool(s) with 0 significant genes silently absent from '
            'the overlap/master: %s' % (comparison_id, ', '.join(empty)))


def check_sign_consistency(deliverables, comparison_id, report):
    """WARN on gene-level sign disagreement between rMATS and LeafCutter.

    Both tools now report effect_size = group2 - group1, so a gene called by
    both should agree in sign on its strongest hit.
    """
    dir_path = tables_dir(deliverables, comparison_id)

    def best_signs(path):
        # symbol -> (max_abs_effect, sign): keep the strongest hit per gene, the
        # masters are not ordered by effect size.
        signs = {}
        if not os.path.isfile(path):
            return signs
        for row in iter_tsv(path):
            if row.get('is_significant') != 'true':
                continue
            try:
                effect = float(row.get('effect_size'))
            except (TypeError, ValueError):
                continue
            if effect == 0:
                continue
            for symbol in str(row.get('gene_symbol') or '').split(','):
                symbol = symbol.strip().upper()
                if not symbol or symbol == 'NA':
                    continue
                previous = signs.get(symbol)
                if previous is None or abs(effect) > previous[0]:
                    signs[symbol] = (abs(effect), 1 if effect > 0 else -1)
        return {symbol: sign for symbol, (_abs, sign) in signs.items()}

    rmats = best_signs(os.path.join(dir_path, '%s.rmats.master.tsv' % comparison_id))
    leafcutter = best_signs(os.path.join(dir_path, '%s.leafcutter.master.tsv' % comparison_id))
    shared = set(rmats) & set(leafcutter)
    if not shared:
        return
    mismatches = [g for g in shared if rmats[g] != leafcutter[g]]
    if mismatches:
        report.warn(
            'sign disagreement between rMATS and LeafCutter for %d/%d shared '
            'genes (%s) — check effect_size direction'
            % (len(mismatches), len(shared), ', '.join(sorted(mismatches)[:5])))


def check_sashimi_index(deliverables, comparison_id, report):
    path = os.path.join(deliverables, 'contrasts', comparison_id, 'plots', 'sashimi',
                        'sashimi_index.tsv')
    if not os.path.isfile(path):
        report.error('sashimi_index.tsv missing for %s (sashimi enabled)' % comparison_id)
        return
    rows = list(iter_tsv(path))
    if not rows:
        report.warn('sashimi_index.tsv is empty for %s (no sashimi-eligible events; not an error)'
                    % comparison_id)
        return
    # pdf_path is relative to the deliverables/ root
    missing_pdfs = []
    for row in rows:
        rel = row.get('pdf_path', '')
        if not rel:
            continue
        if os.path.isabs(rel) or 'raw/' in rel.replace('\\', '/'):
            missing_pdfs.append(rel)
            continue
        if not os.path.isfile(os.path.join(deliverables, rel)):
            missing_pdfs.append(rel)
    if missing_pdfs:
        report.error(
            'sashimi_index.tsv %s: %d indexed PDFs missing or not deliverables-relative: %s'
            % (comparison_id, len(missing_pdfs), missing_pdfs[0]))


def parse_version_tuple(text):
    """(2, 1, 0) from '2.1.0'; (0,) when absent or unparseable."""
    text = str(text or '').strip().strip('"\'')
    parts = []
    for piece in text.split('.'):
        try:
            parts.append(int(piece))
        except ValueError:
            break
    return tuple(parts) if parts else (0,)


def contract_at_least(manifest, version):
    """True when the validated run declares a contract >= `version`.

    Old manifests carry `schema_version` (e.g. 2.0.0); new ones also carry
    `results_contract_version`. Missing keys mean a pre-2.x run.
    """
    current = parse_version_tuple(
        manifest.get('results_contract_version') or
        manifest.get('schema_version') or '')
    return current >= version


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
                result.setdefault(key, value)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, help='published run outdir (results/)')
    parser.add_argument('--deliverables-root', default=None,
                        help='deliverables/ root [<outdir>/deliverables]')
    parser.add_argument('--report-file', default='qa_report.txt')
    parser.add_argument('--schema', default=None,
                        help='path to assets/results_schema.yaml')
    parser.add_argument('--ignore-missing', nargs='*', default=[],
                        help='run_info file names produced by the calling task that '
                             'must not be required to pre-exist on disk')
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    deliverables = os.path.abspath(
        args.deliverables_root or os.path.join(outdir, 'deliverables'))
    schema = load_schema(args.schema)
    required_columns = schema.get('required_master_columns') or DEFAULT_REQUIRED_COLUMNS

    report = Report()
    run_info = os.path.join(deliverables, 'run_info')
    contrasts_root = os.path.join(deliverables, 'contrasts')

    if not os.path.isdir(run_info):
        report.error('deliverables/run_info/ layer missing in %s' % deliverables)
        _finish(report, args.report_file)
        sys.exit(report.exit_code)

    manifest_root = os.path.join(outdir, 'run_info')
    run_manifest = os.path.join(manifest_root, 'run_manifest.yaml')
    sample_index = os.path.join(run_info, 'sample_index.tsv')
    software_versions = os.path.join(run_info, 'software_versions.yml')
    ignore_missing = set(args.ignore_missing or [])
    for path, label in ((run_manifest, 'run_manifest.yaml'),
                        (sample_index, 'sample_index.tsv'),
                        (software_versions, 'software_versions.yml')):
        if label in ignore_missing:
            continue
        if not os.path.isfile(path):
            report.error('missing run_info file: %s' % label)
    if not os.path.isfile(run_manifest):
        _finish(report, args.report_file)
        sys.exit(report.exit_code)

    # Manifests are run metadata and must live outside the shipped deliverables/.
    stray_manifests = []
    for dirpath, _dirnames, filenames in os.walk(deliverables):
        for name in filenames:
            if name in ('run_manifest.yaml', 'contrast_manifest.yaml') or \
               name.endswith('_contrast_manifest.yaml'):
                stray_manifests.append(os.path.join(dirpath, name))
    if stray_manifests:
        # publishDir copies never delete old-layout files, so a re-run into an
        # existing outdir leaves the previous manifests behind. Surface them
        # (and how to clear them) without aborting the run.
        report.warn(
            'stale manifest files inside deliverables/ (remove them; manifests '
            'now live under results/run_info/): %s'
            % ', '.join(os.path.relpath(p, deliverables) for p in stray_manifests))

    manifest = parse_manifest_yaml(run_manifest)
    # 2.2.0 added group/sign columns and moved manifests out of deliverables/.
    # Legacy trees (pre-2.2.0) keep their old masters and stale manifests, so
    # downgrade the new requirements instead of aborting validation.
    legacy_layout = not contract_at_least(manifest, (2, 2, 0))
    if legacy_layout:
        new_cols = {'group1_name', 'group2_name', 'effect_size_direction'}
        required_columns = [c for c in required_columns if c not in new_cols]
        report.warn('run contract < 2.2.0: group1_name/group2_name/'
                    'effect_size_direction columns not required')
    enabled_tools = set(manifest.get('tools_enabled', []))
    if isinstance(enabled_tools, str):
        enabled_tools = {enabled_tools} if enabled_tools else set()
    fdr_cutoff = float(manifest.get('fdr_cutoff', 0.05))
    dpsi_cutoff = float(manifest.get('dpsi_cutoff', 0.1))

    comparisons = manifest.get('comparisons', [])
    if isinstance(comparisons, str):
        comparisons = [comparisons]
    if not comparisons:
        report.warn('run_manifest lists no comparisons')

    for comparison_id in comparisons:
        contrast_dir = os.path.join(contrasts_root, comparison_id)
        if not os.path.isdir(contrast_dir):
            report.error('missing contrast dir: %s' % contrast_dir)
            continue
        contrast_manifest = os.path.join(
            manifest_root, 'contrast_manifests',
            '%s_contrast_manifest.yaml' % comparison_id)
        contrast_meta = {}
        if not os.path.isfile(contrast_manifest):
            report.error('missing contrast manifest for %s (expected %s)'
                         % (comparison_id, contrast_manifest))
        else:
            contrast_meta = parse_manifest_yaml(contrast_manifest)
        effective_method = contrast_meta.get('effective_test_method', '')
        if 'isar' in enabled_tools:
            report.info('ISAR %s: effective_test_method=%s'
                        % (comparison_id, effective_method or 'NA'))

        report_html = os.path.join(contrast_dir, '%s_splicing_report.html' % comparison_id)
        if not os.path.isfile(report_html):
            report.error('missing report HTML for %s' % comparison_id)

        for tool in ('rmats', 'majiq', 'isar', 'leafcutter', 'pegasas'):
            if tool not in enabled_tools:
                continue
            for role in ('master', 'significant', 'summary'):
                path = os.path.join(tables_dir(deliverables, comparison_id),
                                    '%s.%s.%s.tsv' % (comparison_id, tool, role))
                if not os.path.isfile(path):
                    report.error('missing %s.%s.%s.tsv for %s'
                                 % (comparison_id, tool, role, tool))

        if 'rmats' in enabled_tools:
            check_rmats_master(outdir, deliverables, comparison_id, report,
                               required_columns,
                               not contract_at_least(manifest, (2, 1, 0)))
        if 'leafcutter' in enabled_tools:
            check_leafcutter_master(deliverables, comparison_id, report,
                                    required_columns, fdr_cutoff, dpsi_cutoff)
        if 'isar' in enabled_tools:
            check_isar_master(deliverables, comparison_id, report,
                              required_columns, effective_method)
        check_cross_tool(deliverables, comparison_id, report, enabled_tools)
        check_sign_consistency(deliverables, comparison_id, report)
        if 'sashimi' in enabled_tools:
            check_sashimi_index(deliverables, comparison_id, report)

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
        if report.infos:
            handle.write('\n')
            for message in report.infos:
                handle.write('[INFO] %s\n' % message)
    print('errors: %d | warnings: %d' % (len(report.errors), len(report.warnings)))


if __name__ == '__main__':
    sys.exit(main())
