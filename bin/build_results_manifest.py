#!/usr/bin/env python3
"""
Build the publication metadata for the shippable deliverables layer:

  run mode:       run_manifest.yaml, sample_index.tsv, sample_paths.tsv, README.md
  contrast mode:  contrast_manifest.yaml per comparison

Run mode input is a base64-encoded JSON payload with:
  samples:  [{sample_id, condition, replicate, comparison_id, group, bam, bai, salmon_dir}]
  pipeline: {pipeline_name, pipeline_version, nextflow_version, session_id, execution_date,
             outdir, results_contract_version (schema_version), publish_level, fdr_cutoff,
             dpsi_cutoff, strandedness, read_length, organism, genome_build, isar_test_method,
             rmats_cstat, tools: [{name, enabled, container}], comparisons: [..],
             tool_reliability: [{tool, status, note}]}

Contrast mode scans the published layout on disk (deterministic), so the
manifest always reflects what actually exists.

sample_index.tsv is the shipped traceability table and carries no absolute
paths; the cluster paths needed for audit live in sample_paths.tsv, which is
published only under raw/_internal/.

YAML is emitted with a minimal hand-rolled writer (no third-party deps).
"""

import argparse
import base64
import csv
import json
import os
import re
import sys


def yaml_scalar(value):
    """Render a scalar safely as YAML (quoted when not a plain number/bool)."""
    if value is None:
        return 'null'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if text == '':
        return '""'
    if text in {'null', 'true', 'false'} or (text[0].isdigit() and not _is_number(text)):
        return json.dumps(text)
    # Quote anything that could break flow lists ([a, b]) or scalar contexts:
    # commas/brackets/braces would split the value otherwise
    if text[0] in '"\'-!&*?|>%@` ,' or any(c in text for c in ':#\n\t,[]{}'):
        return json.dumps(text)
    return text


def _is_number(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def yaml_block(name, values, indent=0):
    """Serialize a list of (key, value) pairs as a nested YAML mapping."""
    pad = ' ' * indent
    sub = ' ' * (indent + 2)
    lines = [f'{pad}{name}:']
    for key, value in values:
        if isinstance(value, list):
            lines.append(f'{sub}{key}: [{", ".join(yaml_scalar(v) for v in value)}]')
        elif isinstance(value, dict):
            lines.append(f'{sub}{key}:')
            for k2, v2 in value.items():
                lines.append(f'{sub}  {k2}: {yaml_scalar(v2)}')
        else:
            lines.append(f'{sub}{key}: {yaml_scalar(value)}')
    return '\n'.join(lines)


def load_payload(args):
    if getattr(args, 'payload_file', None):
        with open(args.payload_file, 'r', encoding='utf-8') as handle:
            return json.load(handle)
    with open(args.payload_base64, 'r', encoding='utf-8') as handle:
        return json.loads(base64.b64decode(handle.read().strip()).decode('utf-8'))


def build_run_manifest(pipeline):
    critical = {
        key: pipeline.get(key)
        for key in ('fdr_cutoff', 'dpsi_cutoff', 'strandedness', 'read_length',
                    'organism', 'genome_build', 'isar_test_method', 'rmats_cstat')
    }
    tools = pipeline.get('tools', [])
    enabled = [t['name'] for t in tools if t.get('enabled')]
    disabled = [t['name'] for t in tools if not t.get('enabled')]

    lines = []
    lines.append('# Run manifest — shippable publication contract (schema_version ' +
                 str(pipeline.get('results_contract_version', 'unknown')) + ')')
    lines.append(yaml_block('run_manifest', [
        ('pipeline_name', pipeline.get('pipeline_name')),
        ('pipeline_version', pipeline.get('pipeline_version')),
        ('schema_version', pipeline.get('results_contract_version')),
        ('nextflow_version', pipeline.get('nextflow_version')),
        ('session_id', pipeline.get('session_id')),
        ('execution_date', pipeline.get('execution_date')),
        ('publish_level', pipeline.get('publish_level')),
        ('shippable_root', 'deliverables/'),
        ('paths_relative_to', 'shippable_root'),
    ]))
    lines.append(yaml_block('params_criticos', list(critical.items())))
    lines.append('tools_enabled: [%s]' % ', '.join(yaml_scalar(v) for v in enabled))
    lines.append('tools_disabled: [%s]' % ', '.join(yaml_scalar(v) for v in disabled))
    lines.append('comparisons: [%s]' % ', '.join(yaml_scalar(v) for v in pipeline.get('comparisons', [])))
    known = pipeline.get('tool_reliability', [])
    if known:
        lines.append('known_issues:')
        for item in known:
            lines.append('- tool: %s' % yaml_scalar(item.get('tool')))
            lines.append('  status: %s' % yaml_scalar(item.get('status')))
            lines.append('  note: %s' % yaml_scalar(item.get('note') or ''))
    with open('run_manifest.yaml', 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


def build_sample_index(samples):
    """Shipped traceability table — no cluster paths."""
    columns = ['sample_id', 'condition', 'replicate', 'comparison_id', 'group']
    with open('sample_index.tsv', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in samples:
            writer.writerow({
                'sample_id': row.get('sample_id', ''),
                'condition': row.get('condition', ''),
                'replicate': row.get('replicate', ''),
                'comparison_id': row.get('comparison_id', ''),
                'group': row.get('group', ''),
            })


def build_sample_paths(samples):
    """Internal audit table with the absolute input paths (raw/_internal only)."""
    columns = ['sample_id', 'comparison_id', 'group', 'bam_path', 'bai_path', 'salmon_dir']
    with open('sample_paths.tsv', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in samples:
            writer.writerow({
                'sample_id': row.get('sample_id', ''),
                'comparison_id': row.get('comparison_id', ''),
                'group': row.get('group', ''),
                'bam_path': row.get('bam', ''),
                'bai_path': row.get('bai', ''),
                'salmon_dir': row.get('salmon_dir', ''),
            })


def build_readme(pipeline):
    tools = pipeline.get('tools', [])
    enabled = [t['name'] for t in tools if t.get('enabled')]
    comparisons = pipeline.get('comparisons', [])
    lines = []
    lines.append('# alternative-splicing analysis results')
    lines.append('')
    lines.append('Pipeline `%s` v%s (schema_version %s).'
                 % (pipeline.get('pipeline_name', 'alternative-splicing-nextflow'),
                    pipeline.get('pipeline_version', 'unknown'),
                    pipeline.get('results_contract_version', 'unknown')))
    lines.append('')
    lines.append('This folder (`deliverables/`) is the compact shipment: rename it,')
    lines.append('zip it and send it as-is. Every internal path is relative to this')
    lines.append('folder, so it keeps working after renaming or moving it.')
    lines.append('')
    lines.append('`results/raw/` (native per-tool output) is audit evidence kept *outside*')
    lines.append('this folder and is not needed to interpret the results.')
    lines.append('')
    lines.append('## Layout')
    lines.append('')
    lines.append('```')
    lines.append('deliverables/')
    lines.append('  README.md                     # this file')
    lines.append('  run_info/                     # what was run and how to read the tables')
    lines.append('    run_manifest.yaml           # pipeline, params, tools, known issues')
    lines.append('    sample_index.tsv            # sample x comparison x group (no cluster paths)')
    lines.append('    software_versions.yml       # consolidated tool versions')
    lines.append('    qa_report.txt               # results-contract validation output')
    lines.append('    pipeline_info/              # Nextflow trace, timeline, report, DAG')
    lines.append('  contrasts/<comparison_id>/    # one folder per comparison')
    lines.append('    <comparison_id>_splicing_report.html')
    lines.append('    contrast_manifest.yaml      # every deliverable of this comparison')
    lines.append('    tables/                     # <id>.<tool>.{master,significant,summary}.tsv')
    lines.append('    plots/sashimi/              # sashimi_index.tsv + <EVENT_TYPE>/*.pdf')
    lines.append('  cross_contrast/pegasas/       # pathway summary, heatmap, UpSet')
    lines.append('```')
    lines.append('')
    lines.append('## How to read the tables')
    lines.append('')
    lines.append('Every `<tool>.master.tsv` starts with the same standard columns:')
    lines.append('')
    lines.append('| column | meaning |')
    lines.append('| --- | --- |')
    lines.append('| comparison_id | contrast this row belongs to |')
    lines.append('| tool | rmats, majiq, isar, leafcutter or pegasas |')
    lines.append('| feature_type, feature_id | event/isoform/intron/pathway identifier |')
    lines.append('| gene_id, gene_symbol | gene annotation |')
    lines.append('| effect_size, effect_size_type | metric and its type (delta_psi, KS_statistic, ...) |')
    lines.append('| pvalue | nominal p-value when the tool provides one |')
    lines.append('| padj, padj_method | adjusted p-value and how it was computed (BH, rmats_cstat, none, ...) |')
    lines.append('| is_significant | whether the row passes the report significance rule |')
    lines.append('| significance_rule | the exact rule applied |')
    lines.append('| source_file | native file the row was derived from |')
    lines.append('')
    lines.append('Tool-native columns are appended after the standard block and are never')
    lines.append('dropped. `significant.tsv` is the subset with `is_significant == true`;')
    lines.append('`summary.tsv` is a long-format count/effect-size recap.')
    lines.append('')
    lines.append('rMATS rows carry `event_locus` (`chr:start-end`, 1-based inclusive — paste')
    lines.append('it straight into IGV) and `event_coords` (native 0-based coordinates named')
    lines.append('per event type, e.g. `long=…;short=…;flanking=…` for A5SS/A3SS). The legacy')
    lines.append('`exon_start_0base…downstream_ee` columns follow the SE event shape and are')
    lines.append('empty for non-SE event types; use `event_locus`/`event_coords` instead.')
    lines.append('')
    lines.append('Prefer `cross_tool.gene_summary.tsv` for a quick gene list: one row per gene')
    lines.append('with the number of tools that called it significant.')
    lines.append('')
    lines.append('All splicing calls in this folder are computational predictions from')
    lines.append('short-read RNA-seq. Top candidates intended for downstream biological')
    lines.append('interpretation or publication should be orthogonally validated')
    lines.append('(e.g., RT-PCR/qPCR across the relevant junctions, or long-read')
    lines.append('sequencing) before being reported as confirmed splicing events.')
    lines.append('')
    lines.append('## This run')
    lines.append('')
    lines.append('- Tools enabled: %s' % (', '.join(enabled) if enabled else 'none'))
    lines.append('- Comparisons: %s' % (', '.join(comparisons) if comparisons else 'none'))
    lines.append('')
    lines.append('Significance thresholds and per-tool reliability caveats are in')
    lines.append('`run_info/run_manifest.yaml` and each `contrast_manifest.yaml`. Respect the')
    lines.append('`known_issues` notes before using a tool table as authoritative.')
    lines.append('')
    with open('README.md', 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(lines))


def group_counts_from_index(deliverables_root, comparison_id):
    """Count samples per group from the published run_info/sample_index.tsv."""
    counts = {'group1_n': '', 'group2_n': ''}
    index_path = os.path.join(deliverables_root, 'run_info', 'sample_index.tsv')
    if not os.path.isfile(index_path):
        return counts
    n1 = n2 = 0
    with open(index_path, 'r', encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row.get('comparison_id') != comparison_id:
                continue
            group = str(row.get('group', '')).strip()
            if group == '1':
                n1 += 1
            elif group == '2':
                n2 += 1
    if n1:
        counts['group1_n'] = n1
    if n2:
        counts['group2_n'] = n2
    return counts


def build_contrast_manifest(comparison_id, group1_name, group2_name,
                            deliverables_root, results_root, pipeline):
    """Scan the published deliverables/ for one comparison and record it.

    Every shippable path is relative to the deliverables/ root (so the manifest
    stays valid after renaming/moving the folder). Native outputs under
    results/raw/ are listed in a separate, explicitly audit-only section and
    are never referenced as shippable.
    """
    cdir = os.path.join(deliverables_root, 'contrasts', comparison_id)
    tables_dir = os.path.join(cdir, 'tables')
    data_tables = {}
    if os.path.isdir(tables_dir):
        pattern = re.compile(
            r'^%s\.([a-z_]+)\.(master|significant|summary)\.tsv$' % re.escape(comparison_id))
        for name in sorted(os.listdir(tables_dir)):
            match = pattern.match(name)
            if not match or match.group(1) == 'cross_tool':
                continue
            role = match.group(2)
            data_tables.setdefault(match.group(1), {})[role] = os.path.join(
                'contrasts', comparison_id, 'tables', name).replace(os.sep, '/')

    plots = {}
    sashimi_index = os.path.join(cdir, 'plots', 'sashimi', 'sashimi_index.tsv')
    if os.path.isfile(sashimi_index):
        plots['sashimi'] = {'index': os.path.join(
            'contrasts', comparison_id, 'plots', 'sashimi',
            'sashimi_index.tsv').replace(os.sep, '/')}

    cross_tool = {}
    cross_master = os.path.join(
        tables_dir, '%s.cross_tool.master.tsv' % comparison_id)
    if os.path.isfile(cross_master):
        cross_tool['master'] = os.path.join(
            'contrasts', comparison_id, 'tables',
            '%s.cross_tool.master.tsv' % comparison_id).replace(os.sep, '/')
    cross_gene = os.path.join(
        tables_dir, '%s.cross_tool.gene_summary.tsv' % comparison_id)
    if os.path.isfile(cross_gene):
        cross_tool['gene_summary'] = os.path.join(
            'contrasts', comparison_id, 'tables',
            '%s.cross_tool.gene_summary.tsv' % comparison_id).replace(os.sep, '/')

    report_path = None
    report_html = os.path.join(cdir, '%s_splicing_report.html' % comparison_id)
    if os.path.isfile(report_html):
        report_path = os.path.join(
            'contrasts', comparison_id,
            '%s_splicing_report.html' % comparison_id).replace(os.sep, '/')

    # Per-tool raw dirs (audit-only, outside deliverables/): existence reflects
    # publish_level. Paths are relative to results/ (deliverables/ parent).
    raw_dirs = {}
    for tool_dir in ('rmats', 'majiq', 'isar', 'leafcutter', 'sashimi', 'pegasas'):
        tool_path = os.path.join(results_root, 'raw', tool_dir, comparison_id)
        if os.path.isdir(tool_path):
            raw_dirs[tool_dir] = os.path.join('raw', tool_dir, comparison_id)\
                .replace(os.sep, '/')

    thresholds = {
        'fdr_cutoff': pipeline.get('fdr_cutoff'),
        'dpsi_cutoff': pipeline.get('dpsi_cutoff'),
        'majiq_probability_threshold': pipeline.get('majiq_probability_threshold'),
        'majiq_dpsi_cutoff': pipeline.get('majiq_dpsi_cutoff'),
        'rmats_rule': 'padj <= fdr_cutoff & |inc_level_difference| >= dpsi_cutoff',
    }
    for key, value in thresholds.items():
        if value is None:
            thresholds[key] = ''

    group_counts = group_counts_from_index(deliverables_root, comparison_id)

    lines = []
    lines.append('# Contrast deliverable manifest - generated by build_results_manifest.py')
    lines.append('# All shippable paths are relative to the deliverables/ root.')
    lines.append(yaml_block('manifest', [
        ('comparison_id', comparison_id),
        ('group1_name', group1_name),
        ('group2_name', group2_name),
        ('group1_n', group_counts['group1_n']),
        ('group2_n', group_counts['group2_n']),
        ('report_html_path', report_path or ''),
    ]))
    lines.append(yaml_block('thresholds_used', list(thresholds.items())))

    if data_tables:
        lines.append('data_tables:')
        for tool in sorted(data_tables):
            lines.append('  %s:' % tool)
            for role, path in sorted(data_tables[tool].items()):
                lines.append('    %s: %s' % (role, yaml_scalar(path)))
    if cross_tool:
        lines.append('cross_tool:')
        for role, path in sorted(cross_tool.items()):
            lines.append('  %s: %s' % (role, yaml_scalar(path)))
    if plots:
        lines.append('plots:')
        for area, items in plots.items():
            lines.append('  %s:' % area)
            for key, path in items.items():
                lines.append('    %s: %s' % (key, yaml_scalar(path)))
    if raw_dirs:
        lines.append('')
        lines.append('# Native tool output kept for audit. NOT shipped with deliverables/;')
        lines.append('# these paths are relative to results/ (the parent of deliverables/).')
        lines.append('raw_outputs_audit_only:')
        for tool in sorted(raw_dirs):
            lines.append('  %s: %s' % (tool, yaml_scalar(raw_dirs[tool])))

    # tools enabled for this run but with no published deliverable tables
    tools_record = pipeline.get('tools', [])
    enabled = [t.get('name') for t in tools_record if t.get('enabled')]
    missing = [t for t in enabled
               if t not in data_tables and t not in ('sashimi',)]
    if missing:
        lines.append('skipped_tools: [%s]' % ', '.join(yaml_scalar(t) for t in sorted(missing)))

    known = pipeline.get('tool_reliability', [])
    if known:
        lines.append('known_issues:')
        for item in known:
            lines.append('- tool: %s' % yaml_scalar(item.get('tool')))
            lines.append('  status: %s' % yaml_scalar(item.get('status')))
            lines.append('  note: %s' % yaml_scalar(item.get('note') or ''))

    with open('contrast_manifest.yaml', 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=['run', 'contrast_disk'])
    parser.add_argument('--payload-file', help='JSON payload file (alternative to --payload-base64)')
    parser.add_argument('--payload-base64', help='base64-encoded JSON payload')
    args = parser.parse_args()

    payload = load_payload(args)
    if args.mode == 'run':
        build_run_manifest(payload.get('pipeline', {}))
        build_sample_index(payload.get('samples', []))
        build_sample_paths(payload.get('samples', []))
        build_readme(payload.get('pipeline', {}))
    else:
        pipeline = payload.get('pipeline', {})
        results_root = payload.get('outdir') or pipeline.get('outdir')
        deliverables_root = payload.get('deliverables_root') or os.path.join(
            results_root, 'deliverables')
        build_contrast_manifest(
            payload['comparison_id'],
            payload.get('group1_name'),
            payload.get('group2_name'),
            deliverables_root,
            results_root,
            pipeline
        )
    print('[build_results_manifest] %s manifest(s) written to %s' % (args.mode, os.getcwd()))


if __name__ == '__main__':
    sys.exit(main())
