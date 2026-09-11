#!/usr/bin/env python3
"""
Build the publication metadata for the deliverables layer:
  - run mode:         run_manifest.yaml, sample_index.tsv, tools_matrix.tsv
  - contrast mode:    contrast_manifest.yaml per comparison

Run mode input is a base64-encoded JSON payload with:
  samples:  [{sample_id, condition, replicate, comparison_id, group, bam, bai, salmon_dir}]
  pipeline: {pipeline_name, pipeline_version, nextflow_version, session_id, execution_date,
             outdir, results_contract_version, publish_level, fdr_cutoff, dpsi_cutoff,
             strandedness, read_length, organism, genome_build, isar_test_method, rmats_cstat,
             tools: [{name, enabled, container}], comparisons: [..],
             tool_reliability: [{tool, status, note}]}

Contrast mode input carries:
  comparison_id
  entries: [{kind, tool, path, description}]  (file paths relative to deliverables root or outdir)
  and the same pipeline block as run mode.

All outputs are written to the current working directory.
YAML is emitted with a minimal hand-rolled writer (no third-party deps).
"""

import argparse
import base64
import json
import os
import re
import sys


TRUE_VALUES = {'true', 'yes', 'y', '1'}


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
    lines.append('# Publication contract for the deliverables layer (results_contract_version ' +
                 str(pipeline.get('results_contract_version', 'unknown')) + ')')
    lines.append(yaml_block('run_manifest', [
        ('pipeline_name', pipeline.get('pipeline_name')),
        ('pipeline_version', pipeline.get('pipeline_version')),
        ('nextflow_version', pipeline.get('nextflow_version')),
        ('session_id', pipeline.get('session_id')),
        ('execution_date', pipeline.get('execution_date')),
        ('outdir', pipeline.get('outdir')),
        ('results_contract_version', pipeline.get('results_contract_version')),
        ('publish_level', pipeline.get('publish_level')),
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
    import csv
    columns = ['sample_id', 'condition', 'replicate', 'comparison_id',
               'group', 'bam_path', 'bai_path', 'salmon_dir']
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
                'bam_path': row.get('bam', ''),
                'bai_path': row.get('bai', ''),
                'salmon_dir': row.get('salmon_dir', ''),
            })


def build_tools_matrix(pipeline):
    import csv
    columns = ['tool', 'enabled', 'container', 'deliverables']
    with open('tools_matrix.tsv', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for tool in pipeline.get('tools', []):
            writer.writerow({
                'tool': tool.get('name', ''),
                'enabled': str(bool(tool.get('enabled'))).lower(),
                'container': tool.get('container', ''),
                'deliverables': 'master,significant,summary' if tool.get('enabled') else '',
            })


def group_counts_from_index(outdir, comparison_id):
    """Count samples per group from the published sample_index.tsv."""
    import csv
    counts = {'group1_n': '', 'group2_n': ''}
    index_path = os.path.join(outdir, 'deliverables', 'metadata', 'sample_index.tsv')
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


def build_contrast_manifest(comparison_id, group1_name, group2_name, outdir, pipeline):
    """Scan the published deliverables dir for one comparison and record them.

    Every file is listed with its outdir-relative path, so the manifest stays
    valid after the run (mirrors what the report Data exports section shows).
    """
    cdir = os.path.join(outdir, 'deliverables', 'contrasts', comparison_id)
    data_tables_dir = os.path.join(cdir, 'data_tables')
    data_tables = {}
    if os.path.isdir(data_tables_dir):
        pattern = re.compile(r'^%s\.([a-z_]+)\.(master|significant|summary)\.tsv$' % re.escape(comparison_id))
        for name in sorted(os.listdir(data_tables_dir)):
            match = pattern.match(name)
            if not match or match.group(1) == 'cross_tool':
                continue
            role = match.group(2)
            data_tables.setdefault(match.group(1), {})[role] = os.path.join(
                'deliverables', 'contrasts', comparison_id, 'data_tables', name)\
                .replace(os.sep, '/')

    plots = {}
    sashimi_index = os.path.join(cdir, 'plots', 'sashimi', 'sashimi_index.tsv')
    if os.path.isfile(sashimi_index):
        plots['sashimi'] = {'index': os.path.join(
            'deliverables', 'contrasts', comparison_id, 'plots', 'sashimi',
            'sashimi_index.tsv').replace(os.sep, '/')}

    cross_tool = None
    cross_path = os.path.join(data_tables_dir,
                              '%s.cross_tool.master.tsv' % comparison_id)
    if os.path.isfile(cross_path):
        cross_tool = {'master': os.path.join(
            'deliverables', 'contrasts', comparison_id, 'data_tables',
            '%s.cross_tool.master.tsv' % comparison_id).replace(os.sep, '/')}

    report_path = None
    report_html = os.path.join(outdir, 'report',
                               '%s_splicing_report.html' % comparison_id)
    if os.path.isfile(report_html):
        report_path = os.path.join('report',
                                   '%s_splicing_report.html' % comparison_id)\
            .replace(os.sep, '/')

    # per-tool RAW dirs that were published (existence reflects publish_level)
    raw_dirs = {}
    for tool_dir in ('rmats', 'majiq', 'isoformswitchr', 'leafcutter',
                     'sashimi', 'pegasas'):
        tool_path = os.path.join(outdir, tool_dir, comparison_id)
        if os.path.isdir(tool_path):
            raw_dirs[tool_dir] = os.path.join(tool_dir, comparison_id)\
                .replace(os.sep, '/')

    thresholds = {
        'fdr_cutoff': pipeline.get('fdr_cutoff'),
        'dpsi_cutoff': pipeline.get('dpsi_cutoff'),
        'majiq_probability_threshold': pipeline.get('majiq_probability_threshold'),
        'majiq_dpsi_cutoff': pipeline.get('majiq_dpsi_cutoff'),
        'rmats_rule': 'fdr <= fdr_cutoff & |inc_level_difference| >= dpsi_cutoff',
    }
    for key, value in thresholds.items():
        if value is None:
            thresholds[key] = ''

    group_counts = group_counts_from_index(outdir, comparison_id)

    lines = []
    lines.append('# Contrast deliverable manifest - generated by build_results_manifest.py')
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
        lines.append(yaml_block('cross_tool', list(cross_tool.items()), indent=0))
    if plots:
        lines.append('plots:')
        for area, items in plots.items():
            lines.append('  %s:' % area)
            for key, path in items.items():
                lines.append('    %s: %s' % (key, yaml_scalar(path)))
    if raw_dirs:
        lines.append('raw_outputs:')
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
        build_tools_matrix(payload.get('pipeline', {}))
    else:
        pipeline = payload.get('pipeline', {})
        build_contrast_manifest(
            payload['comparison_id'],
            payload.get('group1_name'),
            payload.get('group2_name'),
            payload.get('outdir') or pipeline.get('outdir'),
            pipeline
        )
    print('[build_results_manifest] %s manifest(s) written to %s' % (args.mode, os.getcwd()))


if __name__ == '__main__':
    sys.exit(main())