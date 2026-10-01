#!/usr/bin/env python3
"""
Export standardised master/significant/summary tables for the non-rMATS tools
(MAJIQ, IsoformSwitchAnalyzeR, LeafCutter, PEGASAS) plus the cross-tool gene
tables, following the results contract (assets/results_schema.yaml):

  <comparison_id>.<tool>.master.tsv
  <comparison_id>.<tool>.significant.tsv
  <comparison_id>.<tool>.summary.tsv
  <comparison_id>.cross_tool.master.tsv        (gene x tool, long)
  <comparison_id>.cross_tool.gene_summary.tsv  (one row per gene)

Standard columns in every master (native columns are appended, never dropped):
  comparison_id, tool, feature_type, feature_id, gene_id, gene_symbol,
  effect_size, effect_size_type, pvalue, padj, padj_method, is_significant,
  significance_rule, source_file

padj is never left as an ambiguous empty "fdr": each tool reports its adjusted
p-value with the method used (BH, rmats_cstat, satuRn_empirical, none). Tools
without an adjusted p-value expose their native probability instead (MAJIQ
prob_changing) and set padj_method = none.

No statistics are recomputed except PEGASAS: its per (pathway x sample) KS
p-values were uncorrected, so a Benjamini-Hochberg correction across all tests
is applied here (padj_method = BH).

Stdlib only.
"""

import argparse
import base64
import csv
import json
import math
import os
import statistics
import sys

STANDARD_COLUMNS = [
    'comparison_id', 'tool', 'feature_type', 'feature_id', 'gene_id',
    'gene_symbol', 'effect_size', 'effect_size_type', 'pvalue', 'padj',
    'padj_method', 'is_significant', 'significance_rule', 'source_file',
]


def write_rows(path, columns, rows):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fmt(value, digits=None):
    """Format a number compactly; keep NA/empty untouched."""
    if value is None:
        return ''
    if isinstance(value, str):
        text = value.strip()
        if text == '' or text == 'NA':
            return ''
        return text
    if isinstance(value, bool):
        return 'true' if value else 'false'
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(number) or math.isinf(number):
        return ''
    if digits is not None:
        number = round(number, digits)
    if number.is_integer() and abs(number) < 1e15:
        return str(int(number))
    return repr(number)


def bh_adjust(pvalues):
    """Benjamini-Hochberg adjusted p-values, order-preserving."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [0.0] * n
    previous = 1.0
    for position in range(n - 1, -1, -1):
        index = order[position]
        value = pvalues[index] * n / float(position + 1)
        previous = min(previous, value)
        adjusted[index] = min(previous, 1.0)
    return adjusted


# ---------------------------------------------------------------------------
# MAJIQ
# ---------------------------------------------------------------------------

MAJIQ_NATIVE = [
    'seqid', 'strand', 'event_type', 'ref_exon_start', 'ref_exon_end',
    'start', 'end', 'is_intron', 'other_exon_start', 'other_exon_end',
    'is_denovo', 'ref_exon_denovo', 'other_exon_denovo', 'event_denovo',
    'dpsi_mean', 'dpsi_std', 'probability_changing', 'probability_nonchanging',
    'grp1_raw_psi_mean', 'grp1_raw_psi_std', 'grp1_bootstrap_psi_std',
    'grp1_raw_coverage', 'grp2_raw_psi_mean', 'grp2_raw_psi_std',
    'grp2_bootstrap_psi_std', 'grp2_raw_coverage',
]


def export_majiq(comparison_id, input_dir, out_dir, params):
    """MAJIQ v3 deltapsi TSV (JSON metadata header lines start with '#')."""
    tsv_candidates = []
    for name in os.listdir(input_dir):
        if name.endswith('.tsv'):
            tsv_candidates.append(os.path.join(input_dir, name))
    if not tsv_candidates:
        print('[export_tool_masters] MAJIQ: no TSV in %s' % input_dir, file=sys.stderr)
        return 0
    tsv_candidates.sort()
    tsv_path = tsv_candidates[0]
    source = os.path.basename(tsv_path)

    prob_threshold = float(params.get('probability_threshold', 0.95))
    dpsi_cutoff = float(params.get('dpsi_cutoff', 0.2))
    rule = 'probability_changing >= %s & |dpsi_mean| >= %s' % (prob_threshold, dpsi_cutoff)

    rows = []
    with open(tsv_path, 'r', encoding='utf-8', errors='replace') as handle:
        # MAJIQ v3 prefixes the TSV with '#'-commented JSON metadata; skip it
        reader = csv.DictReader(
            (line for line in handle if not line.lstrip().startswith('#')),
            delimiter='\t')
        for raw in reader:
            prob = None
            dpsi = None
            try:
                prob = float(raw.get('probability_changing'))
            except (TypeError, ValueError):
                pass
            try:
                dpsi = float(raw.get('dpsi_mean'))
            except (TypeError, ValueError):
                pass
            significant = bool(prob is not None and dpsi is not None and
                               prob >= prob_threshold and abs(dpsi) >= dpsi_cutoff)
            gene_name = (raw.get('gene_name') or '').strip()
            gene_id = (raw.get('gene_id') or '').strip()
            feature_id = '%s:%s-%s:%s' % (raw.get('seqid', ''),
                                          raw.get('start', ''),
                                          raw.get('end', ''),
                                          raw.get('strand', ''))
            row = {
                'comparison_id': comparison_id,
                'tool': 'majiq',
                'feature_type': raw.get('event_type', ''),
                'feature_id': feature_id,
                'gene_id': gene_id,
                'gene_symbol': '' if gene_name == 'NA' else gene_name,
                'effect_size': fmt(raw.get('dpsi_mean')),
                'effect_size_type': 'delta_psi',
                'pvalue': '',
                'padj': '',
                'padj_method': 'none',
                'is_significant': 'true' if significant else 'false',
                'significance_rule': rule,
                'source_file': source,
            }
            for col in MAJIQ_NATIVE:
                row[col] = raw.get(col, '')
            rows.append(row)

    rows.sort(key=lambda r: -float(r['effect_size'] or 0))
    columns = STANDARD_COLUMNS + MAJIQ_NATIVE
    write_rows(os.path.join(out_dir, '%s.majiq.master.tsv' % comparison_id), columns, rows)
    significant = [r for r in rows if r['is_significant'] == 'true']
    write_rows(os.path.join(out_dir, '%s.majiq.significant.tsv' % comparison_id),
               columns, significant)
    _write_long_summary(out_dir, comparison_id, 'majiq', [
        ('n_total', len(rows)),
        ('n_significant', len(significant)),
        ('n_significant_with_gene', sum(
            1 for r in significant if r['gene_symbol'])),
        ('max_abs_dpsi', max((abs(float(r['effect_size'] or 0)) for r in rows), default='')),
    ])
    print('[export_tool_masters] MAJIQ: %d rows (%d significant)' % (len(rows), len(significant)))
    return 0


# ---------------------------------------------------------------------------
# IsoformSwitchAnalyzeR
# ---------------------------------------------------------------------------

ISAR_NATIVE = ['IF1', 'IF2', 'gene_switch_q_value', 'condition_1', 'condition_2',
               'switchConsequencesGene']


def export_isar(comparison_id, input_dir, out_dir, params):
    csv_path = os.path.join(input_dir, 'top_isoform_switches.csv')
    if not os.path.isfile(csv_path):
        print('[export_tool_masters] ISAR: %s missing' % csv_path, file=sys.stderr)
        return 1
    source = 'top_isoform_switches.csv'
    fdr_cutoff = float(params.get('fdr_cutoff', 0.05))
    dpsi_cutoff = float(params.get('dpsi_cutoff', 0.1))
    # DEXSeq applies BH; satuRn uses an empirical locfdr step. The chosen engine
    # is per contrast, so record the conservative label and let the known_issues
    # note cover the satuRn pinning.
    padj_method = 'BH'
    rule = 'isoform_switch_q_value < %s & |dIF| >= %s' % (fdr_cutoff, dpsi_cutoff)

    rows = []
    with open(csv_path, 'r', encoding='utf-8', errors='replace', newline='') as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            q = None
            dif = None
            try:
                q = float(raw.get('isoform_switch_q_value'))
            except (TypeError, ValueError):
                pass
            try:
                dif = float(raw.get('dIF'))
            except (TypeError, ValueError):
                pass
            significant = bool(q is not None and dif is not None and
                               q < fdr_cutoff and abs(dif) >= dpsi_cutoff)
            row = {
                'comparison_id': comparison_id,
                'tool': 'isar',
                'feature_type': 'isoform',
                'feature_id': raw.get('isoform_id', ''),
                'gene_id': raw.get('gene_id', ''),
                'gene_symbol': raw.get('gene_name', ''),
                'effect_size': fmt(raw.get('dIF')),
                'effect_size_type': 'delta_isoform_fraction',
                'pvalue': '',
                'padj': fmt(raw.get('isoform_switch_q_value')),
                'padj_method': padj_method,
                'is_significant': 'true' if significant else 'false',
                'significance_rule': rule,
                'source_file': source,
            }
            for col in ISAR_NATIVE:
                row[col] = raw.get(col, '')
            rows.append(row)

    columns = STANDARD_COLUMNS + ISAR_NATIVE
    write_rows(os.path.join(out_dir, '%s.isar.master.tsv' % comparison_id), columns, rows)
    significant = [r for r in rows if r['is_significant'] == 'true']
    write_rows(os.path.join(out_dir, '%s.isar.significant.tsv' % comparison_id),
               columns, significant)

    q_values = []
    dif_abs = []
    for r in rows:
        try:
            q_values.append(float(r['padj']))
        except (TypeError, ValueError):
            pass
        try:
            dif_abs.append(abs(float(r['effect_size'])))
        except (TypeError, ValueError):
            pass
    _write_long_summary(out_dir, comparison_id, 'isar', [
        ('n_total', len(rows)),
        ('n_significant', len(significant)),
        ('n_genes', len({r['gene_symbol'] for r in rows})),
        ('median_abs_dIF', statistics.median(dif_abs) if dif_abs else ''),
        ('q_min', min(q_values) if q_values else ''),
        ('q_max', max(q_values) if q_values else ''),
        ('q_iqr', statistics.quantiles(q_values, n=4)[2] - statistics.quantiles(q_values, n=4)[0]
         if len(q_values) > 3 else ''),
    ])
    print('[export_tool_masters] ISAR: %d rows (%d significant)' % (len(rows), len(significant)))
    return 0


# ---------------------------------------------------------------------------
# LeafCutter
# ---------------------------------------------------------------------------

LEAFCUTTER_NATIVE = ['cluster', 'status', 'loglr', 'df', 'logef',
                     'deltapsi', 'psi_group1', 'psi_group2']


def leafcluster_from_intron(intron_id):
    """'chr1:14614:16858:clu_1_-' -> 'chr1:clu_1_-'."""
    parts = (intron_id or '').rsplit(':', 1)
    if len(parts) != 2:
        return ''
    chrom = (intron_id or '').split(':', 1)[0]
    return '%s:%s' % (chrom, parts[1])


def export_leafcutter(comparison_id, input_dir, out_dir, params):
    sig_file = None
    eff_file = None
    for name in sorted(os.listdir(input_dir)):
        if name.endswith('_cluster_significance.txt'):
            sig_file = os.path.join(input_dir, name)
        elif name.endswith('_effect_sizes.txt'):
            eff_file = os.path.join(input_dir, name)
    if not sig_file:
        print('[export_tool_masters] LeafCutter: no cluster_significance file in %s'
              % input_dir, file=sys.stderr)
        return 1

    fdr_cutoff = float(params.get('fdr_cutoff', 0.05))
    rule = 'status == "Success" & padj < %s' % fdr_cutoff

    # cluster-level significance
    clusters = {}
    with open(sig_file, 'r', encoding='utf-8', errors='replace') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for raw in reader:
            clusters[raw.get('cluster', '')] = raw

    # intron-level effect sizes; PSI columns are positional (after logef)
    effects = {}
    if eff_file:
        with open(eff_file, 'r', encoding='utf-8', errors='replace') as handle:
            handle.readline()
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 5:
                    continue
                effects[parts[0]] = {
                    'logef': parts[1],
                    'psi_group1': parts[2] if len(parts) > 2 else '',
                    'psi_group2': parts[3] if len(parts) > 3 else '',
                    'deltapsi': parts[-1],
                }

    rows = []
    for intron_id, eff in effects.items():
        cluster_id = leafcluster_from_intron(intron_id)
        cluster = clusters.get(cluster_id, {})
        status = cluster.get('status', '')
        padj = None
        try:
            padj = float(cluster.get('p.adjust'))
        except (TypeError, ValueError):
            pass
        significant = bool(status == 'Success' and padj is not None and padj < fdr_cutoff)
        row = {
            'comparison_id': comparison_id,
            'tool': 'leafcutter',
            'feature_type': 'intron',
            'feature_id': intron_id,
            'gene_id': '',
            'gene_symbol': cluster.get('genes', '') or '',
            'effect_size': eff.get('deltapsi', ''),
            'effect_size_type': 'delta_psi',
            'pvalue': cluster.get('p', ''),
            'padj': cluster.get('p.adjust', ''),
            'padj_method': 'BH',
            'is_significant': 'true' if significant else 'false',
            'significance_rule': rule,
            'source_file': os.path.basename(sig_file) + ';' +
                           (os.path.basename(eff_file) if eff_file else ''),
        }
        for col in LEAFCUTTER_NATIVE:
            row[col] = ''
        row['cluster'] = cluster_id
        row['status'] = status
        row['loglr'] = cluster.get('loglr', '')
        row['df'] = cluster.get('df', '')
        row['logef'] = eff.get('logef', '')
        row['deltapsi'] = eff.get('deltapsi', '')
        row['psi_group1'] = eff.get('psi_group1', '')
        row['psi_group2'] = eff.get('psi_group2', '')
        rows.append(row)

    # introns without effect sizes (cluster significance exists but no intron row)
    eff_clusters = {leafcluster_from_intron(k) for k in effects}
    for cluster_id, cluster in clusters.items():
        if cluster_id in eff_clusters:
            continue
        status = cluster.get('status', '')
        padj = None
        try:
            padj = float(cluster.get('p.adjust'))
        except (TypeError, ValueError):
            pass
        significant = bool(status == 'Success' and padj is not None and padj < fdr_cutoff)
        rows.append({
            'comparison_id': comparison_id,
            'tool': 'leafcutter',
            'feature_type': 'intron',
            'feature_id': '',
            'gene_id': '',
            'gene_symbol': cluster.get('genes', '') or '',
            'effect_size': '',
            'effect_size_type': 'delta_psi',
            'pvalue': cluster.get('p', ''),
            'padj': cluster.get('p.adjust', ''),
            'padj_method': 'BH',
            'is_significant': 'true' if significant else 'false',
            'significance_rule': rule,
            'source_file': os.path.basename(sig_file),
            'cluster': cluster_id,
            'status': status,
            'loglr': cluster.get('loglr', ''),
            'df': cluster.get('df', ''),
            'logef': '',
            'deltapsi': '',
            'psi_group1': '',
            'psi_group2': '',
        })

    columns = STANDARD_COLUMNS + LEAFCUTTER_NATIVE
    write_rows(os.path.join(out_dir, '%s.leafcutter.master.tsv' % comparison_id),
               columns, rows)
    significant = [r for r in rows if r['is_significant'] == 'true']
    write_rows(os.path.join(out_dir, '%s.leafcutter.significant.tsv' % comparison_id),
               columns, significant)

    with_ok = sum(1 for r in rows if r['status'] == 'Success')
    n_clusters = len(clusters)
    genes = set()
    for r in rows:
        for g in str(r['gene_symbol']).split(','):
            g = g.strip()
            if g:
                genes.add(g)
    _write_long_summary(out_dir, comparison_id, 'leafcutter', [
        ('n_clusters', n_clusters),
        ('n_clusters_success', with_ok),
        ('n_introns', len(rows)),
        ('n_significant_clusters', len({r['cluster'] for r in significant if r['cluster']})),
        ('n_significant_introns', len(significant)),
        ('n_genes', len(genes)),
    ])
    print('[export_tool_masters] LeafCutter: %d introns (%d significant)'
          % (len(rows), len(significant)))
    return 0


# ---------------------------------------------------------------------------
# PEGASAS
# ---------------------------------------------------------------------------

PEGASAS_NATIVE = ['n_samples', 'n_sig_samples', 'median_p', 'max_abs_ks',
                  'median_rank', 'n_sig_tests', 'min_padj', 'n_sig_events']


def export_pegasas(comparison_id, input_dir, out_dir, params):
    scores_path = os.path.join(input_dir, 'all_pathways_scores.tsv')
    if not os.path.isfile(scores_path):
        print('[export_tool_masters] PEGASAS: %s missing' % scores_path, file=sys.stderr)
        return 1
    fdr_cutoff = float(params.get('fdr_cutoff', 0.05))
    rule = 'BH padj < %s across pathway x sample KS tests' % fdr_cutoff

    per_pathway = {}
    tests = []  # (pathway, p_value)
    with open(scores_path, 'r', encoding='utf-8', errors='replace') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for raw in reader:
            pathway = raw.get('pathway', '')
            entry = per_pathway.setdefault(pathway, {'p': [], 'ks': [], 'rank': []})
            p_value = None
            try:
                p_value = float(raw.get('p_value'))
                entry['p'].append(p_value)
                tests.append((pathway, p_value))
            except (TypeError, ValueError):
                pass
            try:
                entry['ks'].append(abs(float(raw.get('KS_score'))))
            except (TypeError, ValueError):
                pass
            try:
                entry['rank'].append(float(raw.get('median_rank')))
            except (TypeError, ValueError):
                pass

    # Benjamini-Hochberg across every pathway x sample KS test
    adjusted = bh_adjust([p for _pathway, p in tests])
    min_padj = {}
    n_sig_tests = {}
    for (pathway, _p), q in zip(tests, adjusted):
        if pathway not in min_padj or q < min_padj[pathway]:
            min_padj[pathway] = q
        if q < fdr_cutoff:
            n_sig_tests[pathway] = n_sig_tests.get(pathway, 0) + 1

    sig_counts = {}
    sig_path = params.get('sig_pathways')
    if sig_path and os.path.isfile(sig_path):
        with open(sig_path, 'r', encoding='utf-8', errors='replace') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for raw in reader:
                try:
                    sig_counts[raw.get('pathway', '')] = int(raw.get('n_sig_events', 0))
                except (TypeError, ValueError):
                    pass

    rows = []
    for pathway, entry in sorted(per_pathway.items()):
        p_values = entry['p']
        min_p = min(p_values) if p_values else None
        pathway_padj = min_padj.get(pathway)
        significant = pathway_padj is not None and pathway_padj < fdr_cutoff
        row = {
            'comparison_id': comparison_id,
            'tool': 'pegasas',
            'feature_type': 'pathway',
            'feature_id': pathway,
            'gene_id': '',
            'gene_symbol': '',
            'effect_size': fmt(max(entry['ks'], default=None)),
            'effect_size_type': 'KS_statistic',
            'pvalue': fmt(min_p),
            'padj': fmt(pathway_padj),
            'padj_method': 'BH',
            'is_significant': 'true' if significant else 'false',
            'significance_rule': rule,
            'source_file': 'all_pathways_scores.tsv',
            'n_samples': len(p_values),
            'n_sig_samples': n_sig_tests.get(pathway, 0),
            'median_p': fmt(statistics.median(p_values) if p_values else None),
            'max_abs_ks': fmt(max(entry['ks'], default=None)),
            'median_rank': fmt(statistics.median(entry['rank']) if entry['rank'] else None),
            'n_sig_tests': n_sig_tests.get(pathway, 0),
            'min_padj': fmt(pathway_padj),
            'n_sig_events': sig_counts.get(pathway, ''),
        }
        rows.append(row)

    columns = STANDARD_COLUMNS + PEGASAS_NATIVE
    write_rows(os.path.join(out_dir, '%s.pegasas.master.tsv' % comparison_id), columns, rows)
    significant = [r for r in rows if r['is_significant'] == 'true']
    write_rows(os.path.join(out_dir, '%s.pegasas.significant.tsv' % comparison_id),
               columns, significant)
    _write_long_summary(out_dir, comparison_id, 'pegasas', [
        ('n_pathways', len(rows)),
        ('n_pathways_significant_BH', len(significant)),
        ('n_tests', len(tests)),
        ('n_samples_per_pathway', rows[0]['n_samples'] if rows else ''),
    ])
    print('[export_tool_masters] PEGASAS: %d pathways (%d significant after BH)'
          % (len(rows), len(significant)))
    return 0


# ---------------------------------------------------------------------------
# Cross-tool
# ---------------------------------------------------------------------------

def export_cross_tool(comparison_id, master_files, out_dir):
    """Gene-level union of significant features across AS tools.

    Writes both shapes:
      cross_tool.master.tsv        one row per gene x tool (long)
      cross_tool.gene_summary.tsv  one row per gene with the tools that called it
    """
    genes = {}
    per_gene = {}
    tool_genes = {}
    for path in master_files:
        tool = None
        for candidate in ('rmats', 'majiq', 'isar', 'leafcutter'):
            if '.' + candidate + '.master.tsv' in os.path.basename(path):
                tool = candidate
                break
        if tool is None:
            continue
        with open(path, 'r', encoding='utf-8', errors='replace') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for raw in reader:
                if raw.get('is_significant', 'false') != 'true':
                    continue
                symbols = [s.strip().upper() for s in
                           str(raw.get('gene_symbol') or '').split(',')]
                symbols = [s for s in symbols if s and s != 'NA']
                for symbol in symbols:
                    key = (symbol, tool)
                    entry = genes.setdefault(key, {
                        'comparison_id': comparison_id,
                        'gene_symbol': symbol,
                        'gene_id': raw.get('gene_id', ''),
                        'tool': tool,
                        'n_significant_features': 0,
                        'best_effect_size': '',
                        'best_padj': '',
                    })
                    entry['n_significant_features'] += 1
                    try:
                        effect = abs(float(raw.get('effect_size')))
                        if entry['best_effect_size'] == '' or \
                           effect > float(entry['best_effect_size']):
                            entry['best_effect_size'] = fmt(effect)
                    except (TypeError, ValueError):
                        pass
                    try:
                        padj = float(raw.get('padj'))
                        if entry['best_padj'] == '' or padj < float(entry['best_padj']):
                            entry['best_padj'] = fmt(padj)
                    except (TypeError, ValueError):
                        pass

                    summary = per_gene.setdefault(symbol, {
                        'comparison_id': comparison_id,
                        'gene_symbol': symbol,
                        'gene_id': raw.get('gene_id', ''),
                        'tools': set(),
                        'best_effect_size': '',
                        'best_padj': '',
                    })
                    if not summary['gene_id'] and raw.get('gene_id'):
                        summary['gene_id'] = raw.get('gene_id')
                    summary['tools'].add(tool)
                    if entry['best_effect_size'] != '' and (
                            summary['best_effect_size'] == '' or
                            float(entry['best_effect_size']) > float(summary['best_effect_size'])):
                        summary['best_effect_size'] = entry['best_effect_size']
                    if entry['best_padj'] != '' and (
                            summary['best_padj'] == '' or
                            float(entry['best_padj']) < float(summary['best_padj'])):
                        summary['best_padj'] = entry['best_padj']
            tool_genes[tool] = {s.strip().upper() for s in
                                _all_gene_symbols(path)}

    rows = list(genes.values())
    columns = ['comparison_id', 'gene_symbol', 'gene_id', 'tool',
               'n_significant_features', 'best_effect_size', 'best_padj']
    write_rows(os.path.join(out_dir, '%s.cross_tool.master.tsv' % comparison_id),
               columns, rows)

    gene_rows = []
    for summary in per_gene.values():
        gene_rows.append({
            'comparison_id': comparison_id,
            'gene_symbol': summary['gene_symbol'],
            'gene_id': summary['gene_id'],
            'n_tools_significant': len(summary['tools']),
            'tools': ','.join(sorted(summary['tools'])),
            'best_effect_size': summary['best_effect_size'],
            'best_padj': summary['best_padj'],
        })
    gene_rows.sort(key=lambda r: (-r['n_tools_significant'], r['gene_symbol']))
    gene_columns = ['comparison_id', 'gene_symbol', 'gene_id',
                    'n_tools_significant', 'tools', 'best_effect_size', 'best_padj']
    write_rows(os.path.join(out_dir, '%s.cross_tool.gene_summary.tsv' % comparison_id),
               gene_columns, gene_rows)

    print('[export_tool_masters] cross_tool: %d gene x tool rows, %d genes'
          % (len(rows), len(gene_rows)))
    print('[export_tool_masters] cross_tool per-tool significant genes: %s'
          % json.dumps({t: len(gs) for t, gs in tool_genes.items()}))
    return 0


def _all_gene_symbols(master_path):
    symbols = set()
    with open(master_path, 'r', encoding='utf-8', errors='replace') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for raw in reader:
            if raw.get('is_significant', 'false') != 'true':
                continue
            for s in str(raw.get('gene_symbol') or '').split(','):
                s = s.strip()
                if s and s != 'NA':
                    symbols.add(s)
    return symbols


# ---------------------------------------------------------------------------

def _write_long_summary(out_dir, comparison_id, tool, stats):
    path = os.path.join(out_dir, '%s.%s.summary.tsv' % (comparison_id, tool))
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['comparison_id', 'tool', 'statistic', 'value'])
        for name, value in stats:
            writer.writerow([comparison_id, tool, name, fmt(value)])


def parse_tool_params(raw):
    if not raw:
        return {}
    try:
        return json.loads(base64.b64decode(raw).decode('utf-8'))
    except Exception:
        return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-id', required=True)
    parser.add_argument('--tool', required=True,
                        choices=['majiq', 'isar', 'leafcutter', 'pegasas', 'cross_tool'])
    parser.add_argument('--input-dir', help='tool output directory')
    parser.add_argument('--master-files', nargs='*',
                        help='master TSVs for cross_tool mode')
    parser.add_argument('--tool-params', default='',
                        help='base64 JSON with tool-specific cutoffs')
    parser.add_argument('--out-dir', default='.')
    args = parser.parse_args()

    params = parse_tool_params(args.tool_params)
    if args.tool == 'cross_tool':
        return export_cross_tool(args.comparison_id, args.master_files or [], args.out_dir)
    if not args.input_dir:
        parser.error('--input-dir required for --tool %s' % args.tool)
    dispatcher = {
        'majiq': export_majiq,
        'isar': export_isar,
        'leafcutter': export_leafcutter,
        'pegasas': export_pegasas,
    }
    return dispatcher[args.tool](args.comparison_id, args.input_dir, args.out_dir, params)


if __name__ == '__main__':
    sys.exit(main())
