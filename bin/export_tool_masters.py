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
import re
import statistics
import sys

STANDARD_COLUMNS = [
    'comparison_id', 'tool', 'feature_type', 'feature_id', 'gene_id',
    'gene_symbol', 'group1_name', 'group2_name', 'effect_size_direction',
    'effect_size', 'effect_size_type', 'pvalue', 'padj',
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


def mannwhitney_p(group1, group2):
    """Two-sided Mann-Whitney U p-value (normal approximation, tie-averaged).

    Returns None when either group has fewer than 4 values, so small contrasts
    do not get a spuriously precise activity difference.
    """
    n1, n2 = len(group1), len(group2)
    if n1 < 4 or n2 < 4:
        return None
    combined = sorted([(v, 0) for v in group1] + [(v, 1) for v in group2])
    n = n1 + n2
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and combined[j + 1][0] == combined[i][0]:
            j += 1
        average = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = average
        i = j + 1
    rank_sum1 = sum(ranks[k] for k in range(n) if combined[k][1] == 0)
    u1 = rank_sum1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sigma == 0:
        return None
    z = abs(u1 - mu) / sigma
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))))


# ---------------------------------------------------------------------------
# MAJIQ
# ---------------------------------------------------------------------------

MAJIQ_NATIVE = [
    'lsv_id',
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
    group1_name = params.get('group1_name', '')
    group2_name = params.get('group2_name', '')
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
            # An LSV is the reference exon context; a row is one junction within
            # it. Keep lsv_id for LSV-level aggregation and make feature_id the
            # unique per-junction identity (the old seqid:start-end:strand was
            # duplicated across rows of the same LSV).
            ref_exon = '%s-%s' % (raw.get('ref_exon_start', ''), raw.get('ref_exon_end', ''))
            lsv_id = '%s:%s:%s' % (raw.get('seqid', ''), raw.get('strand', ''), ref_exon)
            feature_id = '%s:%s:%s-%s:%s-%s:%s' % (
                raw.get('seqid', ''), raw.get('strand', ''),
                raw.get('start', ''), raw.get('end', ''),
                raw.get('other_exon_start', ''), raw.get('other_exon_end', ''),
                raw.get('is_intron', ''))
            row = {
                'comparison_id': comparison_id,
                'tool': 'majiq',
                'feature_type': raw.get('event_type', ''),
                'feature_id': feature_id,
                'gene_id': gene_id,
                'gene_symbol': '' if gene_name == 'NA' else gene_name,
                'group1_name': group1_name,
                'group2_name': group2_name,
                'effect_size_direction': 'group2_minus_group1',
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
            row['lsv_id'] = lsv_id
            rows.append(row)

    rows.sort(key=lambda r: -float(r['effect_size'] or 0))
    columns = STANDARD_COLUMNS + MAJIQ_NATIVE
    write_rows(os.path.join(out_dir, '%s.majiq.master.tsv' % comparison_id), columns, rows)
    significant = [r for r in rows if r['is_significant'] == 'true']
    write_rows(os.path.join(out_dir, '%s.majiq.significant.tsv' % comparison_id),
               columns, significant)
    _write_long_summary(out_dir, comparison_id, 'majiq', [
        ('n_total', len(rows)),
        ('n_lsvs', len({r['lsv_id'] for r in rows})),
        ('n_significant', len(significant)),
        ('n_significant_lsvs', len({r['lsv_id'] for r in significant})),
        ('n_significant_with_gene', sum(
            1 for r in significant if r['gene_symbol'])),
        ('max_abs_dpsi', max((abs(float(r['effect_size'] or 0)) for r in rows), default='')),
    ])
    print('[export_tool_masters] MAJIQ: %d rows (%d significant)' % (len(rows), len(significant)))
    return 0


# ---------------------------------------------------------------------------
# IsoformSwitchAnalyzeR
# ---------------------------------------------------------------------------

def strip_version(identifier):
    return str(identifier or '').split('.')[0]


def load_transcript_gene_map(gtf_path):
    """ENST (unversioned) -> ENSG from a GTF. Empty dict when unavailable.

    Used to recover the Ensembl gene id for ISAR rows, whose CSV carries the
    gene symbol in gene_id and no ref_gene_id.
    """
    mapping = {}
    if not gtf_path or not os.path.isfile(gtf_path):
        return mapping
    pattern = re.compile(r'gene_id "([^"]+)".*?transcript_id "([^"]+)"')
    with open(gtf_path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            match = pattern.search(fields[8])
            if not match:
                continue
            gene_id = strip_version(match.group(1))
            transcript_id = strip_version(match.group(2))
            if gene_id.startswith('ENSG') and transcript_id:
                mapping[transcript_id] = gene_id
    return mapping


ISAR_NATIVE = ['IF1', 'IF2', 'gene_switch_q_value', 'switchConsequencesGene',
               'effective_test_method']


def export_isar(comparison_id, input_dir, out_dir, params):
    csv_path = os.path.join(input_dir, 'top_isoform_switches.csv')
    if not os.path.isfile(csv_path):
        print('[export_tool_masters] ISAR: %s missing' % csv_path, file=sys.stderr)
        return 1
    source = 'top_isoform_switches.csv'
    fdr_cutoff = float(params.get('fdr_cutoff', 0.05))
    dpsi_cutoff = float(params.get('dpsi_cutoff', 0.1))
    group1_name = params.get('group1_name', '')
    group2_name = params.get('group2_name', '')
    effective_method = str(params.get('effective_test_method', '') or '').strip()
    # DEXSeq applies Benjamini-Hochberg; satuRn uses an empirical locfdr step.
    # Derive the label from the engine actually used for this contrast instead
    # of hardcoding BH.
    if effective_method.lower().startswith('satur'):
        padj_method = 'empirical_FDR'
    elif effective_method.lower() == 'dexseq':
        padj_method = 'BH'
    else:
        padj_method = 'BH'
    rule = 'isoform_switch_q_value <= %s & |dIF| >= %s' % (fdr_cutoff, dpsi_cutoff)
    ensg_map = load_transcript_gene_map(params.get('gtf'))

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
                               q <= fdr_cutoff and abs(dif) >= dpsi_cutoff)
            # The ISAR CSV carries the gene symbol in gene_id and has no
            # ref_gene_id. Recover the Ensembl gene id from the transcript id via
            # the pipeline GTF; leave it empty for unannotated transcripts.
            transcript_id = strip_version(raw.get('isoform_id'))
            gene_id = ensg_map.get(transcript_id, '')
            consequence = str(raw.get('switchConsequencesGene') or '').strip()
            if consequence in ('NA', 'None'):
                consequence = ''
            row = {
                'comparison_id': comparison_id,
                'tool': 'isar',
                'feature_type': 'isoform',
                'feature_id': raw.get('isoform_id', ''),
                'gene_id': gene_id,
                'gene_symbol': str(raw.get('gene_name') or '').strip(),
                'group1_name': group1_name,
                'group2_name': group2_name,
                'effect_size_direction': 'group2_minus_group1',
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
                if col == 'effective_test_method':
                    row[col] = effective_method
                elif col == 'switchConsequencesGene':
                    row[col] = consequence
                else:
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

LEAFCUTTER_NATIVE = ['cluster', 'status', 'cluster_significant', 'loglr', 'df',
                     'logef', 'deltapsi', 'psi_group1', 'psi_group2']


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
    dpsi_cutoff = float(params.get('dpsi_cutoff', 0.1))
    group1_name = params.get('group1_name', '')
    group2_name = params.get('group2_name', '')
    rule = ('status == "Success" & padj_cluster <= %s & '
            '|deltapsi| >= %s') % (fdr_cutoff, dpsi_cutoff)

    # cluster-level significance
    clusters = {}
    with open(sig_file, 'r', encoding='utf-8', errors='replace') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for raw in reader:
            clusters[raw.get('cluster', '')] = raw

    # intron-level effect sizes. PSI columns are named with the two group names
    # by leafcutter_ds.R (columns 3 and 4), in groups_file order; resolve them by
    # header name so the sign and group labels cannot drift with column order.
    effects = {}
    eff_group_names = []
    if eff_file:
        with open(eff_file, 'r', encoding='utf-8', errors='replace') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            eff_group_names = header[2:-1]  # drop intron, logef and deltapsi
            if group1_name and eff_group_names and eff_group_names[0] != group1_name:
                print('[export_tool_masters] LeafCutter: first PSI column is %r but '
                      'group1 is %r — cluster column order is not group1-first '
                      '(%s)' % (eff_group_names[0], group1_name, eff_file),
                      file=sys.stderr)
                return 1
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 5:
                    continue
                record = dict(zip(header, parts))
                effects[parts[0]] = {
                    'logef': record.get('logef', parts[1]),
                    'psi_group1': record.get(group1_name, '') if group1_name else '',
                    'psi_group2': record.get(group2_name, '') if group2_name else '',
                    'deltapsi': record.get('deltapsi', parts[-1]),
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
        cluster_significant = bool(status == 'Success' and padj is not None and padj <= fdr_cutoff)
        dpsi = None
        try:
            dpsi = float(eff.get('deltapsi'))
        except (TypeError, ValueError):
            pass
        # Significance is called per intron, not per cluster: the cluster must
        # pass FDR and the intron must clear the |deltaPSI| cutoff.
        significant = bool(cluster_significant and dpsi is not None and
                           abs(dpsi) >= dpsi_cutoff)
        row = {
            'comparison_id': comparison_id,
            'tool': 'leafcutter',
            'feature_type': 'intron',
            'feature_id': intron_id,
            'gene_id': '',
            'gene_symbol': cluster.get('genes', '') or '',
            'group1_name': group1_name,
            'group2_name': group2_name,
            'effect_size_direction': 'group2_minus_group1',
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
        row['cluster_significant'] = 'true' if cluster_significant else 'false'
        row['loglr'] = cluster.get('loglr', '')
        row['df'] = cluster.get('df', '')
        row['logef'] = eff.get('logef', '')
        row['deltapsi'] = eff.get('deltapsi', '')
        row['psi_group1'] = eff.get('psi_group1', '')
        row['psi_group2'] = eff.get('psi_group2', '')
        rows.append(row)

    columns = STANDARD_COLUMNS + LEAFCUTTER_NATIVE
    write_rows(os.path.join(out_dir, '%s.leafcutter.master.tsv' % comparison_id),
               columns, rows)
    significant = [r for r in rows if r['is_significant'] == 'true']
    write_rows(os.path.join(out_dir, '%s.leafcutter.significant.tsv' % comparison_id),
               columns, significant)

    n_clusters = len(clusters)
    n_clusters_success = len({r['cluster'] for r in rows if r['status'] == 'Success'})
    n_clusters_untested = n_clusters - len({r['cluster'] for r in rows})
    genes = set()
    for r in rows:
        for g in str(r['gene_symbol']).split(','):
            g = g.strip()
            if g:
                genes.add(g)
    _write_long_summary(out_dir, comparison_id, 'leafcutter', [
        ('n_clusters', n_clusters),
        ('n_clusters_success', n_clusters_success),
        ('n_clusters_untested', n_clusters_untested),
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

PEGASAS_NATIVE = ['n_samples', 'n_group1', 'n_group2',
                  'activity_diff', 'activity_diff_pvalue',
                  'ks_score_max', 'ks_score_median', 'ks_min_pvalue',
                  'n_sig_events', 'n_sig_events_perm', 'n_tested', 'min_padj_bh']


def export_pegasas(comparison_id, input_dir, out_dir, params):
    scores_path = os.path.join(input_dir, 'all_pathways_scores.tsv')
    if not os.path.isfile(scores_path):
        print('[export_tool_masters] PEGASAS: %s missing' % scores_path, file=sys.stderr)
        return 1
    group1_name = params.get('group1_name', '')
    group2_name = params.get('group2_name', '')
    # PEGASAS has no boolean significance call. n_sig_events is the BH-adjusted
    # count of event x pathway correlations using the analytic Pearson p-value
    # (BH applied once across the contrast matrix); the permutation p-value has
    # a 1/N_PERMS floor, so its BH count is kept only as n_sig_events_perm.
    rule = ('descriptive: n_sig_events = BH-adjusted event x pathway correlations '
            '(analytic Pearson p, FDR < %s across the contrast matrix); '
            'n_sig_events_perm is the permutation cross-check; no boolean call'
            % params.get('fdr_cutoff', 0.05))

    per_pathway = {}
    with open(scores_path, 'r', encoding='utf-8', errors='replace') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for raw in reader:
            pathway = raw.get('pathway', '')
            entry = per_pathway.setdefault(
                pathway, {'g1': [], 'g2': [], 'all': [], 'p': []})
            group = str(raw.get('group', ''))
            try:
                p_value = float(raw.get('p_value'))
                entry['p'].append(p_value)
            except (TypeError, ValueError):
                pass
            try:
                ks = float(raw.get('KS_score'))
            except (TypeError, ValueError):
                continue
            entry['all'].append(ks)
            if group == group2_name:
                entry['g2'].append(ks)
            elif group == group1_name:
                entry['g1'].append(ks)

    # Group names come from the comparison; activity_diff is undefined when a
    # group has no samples (and its p-value is suppressed for n < 4).
    sig_counts = {}
    sig_path = params.get('sig_pathways')
    if sig_path and os.path.isfile(sig_path):
        with open(sig_path, 'r', encoding='utf-8', errors='replace') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for raw in reader:
                entry = {'n_sig_events': 0, 'n_sig_events_perm': '', 'n_tested': '', 'min_padj': ''}
                try:
                    entry['n_sig_events'] = int(raw.get('n_sig_events', 0))
                except (TypeError, ValueError):
                    pass
                for key in ('n_sig_events_perm', 'n_tested', 'min_padj'):
                    if raw.get(key) not in (None, ''):
                        entry[key] = raw.get(key)
                sig_counts[raw.get('pathway', '')] = entry

    rows = []
    for pathway, entry in sorted(per_pathway.items()):
        g1 = entry['g1']
        g2 = entry['g2']
        all_ks = entry['all']
        activity_diff = None
        if g1 and g2:
            activity_diff = (sum(g2) / len(g2)) - (sum(g1) / len(g1))
        activity_p = mannwhitney_p(g1, g2) if (g1 and g2) else None
        sig_entry = sig_counts.get(pathway, {'n_sig_events': 0, 'n_sig_events_perm': '',
                                             'n_tested': '', 'min_padj': ''})
        n_sig_events = sig_entry['n_sig_events']
        row = {
            'comparison_id': comparison_id,
            'tool': 'pegasas',
            'feature_type': 'pathway',
            'feature_id': pathway,
            'gene_id': '',
            'gene_symbol': '',
            'group1_name': group1_name,
            'group2_name': group2_name,
            'effect_size_direction': 'group2_minus_group1',
            # Descriptive activity difference (mean KS group2 - group1); the KS
            # score itself is never presented as significance.
            'effect_size': fmt(activity_diff),
            'effect_size_type': 'delta_KS',
            'pvalue': '',
            'padj': '',
            'padj_method': 'BH_event_wise',
            # PEGASAS is exploratory: no row is called significant. Rank the
            # pathways by n_sig_events instead.
            'is_significant': 'false',
            'significance_rule': rule,
            'source_file': 'all_pathways_scores.tsv',
            'n_samples': len(all_ks),
            'n_group1': len(g1),
            'n_group2': len(g2),
            'activity_diff': fmt(activity_diff),
            'activity_diff_pvalue': fmt(activity_p),
            'ks_score_max': fmt(max(all_ks) if all_ks else None),
            'ks_score_median': fmt(statistics.median(all_ks) if all_ks else None),
            'ks_min_pvalue': fmt(min(entry['p']) if entry['p'] else None),
            'n_sig_events': n_sig_events,
            'n_sig_events_perm': sig_entry['n_sig_events_perm'],
            'n_tested': sig_entry['n_tested'],
            'min_padj_bh': sig_entry['min_padj'],
        }
        rows.append(row)

    # Rank by the BH-significant event count (the only per-pathway quantity).
    rows.sort(key=lambda r: (-(r['n_sig_events'] if isinstance(r['n_sig_events'], int) else 0),
                             r['feature_id']))
    columns = STANDARD_COLUMNS + PEGASAS_NATIVE
    write_rows(os.path.join(out_dir, '%s.pegasas.master.tsv' % comparison_id), columns, rows)
    # No boolean significance call: significant.tsv is intentionally empty.
    write_rows(os.path.join(out_dir, '%s.pegasas.significant.tsv' % comparison_id),
               columns, [])
    _write_long_summary(out_dir, comparison_id, 'pegasas', [
        ('n_pathways', len(rows)),
        ('n_sig_events_total', sum(r['n_sig_events'] for r in rows
                                   if isinstance(r['n_sig_events'], int))),
        ('top_pathway', rows[0]['feature_id'] if rows else ''),
        ('top_n_sig_events', rows[0]['n_sig_events'] if rows else ''),
        ('n_samples_per_pathway', rows[0]['n_samples'] if rows else ''),
    ])
    print('[export_tool_masters] PEGASAS: %d pathways (ranked by n_sig_events; '
          'no boolean significance)' % len(rows))
    return 0


# ---------------------------------------------------------------------------
# Cross-tool
# ---------------------------------------------------------------------------

def export_cross_tool(comparison_id, master_files, out_dir, params=None):
    """Gene-level union of significant features across AS tools.

    Writes both shapes:
      cross_tool.master.tsv        one row per gene x tool (long)
      cross_tool.gene_summary.tsv  one row per gene with the tools that called it
    """
    params = params or {}
    group1_name = params.get('group1_name', '')
    group2_name = params.get('group2_name', '')

    def detect_tool(path):
        for candidate in ('rmats', 'majiq', 'isar', 'leafcutter'):
            if '.' + candidate + '.master.tsv' in os.path.basename(path):
                return candidate
        return None

    def is_ensg(value):
        return bool(re.match(r'^ENSG\d+', str(value or '').strip()))

    # symbol -> ENSG map built from every row of every master. ISAR/MAJIQ/rMATS
    # carry a real gene_id; LeafCutter only has comma-separated symbols, so the
    # map lets its clusters join the ENSG-keyed overlap instead of falling back
    # to symbol == gene_id.
    symbol_to_ensg = {}
    # Single pass per master: build the symbol->ENSG map from every row and
    # collect significant rows. The map must be complete before LeafCutter
    # (symbol-only) rows are resolved, so processing is deferred.
    significant_rows = []  # (tool, row)
    for path in master_files:
        tool = detect_tool(path)
        if tool is None:
            continue
        with open(path, 'r', encoding='utf-8', errors='replace') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for raw in reader:
                gid = str(raw.get('gene_id') or '').strip()
                if is_ensg(gid):
                    sym = str(raw.get('gene_symbol') or '').split(',')[0].strip().upper()
                    if sym and sym != 'NA':
                        symbol_to_ensg.setdefault(sym, gid)
                if raw.get('is_significant', 'false') == 'true':
                    significant_rows.append((tool, raw))

    genes = {}
    per_gene = {}
    tool_stats = {}

    for tool, raw in significant_rows:
        stats = tool_stats.setdefault(tool, {'significant': set(), 'ok': set()})
        for raw in (raw,):
                gid = str(raw.get('gene_id') or '').strip()
                symbols = [s.strip().upper() for s in
                           str(raw.get('gene_symbol') or '').split(',')]
                symbols = [s for s in symbols if s and s != 'NA']
                reads_ok = raw.get('reads_ok', 'true')
                row_is_ok = reads_ok != 'false'
                for symbol in symbols:
                    if is_ensg(gid):
                        ensg = gid
                    else:
                        ensg = symbol_to_ensg.get(symbol, '')
                    # Key the overlap by ENSG; fall back to the symbol only when
                    # no ENSG could be resolved.
                    gene_key = ensg or symbol
                    key = (gene_key, tool)
                    entry = genes.setdefault(key, {
                        'comparison_id': comparison_id,
                        'gene_symbol': symbol,
                        'gene_id': ensg,
                        'tool': tool,
                        'group1_name': group1_name,
                        'group2_name': group2_name,
                        'effect_size_direction': 'group2_minus_group1',
                        'n_significant_features': 0,
                        'n_ok_features': 0,
                        'best_effect_size': '',
                        'best_padj': '',
                    })
                    if not entry['gene_symbol'] and symbol:
                        entry['gene_symbol'] = symbol
                    if not entry['gene_id'] and ensg:
                        entry['gene_id'] = ensg
                    entry['n_significant_features'] += 1
                    if row_is_ok:
                        entry['n_ok_features'] += 1
                    stats['significant'].add(gene_key)
                    if row_is_ok:
                        stats['ok'].add(gene_key)
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

                    summary = per_gene.setdefault(gene_key, {
                        'comparison_id': comparison_id,
                        'gene_symbol': symbol,
                        'gene_id': ensg,
                        'tools': set(),
                        'best_effect_size': '',
                        'best_padj': '',
                    })
                    if not summary['gene_symbol'] and symbol:
                        summary['gene_symbol'] = symbol
                    if not summary['gene_id'] and ensg:
                        summary['gene_id'] = ensg
                    summary['tools'].add(tool)
                    if entry['best_effect_size'] != '' and (
                            summary['best_effect_size'] == '' or
                            float(entry['best_effect_size']) > float(summary['best_effect_size'])):
                        summary['best_effect_size'] = entry['best_effect_size']
                    if entry['best_padj'] != '' and (
                            summary['best_padj'] == '' or
                            float(entry['best_padj']) < float(summary['best_padj'])):
                        summary['best_padj'] = entry['best_padj']

    rows = list(genes.values())
    columns = ['comparison_id', 'gene_symbol', 'gene_id', 'tool',
               'group1_name', 'group2_name', 'effect_size_direction',
               'n_significant_features', 'n_ok_features',
               'best_effect_size', 'best_padj']
    write_rows(os.path.join(out_dir, '%s.cross_tool.master.tsv' % comparison_id),
               columns, rows)

    gene_rows = []
    for summary in per_gene.values():
        gene_rows.append({
            'comparison_id': comparison_id,
            'gene_symbol': summary['gene_symbol'],
            'gene_id': summary['gene_id'],
            'group1_name': group1_name,
            'group2_name': group2_name,
            'effect_size_direction': 'group2_minus_group1',
            'n_tools_significant': len(summary['tools']),
            'tools': ','.join(sorted(summary['tools'])),
            'best_effect_size': summary['best_effect_size'],
            'best_padj': summary['best_padj'],
        })
    gene_rows.sort(key=lambda r: (-r['n_tools_significant'], r['gene_symbol']))
    gene_columns = ['comparison_id', 'gene_symbol', 'gene_id',
                    'group1_name', 'group2_name', 'effect_size_direction',
                    'n_tools_significant', 'tools', 'best_effect_size', 'best_padj']
    write_rows(os.path.join(out_dir, '%s.cross_tool.gene_summary.tsv' % comparison_id),
               gene_columns, gene_rows)

    print('[export_tool_masters] cross_tool: %d gene x tool rows, %d genes'
          % (len(rows), len(gene_rows)))
    print('[export_tool_masters] cross_tool per-tool significant/ok genes: %s'
          % json.dumps({t: [len(s['significant']), len(s['ok'])]
                        for t, s in tool_stats.items()}))
    return 0


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
        return export_cross_tool(args.comparison_id, args.master_files or [],
                                 args.out_dir, params)
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
