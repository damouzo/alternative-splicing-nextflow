#!/usr/bin/env python3
"""
Build the rMATS CORE deliverables for one comparison:

  <comparison_id>.rmats.master.tsv       all events, all 5 event types
  <comparison_id>.rmats.significant.tsv  events passing the significance rule
  <comparison_id>.rmats.summary.tsv      per event type counts

Significance rule (same as the HTML report):
    is_significant = (FDR <= fdr_cutoff) and (|IncLevelDifference| >= dpsi_cutoff)

Design notes:
  - rMATS 4.3 duplicates the ID column in *.MATS.JC.txt; the first occurrence
    holds the event ID (matches fromGTF.novelSpliceSite.<ET>.txt), the second
    is a duplicate. We keep the first and rename the second.
  - FDR == 0 rows are kept as-is and flagged with fdr_floor_flag (numeric floor
    of --cstat); do not round them silently.
  - Native per-sample PSI vectors are summarised to group means (inc_level_1/2),
    matching what the report's PSI section consumes.

Rows are streamed event-type by event-type to keep peak memory flat on large
runs. Stdlib only.
"""

import argparse
import csv
import os
import statistics
import sys

EVENT_TYPES = ('SE', 'A5SS', 'A3SS', 'MXE', 'RI')

STANDARD_COLUMNS = [
    'comparison_id', 'tool', 'event_type', 'event_id', 'feature_type',
    'feature_id', 'gene_id', 'gene_symbol',
    'chr', 'strand', 'inc_level_1', 'inc_level_2', 'inc_level_difference',
    'effect_size', 'pvalue', 'fdr', 'fdr_floor_flag', 'is_novel_splice_site',
    'is_significant', 'significance_rule', 'source_file',
]

NATIVE_COLUMNS = [
    'exon_start_0base', 'exon_end', 'upstream_es', 'upstream_ee',
    'downstream_es', 'downstream_ee', 'ijc_sample_1_total', 'sjc_sample_1_total',
    'ijc_sample_2_total', 'sjc_sample_2_total', 'inc_form_len', 'skip_form_len',
]


def parse_psi_vector(value):
    """Mean of the comma-separated per-sample PSI vector; None for NA/empty."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == 'NA':
        return None
    parts = text.split(',')
    values = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            continue
    if not values:
        return None
    return round(statistics.mean(values), 6)


def parse_count_vector(value):
    """Sum of IJC/SJC comma-separated counts; None when not countable."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == 'NA':
        return None
    total = 0
    for part in text.split(','):
        try:
            total += int(part)
        except ValueError:
            continue
    return total


def norm_gene_id(value):
    text = (value or '').strip().strip('"')
    return text if text and text != 'NA' else ''


def num(value):
    if value is None:
        return ''
    text = str(value).strip().strip('"')
    if not text or text == 'NA':
        return ''
    try:
        number = float(text)
        if number.is_integer() and abs(number) < 1e15:
            return str(int(number))
        return repr(number)
    except ValueError:
        return text


def iter_jc_rows(comparison_id, rmats_dir, event_type, fdr_cutoff, dpsi_cutoff, rule):
    """Yield (master_row_dict, is_significant) for one MATS.JC file."""
    jc_path = os.path.join(rmats_dir, '%s.MATS.JC.txt' % event_type)
    if not os.path.isfile(jc_path):
        return
    novel_path = os.path.join(rmats_dir, 'fromGTF.novelSpliceSite.%s.txt' % event_type)
    novel_ids = set()
    if os.path.isfile(novel_path):
        with open(novel_path, 'r', encoding='utf-8', errors='replace') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            if 'ID' in header:
                idx = header.index('ID')
                for line in handle:
                    if line.strip():
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) > idx:
                            novel_ids.add(parts[idx])

    with open(jc_path, 'r', encoding='utf-8', errors='replace') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        names = []
        seen = {}
        for col in header:
            if col in seen:
                seen[col] += 1
                names.append(col + '_dup' + str(seen[col]))
            else:
                seen[col] = 0
                names.append(col)
        full = set(names)
        for line in handle:
            if not line.strip():
                continue
            values = line.rstrip('\n').split('\t')
            if len(values) > len(names):
                values = values[:len(names)]
            elif len(values) < len(names):
                values += [''] * (len(names) - len(values))
            row = dict(zip(names, values))

            fdr_float = None
            try:
                fdr_float = float(row.get('FDR'))
            except (TypeError, ValueError):
                pass
            dpsi = None
            try:
                dpsi = float(row.get('IncLevelDifference'))
            except (TypeError, ValueError):
                pass
            is_significant = bool(
                fdr_float is not None and dpsi is not None and
                fdr_float <= fdr_cutoff and abs(dpsi) >= dpsi_cutoff)

            out = {
                'comparison_id': comparison_id,
                'tool': 'rmats',
                'event_type': event_type,
                'event_id': row.get('ID', ''),
                'feature_type': event_type,
                'feature_id': row.get('ID', ''),
                'gene_id': norm_gene_id(row.get('GeneID')),
                'gene_symbol': (row.get('geneSymbol') or '').strip('"'),
                'chr': row.get('chr', ''),
                'strand': row.get('strand', ''),
                'inc_level_1': parse_psi_vector(row.get('IncLevel1')),
                'inc_level_2': parse_psi_vector(row.get('IncLevel2')),
                'inc_level_difference': num(row.get('IncLevelDifference')),
                'effect_size': num(row.get('IncLevelDifference')),
                'pvalue': num(row.get('PValue')),
                'fdr': num(row.get('FDR')),
                'fdr_floor_flag': 'true' if fdr_float == 0.0 else 'false',
                'is_novel_splice_site': 'true' if row.get('ID') in novel_ids else 'false',
                'is_significant': 'true' if is_significant else 'false',
                'significance_rule': rule,
                'source_file': '%s.MATS.JC.txt' % event_type,
            }
            for native in NATIVE_COLUMNS:
                out[native] = ''
            if 'exonStart_0base' in full:
                out['exon_start_0base'] = num(row.get('exonStart_0base'))
            if 'exonEnd' in full:
                out['exon_end'] = num(row.get('exonEnd'))
            if 'upstreamES' in full:
                out['upstream_es'] = num(row.get('upstreamES'))
            if 'upstreamEE' in full:
                out['upstream_ee'] = num(row.get('upstreamEE'))
            if 'downstreamES' in full:
                out['downstream_es'] = num(row.get('downstreamES'))
            if 'downstreamEE' in full:
                out['downstream_ee'] = num(row.get('downstreamEE'))
            if 'IJC_SAMPLE_1' in full:
                out['ijc_sample_1_total'] = parse_count_vector(row.get('IJC_SAMPLE_1')) or ''
            if 'SJC_SAMPLE_1' in full:
                out['sjc_sample_1_total'] = parse_count_vector(row.get('SJC_SAMPLE_1')) or ''
            if 'IJC_SAMPLE_2' in full:
                out['ijc_sample_2_total'] = parse_count_vector(row.get('IJC_SAMPLE_2')) or ''
            if 'SJC_SAMPLE_2' in full:
                out['sjc_sample_2_total'] = parse_count_vector(row.get('SJC_SAMPLE_2')) or ''
            if 'IncFormLen' in full:
                out['inc_form_len'] = num(row.get('IncFormLen'))
            if 'SkipFormLen' in full:
                out['skip_form_len'] = num(row.get('SkipFormLen'))
            yield out, is_significant


def percentile_sorted(values, q):
    if not values:
        return ''
    idx = int(round(q * (len(values) - 1)))
    return values[idx]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-id', required=True)
    parser.add_argument('--rmats-dir', required=True)
    parser.add_argument('--fdr-cutoff', type=float, required=True)
    parser.add_argument('--dpsi-cutoff', type=float, required=True)
    parser.add_argument('--out-dir', default='.')
    args = parser.parse_args()

    rule = 'fdr <= %s & |inc_level_difference| >= %s' % (args.fdr_cutoff, args.dpsi_cutoff)
    prefix = '%s.rmats' % args.comparison_id
    master_path = os.path.join(args.out_dir, prefix + '.master.tsv')
    sig_path = os.path.join(args.out_dir, prefix + '.significant.tsv')
    summary_path = os.path.join(args.out_dir, prefix + '.summary.tsv')

    columns = STANDARD_COLUMNS + NATIVE_COLUMNS
    summary_cols = ['comparison_id', 'tool', 'event_type', 'n_total', 'n_significant',
                    'n_novel_splice_site', 'median_abs_dpsi', 'p90_abs_dpsi']

    n_total = 0
    n_significant = 0
    summary_rows = []

    with open(master_path, 'w', encoding='utf-8', newline='') as master_handle, \
         open(sig_path, 'w', encoding='utf-8', newline='') as sig_handle:
        master_writer = csv.DictWriter(master_handle, fieldnames=columns,
                                       delimiter='\t', lineterminator='\n',
                                       extrasaction='ignore')
        sig_writer = csv.DictWriter(sig_handle, fieldnames=columns,
                                    delimiter='\t', lineterminator='\n',
                                    extrasaction='ignore')
        master_writer.writeheader()
        sig_writer.writeheader()

        for event_type in EVENT_TYPES:
            type_total = 0
            type_sig = 0
            type_novel = 0
            dpsi_abs = []
            for row, is_significant in iter_jc_rows(
                    args.comparison_id, args.rmats_dir, event_type,
                    args.fdr_cutoff, args.dpsi_cutoff, rule):
                master_writer.writerow(row)
                type_total += 1
                if row['is_novel_splice_site'] == 'true':
                    type_novel += 1
                if is_significant:
                    sig_writer.writerow(row)
                    type_sig += 1
                try:
                    dpsi_abs.append(abs(float(row['inc_level_difference'])))
                except (TypeError, ValueError):
                    pass
            if type_total == 0:
                continue
            dpsi_abs.sort()
            summary_rows.append({
                'comparison_id': args.comparison_id,
                'tool': 'rmats',
                'event_type': event_type,
                'n_total': type_total,
                'n_significant': type_sig,
                'n_novel_splice_site': type_novel,
                'median_abs_dpsi': percentile_sorted(dpsi_abs, 0.5),
                'p90_abs_dpsi': percentile_sorted(dpsi_abs, 0.9),
            })
            n_total += type_total
            n_significant += type_sig

    if n_total == 0:
        print('[build_rmats_master] no events found in %s' % args.rmats_dir,
              file=sys.stderr)
        sys.exit(1)

    with open(summary_path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_cols,
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print('[build_rmats_master] %d events (%d significant) | rule: %s'
          % (n_total, n_significant, rule))


if __name__ == '__main__':
    sys.exit(main())