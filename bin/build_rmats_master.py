#!/usr/bin/env python3
"""
Build the rMATS CORE deliverables for one comparison:

  <comparison_id>.rmats.master.tsv       all events, all 5 event types
  <comparison_id>.rmats.significant.tsv  events passing the significance rule
  <comparison_id>.rmats.summary.tsv      per event type counts

Significance rule (same as the HTML report):
    is_significant = (padj <= fdr_cutoff) and (|IncLevelDifference| >= dpsi_cutoff)

Design notes:
  - rMATS 4.3 duplicates the ID column in *.MATS.JC.txt; the first occurrence
    holds the event ID (matches fromGTF.novelSpliceSite.<ET>.txt), the second
    is a duplicate. We keep the first and rename the second.
  - FDR == 0 rows are kept as-is and flagged with fdr_floor_flag (numeric floor
    of --cstat); do not round them silently.
  - Native per-sample PSI vectors are summarised to group means (inc_level_1/2),
    matching what the report's PSI section consumes.
  - Contract 2.1.0 (additive): every row carries event_locus (chr:start-end,
    1-based, IGV-ready) and event_coords (native 0-based coordinates named per
    event type). The six SE-shaped columns (exon_start_0base..downstream_ee)
    are legacy: only SE rows populate them.

Rows are streamed event-type by event-type to keep peak memory flat on large
runs. Stdlib only.
"""

import argparse
import csv
import os
import statistics
import sys

EVENT_TYPES = ('SE', 'A5SS', 'A3SS', 'MXE', 'RI')

# Paste-ready location columns, independent of the event type (contract 2.1.0).
LOCATION_COLUMNS = [
    'event_locus',   # chr:start-end (1-based inclusive) spanning every event coordinate
    'event_coords',  # native rMATS coordinates, 0-based, named per type (EVENT_COORD_FIELDS)
]

# Native coordinate pairs per event type, in rMATS column order.
EVENT_COORD_FIELDS = {
    'SE':   (('exon', 'exonStart_0base', 'exonEnd'),
             ('upstream', 'upstreamES', 'upstreamEE'),
             ('downstream', 'downstreamES', 'downstreamEE')),
    'A5SS': (('long', 'longExonStart_0base', 'longExonEnd'),
             ('short', 'shortES', 'shortEE'),
             ('flanking', 'flankingES', 'flankingEE')),
    'A3SS': (('long', 'longExonStart_0base', 'longExonEnd'),
             ('short', 'shortES', 'shortEE'),
             ('flanking', 'flankingES', 'flankingEE')),
    'MXE':  (('1stExon', '1stExonStart_0base', '1stExonEnd'),
             ('2ndExon', '2ndExonStart_0base', '2ndExonEnd'),
             ('upstream', 'upstreamES', 'upstreamEE'),
             ('downstream', 'downstreamES', 'downstreamEE')),
    'RI':   (('riExon', 'riExonStart_0base', 'riExonEnd'),
             ('upstream', 'upstreamES', 'upstreamEE'),
             ('downstream', 'downstreamES', 'downstreamEE')),
}

STANDARD_COLUMNS = [
    'comparison_id', 'tool', 'event_type', 'event_id', 'feature_type',
    'feature_id', 'gene_id', 'gene_symbol',
    'group1_name', 'group2_name', 'effect_size_direction',
    'chr', 'strand', 'inc_level_1', 'inc_level_2', 'inc_level_difference',
    'effect_size', 'effect_size_type', 'pvalue', 'padj', 'padj_method',
    'fdr_floor_flag', 'is_novel_splice_site',
    'is_significant', 'significance_rule', 'source_file',
]

NATIVE_COLUMNS = [
    # rMATS-native IncLevelDifference (group1 - group2). Kept verbatim so the
    # sign flip into group2 - group1 is auditable; all other columns use the
    # unified convention.
    'inc_level_difference_native',
    # Junction-read support (IJC + SJC) per group and the report_min_reads gate.
    'reads_group1', 'reads_group2', 'reads_ok',
    # annotated | novel_junction | novel_splice_site (from the fromGTF.* files).
    'event_class',
    # canonical assembly chromosome (IGV-able); false for scaffolds/MT.
    'is_igv_supported',
    'ijc_sample_1_total', 'sjc_sample_1_total',
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


def as_int(value):
    """Integer coordinate or None; tolerates float-like strings from rMATS."""
    if value is None:
        return None
    text = str(value).strip().strip('"')
    if not text or text == 'NA':
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


CANONICAL_CHROMS = set(str(i) for i in range(1, 23)) | {'X', 'Y'}


def normalize_chrom(chrom):
    """Map rMATS/LeafCutter 'chr1' to MAJIQ-style '1'; keep MT as 'MT'."""
    value = (chrom or '').strip().strip('"')
    if value[:3].lower() == 'chr':
        value = value[3:]
    upper = value.upper()
    if upper in ('M', 'MT'):
        return 'MT'
    if upper == 'X':
        return 'X'
    if upper == 'Y':
        return 'Y'
    return value


def is_igv_supported(chrom):
    """True for the canonical assembly chromosomes; scaffolds/MT are excluded."""
    return chrom in CANONICAL_CHROMS


def build_event_locus(event_type, row, header_names):
    """(event_locus, event_coords) for one event; ('', '') when unparseable.

    event_locus spans every named coordinate of the event, 1-based inclusive
    (min start + 1 .. max end), so it can be pasted straight into IGV. The
    chromosome name is harmonized to the MAJIQ style (no 'chr' prefix).
    event_coords keeps the raw 0-based values with per-type names.
    """
    pairs = []
    for label, start_col, end_col in EVENT_COORD_FIELDS.get(event_type, ()):
        if start_col not in header_names or end_col not in header_names:
            continue
        start = as_int(row.get(start_col))
        end = as_int(row.get(end_col))
        if start is None or end is None:
            continue
        pairs.append((label, start, end))
    if not pairs:
        return '', ''
    locus = '%d-%d' % (min(pair[1] for pair in pairs) + 1, max(pair[2] for pair in pairs))
    chrom = normalize_chrom(row.get('chr'))
    if chrom:
        locus = '%s:%s' % (chrom, locus)
    return locus, ';'.join('%s=%d-%d' % pair for pair in pairs)


def load_id_set(path):
    """rMATS fromGTF.*.txt ID column as a set (empty set when absent)."""
    ids = set()
    if not os.path.isfile(path):
        return ids
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        if 'ID' not in header:
            return ids
        idx = header.index('ID')
        for line in handle:
            if line.strip():
                parts = line.rstrip('\n').split('\t')
                if len(parts) > idx:
                    ids.add(parts[idx])
    return ids


def iter_jc_rows(comparison_id, rmats_dir, event_type, fdr_cutoff, dpsi_cutoff,
                 rule, group1_name, group2_name, min_reads):
    """Yield (master_row_dict, is_significant) for one MATS.JC file."""
    jc_path = os.path.join(rmats_dir, '%s.MATS.JC.txt' % event_type)
    if not os.path.isfile(jc_path):
        return
    novel_ids = load_id_set(
        os.path.join(rmats_dir, 'fromGTF.novelSpliceSite.%s.txt' % event_type))
    novel_junction_ids = load_id_set(
        os.path.join(rmats_dir, 'fromGTF.novelJunction.%s.txt' % event_type))

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
            # Unified sign: effect_size = group2 - group1. rMATS ships
            # IncLevelDifference as group1 - group2, so negate it and keep the
            # native value in inc_level_difference_native.
            effect_dpsi = -dpsi if dpsi is not None else None
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
                'gene_symbol': '' if (row.get('geneSymbol') or '').strip('"') in ('', 'NA')
                               else (row.get('geneSymbol') or '').strip('"'),
                'group1_name': group1_name,
                'group2_name': group2_name,
                'effect_size_direction': 'group2_minus_group1',
                'chr': normalize_chrom(row.get('chr')),
                'strand': row.get('strand', ''),
                'inc_level_1': parse_psi_vector(row.get('IncLevel1')),
                'inc_level_2': parse_psi_vector(row.get('IncLevel2')),
                'inc_level_difference': num(effect_dpsi),
                'inc_level_difference_native': num(dpsi),
                'effect_size': num(effect_dpsi),
                'effect_size_type': 'delta_psi',
                'pvalue': num(row.get('PValue')),
                'padj': num(row.get('FDR')),
                'padj_method': 'rmats_cstat',
                'fdr_floor_flag': 'true' if fdr_float == 0.0 else 'false',
                'is_novel_splice_site': 'true' if row.get('ID') in novel_ids else 'false',
                'is_significant': 'true' if is_significant else 'false',
                'significance_rule': rule,
                'source_file': '%s.MATS.JC.txt' % event_type,
            }
            for native in LOCATION_COLUMNS + NATIVE_COLUMNS:
                out[native] = ''
            # Native value cleared by the loop above; restore it after.
            out['inc_level_difference_native'] = num(dpsi)
            # Junction-read coverage per group (IJC + SJC) and the report gate.
            ijc1 = parse_count_vector(row.get('IJC_SAMPLE_1')) or 0
            sjc1 = parse_count_vector(row.get('SJC_SAMPLE_1')) or 0
            ijc2 = parse_count_vector(row.get('IJC_SAMPLE_2')) or 0
            sjc2 = parse_count_vector(row.get('SJC_SAMPLE_2')) or 0
            reads1 = ijc1 + sjc1
            reads2 = ijc2 + sjc2
            out['reads_group1'] = reads1
            out['reads_group2'] = reads2
            out['reads_ok'] = 'true' if (reads1 >= min_reads and reads2 >= min_reads) else 'false'
            event_id = row.get('ID')
            if event_id in novel_ids:
                out['event_class'] = 'novel_splice_site'
            elif event_id in novel_junction_ids:
                out['event_class'] = 'novel_junction'
            else:
                out['event_class'] = 'annotated'
            out['is_igv_supported'] = 'true' if is_igv_supported(
                normalize_chrom(row.get('chr'))) else 'false'
            out['event_locus'], out['event_coords'] = build_event_locus(
                event_type, row, full)
            if 'IJC_SAMPLE_1' in full:
                count = parse_count_vector(row.get('IJC_SAMPLE_1'))
                out['ijc_sample_1_total'] = '' if count is None else count
            if 'SJC_SAMPLE_1' in full:
                count = parse_count_vector(row.get('SJC_SAMPLE_1'))
                out['sjc_sample_1_total'] = '' if count is None else count
            if 'IJC_SAMPLE_2' in full:
                count = parse_count_vector(row.get('IJC_SAMPLE_2'))
                out['ijc_sample_2_total'] = '' if count is None else count
            if 'SJC_SAMPLE_2' in full:
                count = parse_count_vector(row.get('SJC_SAMPLE_2'))
                out['sjc_sample_2_total'] = '' if count is None else count
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


def sort_number(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-id', required=True)
    parser.add_argument('--rmats-dir', required=True)
    parser.add_argument('--fdr-cutoff', type=float, required=True)
    parser.add_argument('--dpsi-cutoff', type=float, required=True)
    parser.add_argument('--group1-name', default='')
    parser.add_argument('--group2-name', default='')
    parser.add_argument('--report-min-reads', type=int, default=0,
                        help='reads_ok requires >= this many junction reads (IJC+SJC) per group')
    parser.add_argument('--out-dir', default='.')
    args = parser.parse_args()

    rule = 'padj <= %s & |inc_level_difference| >= %s' % (args.fdr_cutoff, args.dpsi_cutoff)
    prefix = '%s.rmats' % args.comparison_id
    master_path = os.path.join(args.out_dir, prefix + '.master.tsv')
    sig_path = os.path.join(args.out_dir, prefix + '.significant.tsv')
    summary_path = os.path.join(args.out_dir, prefix + '.summary.tsv')

    columns = STANDARD_COLUMNS + LOCATION_COLUMNS + NATIVE_COLUMNS
    summary_cols = ['comparison_id', 'tool', 'event_type', 'n_total', 'n_significant',
                    'n_novel_splice_site', 'median_abs_dpsi', 'p90_abs_dpsi']

    n_total = 0
    n_significant = 0
    summary_rows = []
    sig_rows = []

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
                    args.fdr_cutoff, args.dpsi_cutoff, rule,
                    args.group1_name, args.group2_name, args.report_min_reads):
                master_writer.writerow(row)
                type_total += 1
                if row['is_novel_splice_site'] == 'true':
                    type_novel += 1
                if is_significant:
                    sig_rows.append(row)
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

        # significant.tsv is ordered for reading: most significant first, ties by
        # strongest effect size.
        sig_rows.sort(key=lambda row: (sort_number(row.get('padj'), float('inf')),
                                       -abs(sort_number(row.get('effect_size'), 0.0))))
        for row in sig_rows:
            sig_writer.writerow(row)

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