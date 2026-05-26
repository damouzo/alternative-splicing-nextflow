#!/usr/bin/env python3
"""
Validate samplesheet and check critical requirements for alternative splicing analysis.

This script performs comprehensive validation:
1. Column names and format
2. File existence and readability
3. Condition and replicate structure
4. Chromosome naming consistency (BAM vs GTF)
5. Read length verification via pysam

Critical checks prevent silent pipeline failures.
"""

import sys
import argparse
import csv
from pathlib import Path
from collections import defaultdict, Counter
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate samplesheet for alternative splicing pipeline"
    )
    parser.add_argument(
        "--samplesheet",
        type=str,
        required=True,
        help="Path to samplesheet CSV"
    )
    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="Path to GTF annotation file"
    )
    parser.add_argument(
        "--read_length",
        type=int,
        required=True,
        help="Expected read length (after trimming)"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to validated output CSV"
    )
    return parser.parse_args()


def check_file_exists(filepath, description):
    """Check if file exists and is readable."""
    p = Path(filepath)
    if not p.exists():
        sys.exit(f"ERROR: {description} not found: {filepath}")
    if not p.is_file():
        sys.exit(f"ERROR: {description} is not a file: {filepath}")
    return p


def get_gtf_chr_prefix(gtf_path):
    """Determine if GTF uses 'chr' prefix by checking first sequence name."""
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                chr_name = fields[0]
                return chr_name.startswith('chr')
    sys.exit(f"ERROR: Could not determine chromosome naming from GTF: {gtf_path}")


def get_bam_chr_prefix(bam_path):
    """Check if BAM uses 'chr' prefix in @SQ header lines."""
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            if len(bam.references) == 0:
                sys.exit(f"ERROR: BAM has no reference sequences: {bam_path}")
            first_chr = bam.references[0]
            return first_chr.startswith('chr')
    except Exception as e:
        sys.exit(f"ERROR: Could not read BAM header from {bam_path}: {e}")


def get_bam_read_length(bam_path, n_reads=10000):
    """
    Extract modal read length from BAM by sampling aligned reads.
    Returns the most common read length observed.
    """
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            lengths = []
            for i, read in enumerate(bam.fetch(until_eof=True)):
                if i >= n_reads:
                    break
                if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    lengths.append(read.query_length)
            
            if not lengths:
                sys.exit(f"ERROR: No aligned reads found in BAM: {bam_path}")
            
            # Get modal (most common) read length
            length_counts = Counter(lengths)
            modal_length = length_counts.most_common(1)[0][0]
            return modal_length
    except Exception as e:
        sys.exit(f"ERROR: Could not extract read lengths from {bam_path}: {e}")


def validate_samplesheet(samplesheet_path, gtf_path, expected_read_length):
    """
    Validate samplesheet format and content.
    Returns list of validated sample dicts.
    """
    required_cols = ['sample', 'condition', 'replicate', 'bam', 'bai', 'salmon_dir']
    
    # Read samplesheet
    samples = []
    try:
        with open(samplesheet_path, 'r') as f:
            reader = csv.DictReader(f)
            
            # Check column names
            if not reader.fieldnames:
                sys.exit("ERROR: Samplesheet is empty or malformed")
            
            missing_cols = set(required_cols) - set(reader.fieldnames)
            if missing_cols:
                sys.exit(f"ERROR: Missing required columns: {', '.join(missing_cols)}")
            
            for i, row in enumerate(reader, start=2):  # start=2 because row 1 is header
                # Check for empty fields
                for col in required_cols:
                    if not row[col] or row[col].strip() == '':
                        sys.exit(f"ERROR: Empty value in column '{col}' at row {i}")
                
                samples.append(row)
    
    except FileNotFoundError:
        sys.exit(f"ERROR: Samplesheet not found: {samplesheet_path}")
    except Exception as e:
        sys.exit(f"ERROR: Failed to parse samplesheet: {e}")
    
    if not samples:
        sys.exit("ERROR: Samplesheet contains no data rows")
    
    # Check for duplicate sample IDs
    sample_ids = [s['sample'] for s in samples]
    duplicates = [sid for sid, count in Counter(sample_ids).items() if count > 1]
    if duplicates:
        sys.exit(f"ERROR: Duplicate sample IDs found: {', '.join(duplicates)}")
    
    # Check file existence
    print("Checking file paths...")
    for sample in samples:
        check_file_exists(sample['bam'], f"BAM for sample {sample['sample']}")
        check_file_exists(sample['bai'], f"BAI for sample {sample['sample']}")
        
        salmon_dir = Path(sample['salmon_dir'])
        if not salmon_dir.exists() or not salmon_dir.is_dir():
            sys.exit(f"ERROR: Salmon directory not found for sample {sample['sample']}: {sample['salmon_dir']}")
        
        quant_sf = salmon_dir / 'quant.sf'
        if not quant_sf.exists():
            sys.exit(f"ERROR: quant.sf not found in Salmon directory for sample {sample['sample']}: {quant_sf}")
    
    # Validate conditions
    conditions = set(s['condition'] for s in samples)
    if len(conditions) < 2:
        sys.exit(f"ERROR: Expected at least 2 conditions, found {len(conditions)}: {', '.join(conditions)}")
    print(f"✓ Found {len(conditions)} condition(s): {', '.join(sorted(conditions))}")
    
    # Count replicates per condition
    condition_counts = defaultdict(int)
    for sample in samples:
        condition_counts[sample['condition']] += 1
    
    for condition, count in condition_counts.items():
        if count < 2:
            sys.exit(f"ERROR: Condition '{condition}' has only {count} replicate(s). Minimum 2 required.")
        elif count < 3:
            print(f"WARNING: Condition '{condition}' has only {count} replicates. Minimum 3 recommended for statistical power.", file=sys.stderr)
    
    # Chromosome naming consistency check
    print("Checking chromosome naming consistency...")
    gtf_has_chr = get_gtf_chr_prefix(gtf_path)
    
    # Check first BAM from each condition
    checked_conditions = set()
    for sample in samples:
        if sample['condition'] not in checked_conditions:
            bam_has_chr = get_bam_chr_prefix(sample['bam'])
            
            if gtf_has_chr != bam_has_chr:
                gtf_style = "chr-prefixed (chr1, chr2, ...)" if gtf_has_chr else "bare (1, 2, ...)"
                bam_style = "chr-prefixed (chr1, chr2, ...)" if bam_has_chr else "bare (1, 2, ...)"
                sys.exit(
                    f"ERROR: Chromosome naming mismatch detected!\n"
                    f"  GTF uses {gtf_style}\n"
                    f"  BAM uses {bam_style} (sample: {sample['sample']})\n"
                    f"This will cause rMATS and MAJIQ to detect zero events silently.\n"
                    f"Solution: Use the same reference files (GTF + genome FASTA) in both nf-core/rnaseq and this pipeline."
                )
            
            checked_conditions.add(sample['condition'])
    
    print(f"✓ Chromosome naming is consistent ({'chr-prefixed' if gtf_has_chr else 'bare numbers'})")
    
    # Read length validation
    print("Validating read lengths from BAMs...")
    read_length_tolerance = 10  # bp
    
    for sample in samples:
        observed_length = get_bam_read_length(sample['bam'])
        diff = abs(observed_length - expected_read_length)
        
        if diff > read_length_tolerance:
            sys.exit(
                f"ERROR: Read length mismatch in sample '{sample['sample']}':\n"
                f"  Expected: {expected_read_length} bp (from --read_length parameter)\n"
                f"  Observed: {observed_length} bp (modal length in BAM)\n"
                f"  Difference: {diff} bp (tolerance: {read_length_tolerance} bp)\n"
                f"This will cause incorrect PSI normalization in rMATS and MAJIQ.\n"
                f"Solution: Check actual post-trimming read lengths in nf-core/rnaseq logs and update --read_length parameter."
            )
        
        print(f"  {sample['sample']}: {observed_length} bp (expected {expected_read_length} bp) ✓")
    
    print(f"✓ All samples have consistent read lengths (within {read_length_tolerance} bp tolerance)")
    
    return samples


def write_validated_samplesheet(samples, output_path):
    """Write validated samplesheet to output."""
    fieldnames = ['sample', 'condition', 'replicate', 'bam', 'bai', 'salmon_dir']
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(samples)
    
    print(f"✓ Validated samplesheet written to: {output_path}")


def main():
    args = parse_args()
    
    print("=" * 60)
    print("Alternative Splicing Pipeline - Samplesheet Validation")
    print("=" * 60)
    
    # Validate inputs
    validated_samples = validate_samplesheet(
        args.samplesheet,
        args.gtf,
        args.read_length
    )
    
    # Write output
    write_validated_samplesheet(validated_samples, args.output)
    
    print("=" * 60)
    print(f"✓ Validation complete: {len(validated_samples)} samples passed all checks")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
