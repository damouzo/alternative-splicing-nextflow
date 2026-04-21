#!/usr/bin/env python3
"""
Prepare rMATS input files (b1.txt, b2.txt or single-sample BAM list).

rMATS requires:
- For PREP: single BAM file per run (one-line text file)
- For POST: comma-separated BAM lists for group1 (b1.txt) and group2 (b2.txt)

All paths must be absolute.
"""

import sys
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare rMATS input files"
    )
    parser.add_argument(
        "--bam_files",
        type=str,
        required=True,
        help="Comma-separated list of BAM file paths"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file path"
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=['single', 'list'],
        default='list',
        help="Output mode: 'single' for one BAM per line, 'list' for comma-separated"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Parse BAM file list
    bam_paths = [p.strip() for p in args.bam_files.split(',') if p.strip()]
    
    if not bam_paths:
        sys.exit("ERROR: No BAM files provided")
    
    # Verify all BAM files exist and convert to absolute paths
    absolute_paths = []
    for bam_path in bam_paths:
        p = Path(bam_path).resolve()
        if not p.exists():
            sys.exit(f"ERROR: BAM file not found: {bam_path}")
        if not p.is_file():
            sys.exit(f"ERROR: Not a file: {bam_path}")
        absolute_paths.append(str(p))
    
    # Write output
    with open(args.output, 'w') as f:
        if args.mode == 'single':
            # One BAM per line (for PREP with multiple samples)
            for path in absolute_paths:
                f.write(f"{path}\n")
        else:
            # Comma-separated list (for POST b1.txt/b2.txt)
            f.write(','.join(absolute_paths))
    
    print(f"Created rMATS input file: {args.output}")
    print(f"  Mode: {args.mode}")
    print(f"  BAM files: {len(absolute_paths)}")


if __name__ == "__main__":
    sys.exit(main())
