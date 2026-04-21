#!/usr/bin/env python3
"""
Generate MAJIQ configuration file (majiq.conf).

MAJIQ requires an INI-style config with:
- [info] section: readlen, bamdirs, genome, strandness
- [experiments] section: sample_id = bam_filename mappings

All BAM files must be accessible from directories listed in bamdirs.
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate MAJIQ configuration file"
    )
    parser.add_argument(
        "--samples",
        type=str,
        required=True,
        help="Sample information: sample_id:bam_path pairs, comma-separated"
    )
    parser.add_argument(
        "--read_length",
        type=int,
        required=True,
        help="Read length (after trimming)"
    )
    parser.add_argument(
        "--strandedness",
        type=str,
        required=True,
        choices=['unstranded', 'forward', 'reverse'],
        help="Library strandedness"
    )
    parser.add_argument(
        "--genome_build",
        type=str,
        default='hg38',
        help="Genome build identifier (e.g., hg38, mm10)"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output config file path"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Parse sample information
    # Format: sample1:path1,sample2:path2,...
    samples = []
    bam_dirs = set()
    
    for entry in args.samples.split(','):
        entry = entry.strip()
        if not entry:
            continue
        
        try:
            sample_id, bam_path = entry.split(':', 1)
        except ValueError:
            sys.exit(f"ERROR: Invalid sample entry format: {entry}. Expected 'sample_id:bam_path'")
        
        bam_path = Path(bam_path).resolve()
        
        if not bam_path.exists():
            sys.exit(f"ERROR: BAM file not found: {bam_path}")
        
        # Get directory containing BAM
        bam_dir = bam_path.parent
        bam_dirs.add(str(bam_dir))
        
        samples.append({
            'id': sample_id,
            'bam': bam_path.name,  # Just filename for config
            'dir': str(bam_dir)
        })
    
    if not samples:
        sys.exit("ERROR: No samples provided")
    
    # Map strandedness to MAJIQ format
    strandness_map = {
        'unstranded': 'none',
        'forward': 'forward',
        'reverse': 'reverse'
    }
    majiq_strandness = strandness_map[args.strandedness]
    
    # Generate config file
    with open(args.output, 'w') as f:
        # [info] section
        f.write("[info]\n")
        f.write(f"readlen={args.read_length}\n")
        f.write(f"bamdirs={','.join(sorted(bam_dirs))}\n")
        f.write(f"genome={args.genome_build}\n")
        f.write(f"strandness={majiq_strandness}\n")
        f.write("\n")
        
        # [experiments] section
        f.write("[experiments]\n")
        for sample in samples:
            f.write(f"{sample['id']}={sample['bam']}\n")
    
    print(f"Created MAJIQ config: {args.output}")
    print(f"  Samples: {len(samples)}")
    print(f"  BAM directories: {len(bam_dirs)}")
    print(f"  Strandness: {majiq_strandness}")


if __name__ == "__main__":
    sys.exit(main())
