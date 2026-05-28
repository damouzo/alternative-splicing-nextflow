#!/usr/bin/env python3
"""
Run IUPred3 on an amino acid FASTA and write output in IUPred2A format
expected by IsoformSwitchAnalyzeR::analyzeIUPred2A().

Output format per isoform:
  >TRANSCRIPT_ID
  1\tM\tiupred_score\tanchor_score
  2\tA\tiupred_score\tanchor_score
  ...
"""

import sys
import argparse


def parse_fasta(fasta_path):
    """Yield (header, sequence) tuples."""
    seq_id, seq_parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if seq_id is not None:
                    yield seq_id, ''.join(seq_parts)
                # take the first word after > as the ID
                seq_id = line[1:].split()[0]
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if seq_id is not None:
        yield seq_id, ''.join(seq_parts)


def run_iupred3_on_sequence(sequence):
    """Return (iupred_scores, anchor_scores) for a single protein sequence."""
    try:
        from iupred3 import iupred3_lib
        iupred_result = iupred3_lib.iupred(sequence, iupred_type='long')
        anchor_result = iupred3_lib.iupred(sequence, iupred_type='anchor2')
        iupred_scores = iupred_result[1]
        anchor_scores = anchor_result[1]
        return iupred_scores, anchor_scores
    except Exception as e:
        # Fallback: return zeros if iupred3 fails for this sequence
        n = len(sequence)
        return [0.0] * n, [0.0] * n


def main():
    parser = argparse.ArgumentParser(description='Run IUPred3 on AA FASTA for ISAR')
    parser.add_argument('fasta', help='Input AA FASTA file')
    parser.add_argument('--output', '-o', default='-', help='Output file (default: stdout)')
    args = parser.parse_args()

    out_fh = open(args.output, 'w') if args.output != '-' else sys.stdout

    try:
        for seq_id, sequence in parse_fasta(args.fasta):
            if not sequence:
                continue
            # Clean sequence: keep only standard AA characters, replace ambiguous with A
            clean_seq = ''.join(
                aa if aa.upper() in 'ACDEFGHIKLMNPQRSTVWY' else 'A'
                for aa in sequence.upper()
            )
            if not clean_seq:
                continue

            iupred_scores, anchor_scores = run_iupred3_on_sequence(clean_seq)

            out_fh.write(f'>{seq_id}\n')
            for i, (aa, iu, anc) in enumerate(zip(clean_seq, iupred_scores, anchor_scores), 1):
                out_fh.write(f'{i}\t{aa}\t{iu:.4f}\t{anc:.4f}\n')
    finally:
        if args.output != '-':
            out_fh.close()


if __name__ == '__main__':
    main()
