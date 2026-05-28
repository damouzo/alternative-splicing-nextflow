#!/usr/bin/env python3
"""
Filter rMATS output to a top-N subset for sashimi plot generation.

Selects the top N events per event type ranked by a combined priority score:
  priority = -log10(FDR) * |IncLevelDifference|

Only events passing the FDR and |dPSI| cutoffs are considered. Outputs one
filtered file per event type that passes the selection.
"""

import argparse
import math
import os
import sys


EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]


def priority_score(fdr: float, dpsi: float) -> float:
    """Higher is more significant + larger effect."""
    fdr_clamped = max(fdr, 1e-300)
    return -math.log10(fdr_clamped) * abs(dpsi)


def filter_events(rmats_dir: str, out_dir: str, top_n: int,
                  fdr_cutoff: float, dpsi_cutoff: float) -> None:
    os.makedirs(out_dir, exist_ok=True)
    selected = {}

    for event_type in EVENT_TYPES:
        jc_file = os.path.join(rmats_dir, f"{event_type}.MATS.JC.txt")
        if not os.path.isfile(jc_file):
            continue

        with open(jc_file) as fh:
            header = fh.readline().rstrip("\n").split("\t")

        try:
            fdr_idx  = header.index("FDR")
            dpsi_idx = header.index("IncLevelDifference")
        except ValueError:
            print(f"[WARN] {event_type}: missing FDR or IncLevelDifference column — skipping",
                  file=sys.stderr)
            continue

        candidates = []
        with open(jc_file) as fh:
            fh.readline()  # skip header
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(fdr_idx, dpsi_idx):
                    continue
                try:
                    fdr_val  = float(parts[fdr_idx])
                    dpsi_val = float(parts[dpsi_idx])
                except ValueError:
                    continue

                if fdr_val <= fdr_cutoff and abs(dpsi_val) >= dpsi_cutoff:
                    score = priority_score(fdr_val, dpsi_val)
                    candidates.append((score, line))

        if not candidates:
            continue

        # Sort descending by priority score; take top N
        candidates.sort(key=lambda x: x[0], reverse=True)
        top = candidates[:top_n]

        out_file = os.path.join(out_dir, f"{event_type}.top.txt")
        with open(out_file, "w") as fh:
            # Write header as-is from original file
            with open(jc_file) as src:
                fh.write(src.readline())
            for _, line in top:
                fh.write(line)

        selected[event_type] = len(top)
        print(f"[INFO] {event_type}: selected {len(top)} / {len(candidates)} significant events")

    if not selected:
        print("[WARN] No significant events found for any event type. "
              "Check FDR/dPSI cutoffs.", file=sys.stderr)
        sys.exit(0)

    # Write a summary JSON for the report
    import json
    summary = {"top_n": top_n, "fdr_cutoff": fdr_cutoff, "dpsi_cutoff": dpsi_cutoff,
               "selected": selected}
    with open(os.path.join(out_dir, "filter_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("rmats_dir",   help="Directory with rMATS output files")
    p.add_argument("out_dir",     help="Output directory for filtered event files")
    p.add_argument("--top-n",     type=int,   default=10,   help="Events per event type [10]")
    p.add_argument("--fdr",       type=float, default=0.05, help="FDR cutoff [0.05]")
    p.add_argument("--dpsi",      type=float, default=0.1,  help="|dPSI| cutoff [0.1]")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    filter_events(
        rmats_dir  = args.rmats_dir,
        out_dir    = args.out_dir,
        top_n      = args.top_n,
        fdr_cutoff = args.fdr,
        dpsi_cutoff = args.dpsi,
    )
