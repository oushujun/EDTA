#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict

# Calculate genome coverage for filtered TE families from TE_consistency.py output
# Shujun Ou (shujun.ou.1@gmail.com)
# v0.1: 03/10/2026
# Usage: 
#   awk '$1=="Fam" && $3>10 && $6>0.9 && $7<10' consistency.sum | python3 TE_coverage.py -anno TEanno.out


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate genome coverage for filtered TE families. "
                    "Reads family names from stdin (piped consistency Fam rows).",
        epilog="Example: awk '$1==\"Fam\" && $3>10 && $6>0.9 && $7<10' "
               "consistency.txt | python3 TE_coverage.py -anno TEanno.out")
    parser.add_argument("-anno", required=True,
                        help="TE annotation in RepeatMasker .out format")
    return parser.parse_args()

def main():
    args = parse_args()

    # Step 1: Read family names from stdin (Fam rows from consistency output)
    families = set()
    for line in sys.stdin:
        parts = line.strip().split()
        if not parts:
            continue
        # Accept both raw family names and full Fam rows
        if parts[0] == "Fam" and len(parts) >= 2:
            families.add(parts[1])
        elif not parts[0].startswith(("Ind", "Fam/Ind", "---")):
            families.add(parts[0])

    if not families:
        print("No families found in stdin.", file=sys.stderr)
        sys.exit(1)

    # Step 2: Parse anno file and compute per-family coverage
    fam_bp = defaultdict(int)
    fam_count = defaultdict(int)
    total_anno_bp = 0

    with open(args.anno) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 11:
                continue
            try:
                start = int(parts[5])
                end = int(parts[6])
            except ValueError:
                continue
            bp = end - start + 1
            total_anno_bp += bp
            fam = parts[9]
            if fam in families:
                fam_bp[fam] += bp
                fam_count[fam] += 1

    # Step 3: Output
    filtered_bp = sum(fam_bp.values())
    filtered_count = sum(fam_count.values())

    # Per-family table
    col_fam = max((len(f) for f in fam_bp), default=10)
    col_fam = max(col_fam, len("Family"))
    header = f"{'Family':<{col_fam}}  {'Copies':>8}  {'Coverage_bp':>12}  {'Coverage_Mb':>12}"
    print(header)
    print("-" * len(header))
    for fam in sorted(fam_bp, key=lambda x: fam_bp[x], reverse=True):
        print(f"{fam:<{col_fam}}  {fam_count[fam]:>8}  {fam_bp[fam]:>12}  {fam_bp[fam]/1e6:>12.3f}")

    # Summary
    print("-" * len(header))
    print(f"{'Total':<{col_fam}}  {filtered_count:>8}  {filtered_bp:>12}  {filtered_bp/1e6:>12.3f}")
    print()
    print(f"Filtered families:    {len(families)}")
    print(f"Matched families:     {len(fam_bp)}")
    print(f"Total annotation bp:  {total_anno_bp}")
    print(f"Filtered coverage bp: {filtered_bp}")
    print(f"Fraction of anno:     {filtered_bp/total_anno_bp:.4f}" if total_anno_bp > 0 else "Fraction of anno:     NA")

if __name__ == "__main__":
    main()
