#!/usr/bin/env python3
import sys
import re
import argparse
from collections import defaultdict

# Calculate TE annotation consistency per family
# Elena Zhu, Shujun Ou (shujun.ou.1@gmail.com), and ChatGPT
# v0.1: 08/15/2025
# v0.2: 01/04/2026
# v0.3: 03/10/2026 Add mincov filter, category breakdown (all/nested/redun), flag-based CLI

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate TE annotation consistency per family")
    parser.add_argument("-anno", required=True,
                        help="TE annotation in RepeatMasker .out format")
    parser.add_argument("-stat", required=True,
                        help=".stat file from cleanup_nested.pl")
    parser.add_argument("-out", required=True,
                        help="Output file path")
    parser.add_argument("-mincov", type=float, default=0.95,
                        help="Minimum reciprocal coverage for redun category (default: 0.95)")
    return parser.parse_args()

def parse_stat_file(stat_file, mincov):
    """Parse .stat file and classify entries into nested/redun categories per reference coordinate."""
    # Per-ref counters: {ref_coord: {category: {"cons": N, "incons": N}}}
    stats = defaultdict(lambda: {
        "nested": {"cons": 0, "incons": 0},
        "redun":  {"cons": 0, "incons": 0},
    })
    ref_keys = set()

    with open(stat_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Determine status
            if "Cleaned." in line:
                status = "nested"
            elif "Discarded." in line:
                status = "redun"
            else:
                continue

            # Skip corrupted rows (valid rows have at least 10 fields)
            fields = line.split()
            if len(fields) < 10:
                continue

            # Extract copy type from subject ID (first field: id|type)
            cp = fields[0].split("|")
            if len(cp) < 2:
                continue
            copy_type = cp[1]

            # Extract reference ID and type
            if status == "redun":
                # Discarded: "by {ref_id}|{ref_type};"
                m = re.search(r'by\s+(\S+)\|(\S+?);', line)
            else:
                # Cleaned: "of {ref_id}|{ref_type};"
                m = re.search(r'of\s+(\S+)\|(\S+?);', line)
            if not m:
                continue
            full_ref, ref_type = m.group(1), m.group(2)

            # Extract qcov and scov
            if status == "redun":
                qcov_m = re.search(r'qcov:\s+([\d.]+)', line)
                scov_m = re.search(r'scov:\s+([\d.]+)', line)
            else:
                # Cleaned: qcov from "covering X.XXX of"
                qcov_m = re.search(r'covering\s+([\d.]+)\s+of', line)
                scov_m = re.search(r'scov:\s+([\d.]+)', line)

            # Validate extracted values (guard against corrupted lines from multithreading)
            try:
                qcov = float(qcov_m.group(1)) if qcov_m else 0
                scov = float(scov_m.group(1)) if scov_m else 0
            except ValueError:
                continue
            if qcov > 1.0 or scov > 1.0:
                continue

            # For redun: apply reciprocal coverage filter
            if status == "redun" and mincov > 0:
                if qcov > 0 and scov > 0:
                    if qcov < mincov or scov < mincov:
                        continue

            ref_keys.add(full_ref)
            is_consistent = (copy_type == ref_type)
            if is_consistent:
                stats[full_ref][status]["cons"] += 1
            else:
                stats[full_ref][status]["incons"] += 1

    return stats, ref_keys

def compute_rate(cons, incons):
    total = cons + incons
    return cons / total if total > 0 else None

def format_rate(rate):
    return f"{rate:.3f}" if rate is not None else "NA"

def parse_anno_file(te_file):
    """Parse TEanno.out and return {family: [coordinates]}."""
    copies = defaultdict(list)
    with open(te_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 11:
                continue
            fam = parts[9]
            chr_ = parts[4]
            try:
                start = int(parts[5])
                end = int(parts[6])
            except ValueError:
                continue
            coord = f"{chr_}:{start}..{end}"
            copies[fam].append(coord)
    return copies

def main():
    args = parse_args()

    # Step 1: Parse stat file
    stats, ref_keys = parse_stat_file(args.stat, args.mincov)

    # Step 2: Parse anno file and match to stat entries
    fam_copies = parse_anno_file(args.anno)

    # Step 3: Build per-family, per-coordinate data
    categories = ["all", "nested", "redun"]
    cat_cols = []
    for cat in categories:
        cat_cols += [f"{cat}_Total", f"{cat}_Cons", f"{cat}_Incons", f"{cat}_Rate"]

    column_names = ["Fam/Ind", "Coordinate"] + cat_cols

    # Collect rows
    all_rows = []  # list of (fam, rows_for_fam)
    for fam in sorted(fam_copies.keys()):
        ind_rows = []
        for coord in fam_copies[fam]:
            if coord not in stats:
                continue
            s = stats[coord]
            row = {"Coordinate": coord}
            for cat in categories:
                if cat == "all":
                    cons = s["nested"]["cons"] + s["redun"]["cons"]
                    incons = s["nested"]["incons"] + s["redun"]["incons"]
                else:
                    cons = s[cat]["cons"]
                    incons = s[cat]["incons"]
                total = cons + incons
                rate = compute_rate(cons, incons)
                row[f"{cat}_Total"] = total
                row[f"{cat}_Cons"] = cons
                row[f"{cat}_Incons"] = incons
                row[f"{cat}_Rate"] = rate
            ind_rows.append(row)

        if not ind_rows:
            continue

        # Family summary
        fam_summary = {"Coordinate": fam}
        for cat in categories:
            totals_sum = sum(r[f"{cat}_Total"] for r in ind_rows)
            cons_sum = sum(r[f"{cat}_Cons"] for r in ind_rows)
            incons_sum = sum(r[f"{cat}_Incons"] for r in ind_rows)
            fam_rate = cons_sum / totals_sum if totals_sum > 0 else None
            fam_summary[f"{cat}_Total"] = totals_sum
            fam_summary[f"{cat}_Cons"] = cons_sum
            fam_summary[f"{cat}_Incons"] = incons_sum
            fam_summary[f"{cat}_Rate"] = fam_rate

        all_rows.append((fam, fam_summary, ind_rows))

    # Step 4: Compute column widths
    col_widths = {name: len(name) for name in column_names}

    for fam, fam_summary, ind_rows in all_rows:
        for row_data in [fam_summary] + ind_rows:
            label = "Fam" if row_data is fam_summary else "Ind"
            vals = [label, row_data["Coordinate"]]
            for cat in categories:
                vals += [
                    str(row_data[f"{cat}_Total"]),
                    str(row_data[f"{cat}_Cons"]),
                    str(row_data[f"{cat}_Incons"]),
                    format_rate(row_data[f"{cat}_Rate"]),
                ]
            for i, name in enumerate(column_names):
                col_widths[name] = max(col_widths[name], len(vals[i]))

    # Step 5: Write output
    with open(args.out, 'w') as out:
        # Header
        header = "  ".join(f"{n:<{col_widths[n]}}" for n in column_names)
        out.write(header + "\n")
        out.write("-" * len(header) + "\n")

        for fam, fam_summary, ind_rows in all_rows:
            # Family row
            vals = ["Fam", fam_summary["Coordinate"]]
            for cat in categories:
                vals += [
                    str(fam_summary[f"{cat}_Total"]),
                    str(fam_summary[f"{cat}_Cons"]),
                    str(fam_summary[f"{cat}_Incons"]),
                    format_rate(fam_summary[f"{cat}_Rate"]),
                ]
            out.write("  ".join(
                f"{v:<{col_widths[column_names[i]]}}" for i, v in enumerate(vals)
            ) + "\n")

            # Individual rows
            for row_data in ind_rows:
                vals = ["Ind", row_data["Coordinate"]]
                for cat in categories:
                    vals += [
                        str(row_data[f"{cat}_Total"]),
                        str(row_data[f"{cat}_Cons"]),
                        str(row_data[f"{cat}_Incons"]),
                        format_rate(row_data[f"{cat}_Rate"]),
                    ]
                out.write("  ".join(
                    f"{v:<{col_widths[column_names[i]]}}" for i, v in enumerate(vals)
                ) + "\n")

if __name__ == "__main__":
    main()
