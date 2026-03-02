#!/usr/bin/env python3
import sys
from collections import defaultdict
from statistics import mean

# Calculate TE annotation consistency per family
# Elena Zhu, Shujun Ou (shujun.ou.1@gmail.com), and ChatGPT
# v0.1: 08/15/2025
# v0.2: 01/04/2026

if len(sys.argv) != 3:
    print("Usage: python TE_consistency.py <TEanno.out> <stat_file>")
    sys.exit(1)

te_file = sys.argv[1]
stat_file = sys.argv[2]

# Step 1: Parse .stat file from cleanup_nested.pl
cleaned_cons = defaultdict(int)
discarded_cons = defaultdict(int)
cleaned_incons = defaultdict(int)
discarded_incons = defaultdict(int)
total = defaultdict(int)
ref_keys = set()

with open(stat_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # Determine status
        if "Cleaned." in line:
            status = "cleaned"
        elif "Discarded." in line:
            status = "discarded"
        else:
            continue
        # Extract copy type
        parts = line.split()
        if not parts:
            continue
        cp = parts[0].split("|")
        if len(cp) < 2:
            continue
        copy_type = cp[1]
        # Extract reference coordinate/type
        ref_str = ""
        if "by " in line:
            ref_str = line.split("by ")[1].split(";")[0]
        elif "of " in line:
            ref_str = line.split("of ")[1].split(";")[0]
        else:
            continue
        rf = ref_str.split("|")
        if len(rf) < 2:
            continue
        full_ref, ref_type = rf[0], rf[1]
        ref_keys.add(full_ref)
        total[full_ref] += 1

        if copy_type == ref_type:
            if status == "cleaned":
                cleaned_cons[full_ref] += 1
            else:
                discarded_cons[full_ref] += 1
        else:
            if status == "cleaned":
                cleaned_incons[full_ref] += 1
            else:
                discarded_incons[full_ref] += 1

# Step 2: Compute consistency values (at the 0–1 scale)
def ratio(n, total):
    return (n / total) if total > 0 else 0.0

consistency_table = {}
for ref in sorted(ref_keys):
    tot = total[ref]
    cc, dc = cleaned_cons[ref], discarded_cons[ref]
    ci, di = cleaned_incons[ref], discarded_incons[ref]
    total_consistency = ratio(cc, tot) + ratio(dc, tot)
    consistency_table[ref] = {
        "Total": tot,
        "Cln_Cons": cc,
        "Dis_Cons": dc,
        "Cln_Incons": ci,
        "Dis_Incons": di,
        "Consistency": total_consistency
    }

# Step 3: Parse TEanno.out to count TE families
family_counts = defaultdict(int)
copies = defaultdict(list)

with open(te_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 11:
            continue
        key = f"{parts[9]}#{parts[10]}"
        family_counts[key] += 1

# Step 4: Match coordinates to consistency stats
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
        if coord not in consistency_table:
            continue
        data = consistency_table[coord]
        consistency = round(data["Consistency"], 3)
        copies[fam].append({
            "Coordinate": coord,
            "Consistency": consistency,
            "Total": data["Total"],
            "Cln_Cons": data["Cln_Cons"],
            "Dis_Cons": data["Dis_Cons"],
            "Cln_Incons": data["Cln_Incons"],
            "Dis_Incons": data["Dis_Incons"]
        })

# Step 5: Format output (at 0-1 scale)
column_names = ["Fam/Ind", "Coordinate", "Consistency", "Total",
                "Cln_Cons", "Dis_Cons", "Cln_Incons", "Dis_Incons"]
col_widths = {name: len(name) for name in column_names}

# Adjust widths dynamically
for fam, entries in copies.items():
    for e in entries:
        vals = ["Ind", e["Coordinate"], f"{e['Consistency']:.3f}", str(e["Total"]),
                str(e["Cln_Cons"]), str(e["Dis_Cons"]),
                str(e["Cln_Incons"]), str(e["Dis_Incons"])]
        for i, name in enumerate(column_names):
            col_widths[name] = max(col_widths[name], len(vals[i]))

# Print header
header = "  ".join(f"{n:<{col_widths[n]}}" for n in column_names)
print(header)
print("-" * sum(col_widths.values()))

# Output rows per family
for fam, entries in sorted(copies.items()):
    percentages = [e["Consistency"] for e in entries if isinstance(e["Consistency"], (int, float))]
    totals = {k: 0 for k in ["Total", "Cln_Cons", "Dis_Cons", "Cln_Incons", "Dis_Incons"]}
    for e in entries:
        for key in totals:
            totals[key] += float(e[key])

    avg_cons = f"{mean(percentages):.3f}" if percentages else "NA"
    fam_row = ["Fam", fam, avg_cons,
               f"{round(totals['Total'], 1)}", f"{round(totals['Cln_Cons'], 1)}",
               f"{round(totals['Dis_Cons'], 1)}", f"{round(totals['Cln_Incons'], 1)}",
               f"{round(totals['Dis_Incons'], 1)}"]
    print("  ".join(f"{v:<{col_widths[column_names[i]]}}" for i, v in enumerate(fam_row)))

    for e in entries:
        row = ["Ind", e["Coordinate"], f"{e['Consistency']:.3f}", str(e["Total"]),
               str(e["Cln_Cons"]), str(e["Dis_Cons"]),
               str(e["Cln_Incons"]), str(e["Dis_Incons"])]
        print("  ".join(f"{v:<{col_widths[column_names[i]]}}" for i, v in enumerate(row)))

