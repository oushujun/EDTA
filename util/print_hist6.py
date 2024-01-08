#!/usr/bin/env python3

# Print a histogram on terminal from input data.
# Shujun Ou (shujun.ou.1@gmail.com), facilitated by ChatGPT
# Usage:
# python print_hist6.py -col 6 -binsize 0.2 -file Rice_MSU7.fasta.mod.Helitron.raw.fa.HQ-Rice_MSU7.fasta.mod.LTR.intact.fa.stat|less
# awk '{print $6}' Rice_MSU7.fasta.mod.Helitron.raw.fa.HQ-Rice_MSU7.fasta.mod.LTR.intact.fa.stat | python print_hist6.py -binsize 0.2 | less


import numpy as np
import argparse
import sys

def print_histogram(data, bin_size, max_bin, symbol='*'):
    # Define bin edges with given bin size and a final bin for values > max_bin
    bin_edges = np.arange(0, max_bin + bin_size, bin_size)
    bin_edges = np.append(bin_edges, float('inf'))  # Add an extra bin for values > max_bin

    # Calculate histogram
    counts, _ = np.histogram(data, bins=bin_edges)

    # Display histogram
    for i in range(len(counts)):
        if bin_edges[i] >= max_bin:
            bin_range = f"> {max_bin}"
        else:
            bin_range = f"{bin_edges[i]:.2f} - {bin_edges[i+1]:.2f}"
        print(f"Bin {bin_range}: {symbol * counts[i]}")

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description='Print a histogram on terminal from input data.')
    parser.add_argument('-binsize', type=float, default=0.5, help='Size of each bin. Default is 0.5.')
    parser.add_argument('-binmax', type=float, default=20, help='Maximum value for merging bins. Default is 20.')
    parser.add_argument('-col', type=int, default=1, help='Column number for data input (1-indexed). Default is 1.')
    parser.add_argument('-file', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Data input file. If not provided, reads from stdin.')
    args = parser.parse_args()

    # Read data, skipping non-numeric lines
    input_data = []
    for line in args.file:
        try:
            # Split the line into columns and select the specified column
            col_data = line.strip().split()[args.col - 1]  # Adjust for 0-indexed array
            num = float(col_data)
            input_data.append(num)
        except (ValueError, IndexError):
            continue  # Skip lines with non-numeric data or insufficient columns

    # Convert input data to a numpy array
    data = np.array(input_data)

    # Print histogram
    print_histogram(data, args.binsize, args.binmax)

if __name__ == "__main__":
    main()

