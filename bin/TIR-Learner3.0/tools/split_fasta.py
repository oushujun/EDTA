#!/usr/bin/env python3
# Split fasta tool
# Tianyu Lu (tlu83@wisc.edu)
# 2024-01-20

import argparse
from Bio import SeqIO
import os
from math import ceil
import random


def write_batch(base_filename, batch, file_suffix, file_number):
    filename = f"{base_filename}_{file_suffix}_{file_number}.fa"
    with open(filename, "w") as handle:
        SeqIO.write(batch, handle, "fasta")


def split_randomly(records):
    random.shuffle(records)
    return records


def split_by_seq_num(fasta_file, split_seq_num, randomize=False):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if randomize:
        records = split_randomly(records)
    base_filename = os.path.splitext(fasta_file)[0]
    for i in range(0, len(records), split_seq_num):
        write_batch(base_filename, records[i:i + split_seq_num], f"seq_num_{split_seq_num}", i // split_seq_num + 1)


def split_by_seq_len(fasta_file, split_seq_len, randomize=False):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if randomize:
        records = split_randomly(records)
    base_filename = os.path.splitext(fasta_file)[0]
    batch, current_length, file_number = [], 0, 1
    for record in records:
        if current_length + len(record.seq) > split_seq_len:
            write_batch(base_filename, batch, split_seq_len, file_number)
            batch, current_length = [record], len(record.seq)
            file_number += 1
        else:
            batch.append(record)
            current_length += len(record.seq)
    if batch:
        write_batch(base_filename, batch, f"seq_len_{split_seq_len}", file_number)


def split_by_file_num(fasta_file, split_file_num, randomize=False):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if randomize:
        records = split_randomly(records)
    records_per_file = ceil(len(records) / split_file_num)
    base_filename = os.path.splitext(fasta_file)[0]
    for i in range(0, len(records), records_per_file):
        write_batch(base_filename, records[i:i + records_per_file], f"file_num_{split_file_num}",
                    i // records_per_file + 1)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta_file", type=str, help="Path to the FASTA file", required=True)
    parser.add_argument("-sn", "--split_seq_num", type=int, help="Number of sequences per split file", default=None)
    parser.add_argument("-sl", "--split_seq_len", type=int, help="Maximum total sequence length per split file",
                        default=None)
    parser.add_argument("-fn", "--split_file_num", type=int, help="Number of files to split the FASTA into",
                        default=None)
    parser.add_argument("-r", "--random", action="store_true", help="Randomize sequence order before splitting")
    return parser.parse_args()


def main(fasta_file=None, split_seq_num=None, split_seq_len=None, split_file_num=None, randomize=False):
    if fasta_file is None:
        # Command line execution
        args = parse_arguments()
        fasta_file = args.fasta_file
        split_seq_num = args.split_seq_num
        split_seq_len = args.split_seq_len
        split_file_num = args.split_file_num
        randomize = args.random

    if split_seq_num:
        split_by_seq_num(fasta_file, split_seq_num, randomize)
    elif split_seq_len:
        split_by_seq_len(fasta_file, split_seq_len, randomize)
    elif split_file_num:
        split_by_file_num(fasta_file, split_file_num, randomize)
    else:
        raise SystemExit("No valid splitting criteria provided. Please specify either split_seq_num, split_seq_len, "
                         "or split_file_num.")


if __name__ == "__main__":
    main()
