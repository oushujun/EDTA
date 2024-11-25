# convert minimap2 paf file to blastn format 6 files
#   - each output contains ≤ 50M hits (~5Gb)
#   - multithreading acclerated
#   - keep the sequence input order in the output file(s)
# contributors: Herui Liao, Chris Benson, Shujun Ou

import pandas as pd
import subprocess
import os.path
import shlex
import sys
import math
import re
import concurrent.futures
import tempfile


# matching patterns
nonmatch_pattern = re.compile(r"NM:i:(\d+)")
cigar_pattern = re.compile(r"cg:Z:([A-Za-z0-9]+)")

class QualityCalculations:
    def __init__(self, genome_size=5000000000):
        self._k = 0.1
        self._lambda = 1.58
        self.genome_size = genome_size

    def calc_bitscore(self, alen, nonmatch):
        score = alen - 2 * nonmatch
        return (score * self._lambda - math.log(self._k)) / math.log(2.0)

    def calc_evalue(self, alen, nonmatch):
        score = max(0, alen - 2 * nonmatch)
        return self._k * alen * self.genome_size * math.exp(-self._lambda * score)

    def calc_gap_openings(self, cigar):
        go = 0
        for char in cigar:
            if char == "I" or char == "D":
                go += 1
        return go

# Function to match pattern in either `cg` or `extra` column
def find_cigar(row, qc):
    try:
        match = re.search(cigar_pattern, row['cg'])
    except:
        match = False
    if match:
        return qc.calc_gap_openings(match.group(1))
    try:
        match = re.search(cigar_pattern, row['extra'])
    except:
        match=False
    if match:
        return qc.calc_gap_openings(match.group(1))
    return 0

def process_chunk(df, qc):
    try:
        df["nonmatch"] = df.NM.map(
            lambda x: int(re.search(nonmatch_pattern, x).group(1))
        )
    except:
        df["nonmatch"]=0

    df['gap_openings'] = df.apply(lambda row: find_cigar(row, qc), axis=1)

    df["bitscore"] = [
        qc.calc_bitscore(a, n) for a, n in zip(df["alen"], df["nonmatch"])
    ]
    df["evalue"] = [qc.calc_evalue(a, n) for a, n in zip(df["alen"], df["nonmatch"])]
    df["percent_ident"] = [
        (nmatch / a) * 100 for a, nmatch in zip(df["alen"], df["nmatch"])
    ]
    df = df.round({"bitscore": 3, "percent_ident": 3})
    m = df["strand"] == "-"
    df.loc[m, ["tstart", "tend"]] = (df.loc[m, ["tend", "tstart"]].values)

    blast = df.loc[
        :,
        [
            "qname",
            "tname",
            "percent_ident",
            "alen",
            "nonmatch",
            "gap_openings",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "evalue",
            "bitscore",
        ],
    ]
    blast["qstart"] = blast["qstart"] + 1
    blast.loc[~m, "tstart"] = blast.loc[~m, "tstart"] + 1
    blast.loc[m, "tend"] = blast.loc[m, "tend"] + 1
    blast.loc[:, "tstart"] = blast.loc[:, "tstart"].astype(int)
    blast.loc[:, "tend"] = blast.loc[:, "tend"].astype(int)
    return blast

def parallel_process_chunk(chunk, chunk_idx, qc, output_dir):
    processed_chunk = process_chunk(chunk, qc)
    temp_file = os.path.join(output_dir, f"temp_{chunk_idx}.txt")
    processed_chunk.to_csv(temp_file, sep="\t", header=False, index=False)
    return temp_file

def parallel_process_file(input_file, odir, chunk_size=100000, max_lines=100000000, n_threads=5):
    qc = QualityCalculations()
    output_dir = tempfile.mkdtemp()
    temp_files = []

    # Using ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for chunk_idx, chunk in enumerate(pd.read_csv(
            input_file,
            engine='c',
            sep="\t",
            header=None,
            names=[
                "qname",
                "qlen",
                "qstart",
                "qend",
                "strand",
                "tname",
                "tlen",
                "tstart",
                "tend",
                "nmatch",
                "alen",
                "mapq",
                'NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 'de', 'rl', 'cg', 'extra',
            ],
            dtype={'cg':'str', 'extra':'str'},
            chunksize=chunk_size)
            ):
            futures.append(executor.submit(parallel_process_chunk, chunk, chunk_idx, qc, output_dir))

        # Collecting the results (temp file paths)
        for future in concurrent.futures.as_completed(futures):
            temp_files.append(future.result())


    # Merging the temporary files to produce the final output in the same order as the input
    part_num = 1  # Suffix for output files if they exceed max lines
    line_count = 0  # Line counter
    output_file = os.path.join(odir, os.path.basename(paf_file.replace(".paf", ".out")))
    output_file_base = output_file  # Save the base output filename
    outfile = open(output_file, "w")  # Open the first output file

    # Sort the temp_files list by the chunk index to maintain input order
    temp_files_sorted = sorted(temp_files, key=lambda x: int(re.search(re.compile(r"temp_(\d+).txt"), x).group(1)))
    prev=''
    for temp_file in temp_files_sorted:
        with open(temp_file, "r") as infile:
            for line in infile:
                ele=line.split()
                if line_count >= max_lines:
                    # Close the current file and open a new one with an incremented suffix

                    if ele[0]==prev:
                        outfile.write(line)
                        continue

                    outfile.close()
                    part_num += 1
                    output_file = f"{output_file_base}_part{part_num}"
                    outfile = open(output_file, "w")
                    line_count = 0  # Reset line count for the new file
                outfile.write(line)
                line_count += 1
                prev=ele[0]
        os.remove(temp_file)

    outfile.close()  # Close the last file
    os.rmdir(output_dir)

if __name__ == "__main__":
    paf_file = sys.argv[1]
    odir = sys.argv[2]
    n_threads = int(sys.argv[4]) if '-t' in sys.argv else 5
    chunk_size=100000
    max_lines = 100000000 #~10G/output file
    #max_lines=10
    parallel_process_file(paf_file, odir, chunk_size, max_lines, n_threads)
