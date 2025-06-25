#!/usr/bin/env python3
"""
Parallel HelitronScanner wrapper in Python.

### Original author: Michelle Stitzer, Apr 11, 2018
### Modifier: Shujun Ou (shujun.ou.1@gmail.com), May 1, 2019
### Revised: Tianyu Lu (tianyu@lu.fm), July 26, 2024
### Parallelized: Kenji Gerhardt (kenji.gerhardt@gmail.com), June 24, 2025

Usage:
  python helitron_scanner_parallel.py genome.fasta [--cpu N] [--hsdir PATH]

This script splits the input genome into chunks, runs HelitronScanner in parallel,
combines and adjusts coordinates, then formats and cleans up.
"""
import os
import sys
import re
import argparse
import subprocess
import tempfile
import shutil
from multiprocessing import Pool

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parallel HelitronScanner wrapper"
    )
    parser.add_argument(
        "--genome", help="Path to genome fasta file"
    )
    parser.add_argument(
        "--cpu", type=int, default=4,
        help="Number of parallel jobs (default: 4)"
    )
    parser.add_argument(
        "--hsdir", default=None,
        help="Path to HelitronScanner jar or directory (overrides HELITRONSCANNER_DIR env)"
    )
    return parser.parse_args()

# Determine absolute path to HelitronScanner directory.
def resolve_hsdir(arg_path=None):
    def normalize(p):
        # Expand user (~), vars, and strip quotes
        p = os.path.expanduser(os.path.expandvars(p.strip('"\'')))
        # Resolve symlinks to get actual file or directory
        return os.path.realpath(p)

    def verify_jar(path):
        jar = os.path.join(path, 'HelitronScanner.jar')
        if not os.path.isfile(jar):
            raise FileNotFoundError(f"HelitronScanner.jar not found in: {path}")
        return path

    # 1. CLI argument
    if arg_path:
        raw = arg_path
        p_real = normalize(raw)
        # If it's a directory, check for binary inside
        if os.path.isdir(p_real):
            binpath = os.path.join(p_real, 'HelitronScanner')
            jarpath = os.path.join(p_real, 'HelitronScanner.jar')
            if os.path.isfile(binpath):
                # follow the binary symlink
                p_real = os.path.realpath(binpath)
                p_dir = os.path.dirname(p_real)
            elif os.path.isfile(jarpath):
                p_dir = p_real
            else:
                raise FileNotFoundError(f"Neither HelitronScanner binary nor jar found in: {p_real}")
        else:
            # provided a file: expect binary
            if os.path.basename(p_real) == 'HelitronScanner':
                # resolve binary link
                p_dir = os.path.dirname(p_real)
            else:
                raise FileNotFoundError(f"Provided path is not a HelitronScanner directory or binary: {raw}")
        return verify_jar(p_dir)

    # 2. Find in PATH
    exe_name = 'HelitronScanner'
    exe_path = shutil.which(exe_name)
    if exe_path:
        p_real = normalize(exe_path)
        p_dir = os.path.dirname(p_real)
        return verify_jar(p_dir)

    # 3. Not found
    raise FileNotFoundError(
        "Could not locate HelitronScanner. "
        "Please specify with --hsdir or ensure the executable is in your PATH."
    )

def split_genome(genome_path, out_dir, chunk_size=10_000_000, overlap=50_000):
    """
    Split a fasta genome into chunks of given size with overlap.
    """
    def parse_genome(fasta):
        seqid = None
        buf = []
        with open(fasta) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    if seqid is not None:
                        yield seqid, ''.join(buf)
                    seqid = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line)
        if seqid is not None:
            yield seqid, ''.join(buf)

    os.makedirs(out_dir, exist_ok=True)
    for rec, seq in parse_genome(genome_path):
        length = len(seq)
        start, end = 0, 0 #initialize
        while end < length: #simplify loop check
            end = min(start + chunk_size, length)
            chunk_id = f"{rec}_{start}_{end}"
            with open(os.path.join(out_dir, f"{chunk_id}.fasta"), 'w') as out:
                out.write(f">{chunk_id}\n{seq[start:end]}\n")
            start = end - overlap

def process_chunk(args):
    """
    Run HelitronScanner commands on a single chunk.
    Returns (chunk_path, success, error_message).
    """
    chunk_path, hsdir = args
    prefix = os.path.splitext(chunk_path)[0]
    jar = os.path.join(hsdir, 'HelitronScanner.jar')
    head_lcvs = os.path.join(hsdir, 'TrainingSet', 'head.lcvs')
    tail_lcvs = os.path.join(hsdir, 'TrainingSet', 'tail.lcvs')
    cmds = []
    # direct orientation
    cmds.extend([
        ["java", "-Xmx2g", "-jar", jar, "scanHead", "-lcv_filepath", head_lcvs, "-g", chunk_path, "-buffer_size", "0", "-output", f"{prefix}.head"],
        ["java", "-Xmx2g", "-jar", jar, "scanTail", "-lcv_filepath", tail_lcvs, "-g", chunk_path, "-buffer_size", "0", "-output", f"{prefix}.tail"],
        ["java", "-Xmx2g", "-jar", jar, "pairends", "-head_score", f"{prefix}.head", "-tail_score", f"{prefix}.tail", "-output", f"{prefix}.pairends"],
        ["java", "-Xmx2g", "-jar", jar, "draw", "-pscore", f"{prefix}.pairends", "-g", chunk_path, "-output", f"{prefix}.draw", "-pure_helitron"]
    ])
    # reverse complement
    cmds.extend([
        ["java", "-Xmx2g", "-jar", jar, "scanHead", "-lcv_filepath", head_lcvs, "-g", chunk_path, "-buffer_size", "0", "--rc", "-output", f"{prefix}.rc.head"],
        ["java", "-Xmx2g", "-jar", jar, "scanTail", "-lcv_filepath", tail_lcvs, "-g", chunk_path, "-buffer_size", "0", "--rc", "-output", f"{prefix}.rc.tail"],
        ["java", "-Xmx2g", "-jar", jar, "pairends", "-head_score", f"{prefix}.rc.head", "-tail_score", f"{prefix}.rc.tail", "--rc", "-output", f"{prefix}.rc.pairends"],
        ["java", "-Xmx2g", "-jar", jar, "draw", "-pscore", f"{prefix}.rc.pairends", "-g", chunk_path, "-output", f"{prefix}.draw.rc", "-pure_helitron"]
    ])
    for c in cmds:
        try:
            quiet = subprocess.run(c, check=True)          
        except subprocess.CalledProcessError as e:
            return (chunk_path, False, str(e))
            
    try:
        #Clean as we go
        os.remove(f"{prefix}.head")
        os.remove(f"{prefix}.tail")
        os.remove(f"{prefix}.pairends")
        os.remove(f"{prefix}.rc.head")
        os.remove(f"{prefix}.rc.tail")
        os.remove(f"{prefix}.rc.pairends")
        os.remove(chunk_path)
    except: # if this fails, shutil remove will work later.
        pass
            
    return (chunk_path, True, '')

#Functions for combining chunks from parallel processing into a final output file.
def adjust_file(file, fmt_match):
    with open(file) as fh:
        for line in fh:
            if line.startswith('>'):
                matches = fmt_match.findall(line)[0]
                
                start = int(matches[0])
                
                pos1 = str(int(matches[2]) + start)
                pos2 = str(int(matches[3]) + start)
                
                line = line.replace('_'+matches[0]+'_', '_')
                line = line.replace('_'+matches[1]+'_', '_')
                line = line.replace('_'+matches[2]+'-', '_'+pos1+'-')
                line = line.replace('-'+matches[3]+' [', '-'+pos2+' [')
                
                segs = ''
                if len(matches[4]) > 0:
                    segs = matches[4].split()
                    
                    updated = []
                    for s in segs:
                        subsegs = s.split(':')
                        subsegs[0] = str(int(subsegs[0]) + start)
                        subsegs = ':'.join(subsegs)
                        updated.append(subsegs)
                    
                    segs = ' '.join(updated)
                    
                    line = line.replace(matches[4], segs)
                
                yield line, True
            else:
                yield line, False
        
def collect_inputs(tmpdir):
    forward = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith('.draw.hel.fa')]
    reverse = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith('.draw.rc.hel.fa')]

    return forward, reverse
             
def build_output(outfile, input_files):
    matcher = re.compile(r'_(\d+)_(\d+)_#SUB_(\d+)-(\d+) \[.+Multi_5\'_ends:(.+)?')
    
    headers = set()
    print_sequence = True
    with open(outfile, 'w') as out:
        for f in input_files:
            for line, isheader in adjust_file(f, matcher):
                if isheader:
                    simpleheader = line.split(';')[0] #just in case the ends are different; report only one
                    if not simpleheader in headers:
                        headers.add(simpleheader)
                        print_sequence = True
                        out.write(line)
                    else:
                        print_sequence = False #This helitron was in an overlap
                else:
                    if print_sequence:
                        out.write(line)
    
def combine_and_adjust(genome, tempdir):
    forward_in, reverse_in = collect_inputs(tempdir)
    forward_out, reverse_out = f"{genome}.HelitronScanner.draw.hel.fa", f"{genome}.HelitronScanner.draw.rc.hel.fa"
    
    build_output(forward_out, forward_in)
    build_output(reverse_out, reverse_in)
    
def main():
    args = parse_args()
    genome = args.genome
    cpu = args.cpu
    hsdir = resolve_hsdir(args.hsdir)
    print(f"Using HelitronScanner directory: {hsdir}")

    script_dir = os.path.dirname(os.path.realpath(__file__))
    cwd = os.getcwd()
    tmpdir = tempfile.mkdtemp(dir=cwd)
    print(f"Created temp directory: {tmpdir}")

    print("Splitting genome into chunks...")
    split_genome(genome, tmpdir)

    print("Processing chunks using Python multiprocessing...")
    chunks = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith('.fasta')]
    if not chunks:
        sys.exit("Error: no chunk files to process")
    with Pool(processes=cpu) as pool:
        for chunk, ok, err in pool.imap_unordered(process_chunk, [(c, hsdir) for c in chunks]):
            if ok:
                print(f"Successfully processed {chunk}")
            else:
                print(f"ERROR processing {chunk}: {err}", file=sys.stderr)
                pool.terminate()
                shutil.rmtree(tmpdir)
                sys.exit(1)

    print("Combining results and adjusting coordinates...")
    combine_and_adjust(genome, tmpdir)

    print("Cleaning up temporary files...")
    shutil.rmtree(tmpdir)

    print("HelitronScanner parallel processing completed successfully!")


if __name__ == '__main__':
    main()
