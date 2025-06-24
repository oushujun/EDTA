#!/bin/bash -login
### This script was modified from https://github.com/mcstitzer/maize_v4_TE_annotation/blob/master/helitron/run_helitron_scanner.sh
### Original author: Michelle Stitzer, Apr 11, 2018
### Modifier: Shujun Ou (shujun.ou.1@gmail.com), May 1, 2019
### Revised: Tianyu Lu (tianyu@lu.fm), July 26, 2024
### Parallelized: Kenji Gerhardt (kenji.gerhardt@gmail.com), June 24, 2025
### specify the genome file
GENOME=$1
### the base path of this script
path=$(dirname "$0")
## where to find HelitronScanner.jar
#HSDIR=$path/../bin/HelitronScanner
# Find original directory of bash script, resolving symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
SOURCE=$3
echo $SOURCE
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
HSDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
### preset CPU and max memory
CPU=4
MEMGB=32 #Gb
### allow user to specify CPU number to run HelitronScanner
if [ ! -z "$2" ];
then CPU=$2
fi

dir=$(pwd)
# Create temporary directory
TEMP_DIR=$(mktemp -d -p $dir)
#mkdir TEMP_DIR
echo "Created temp directory: $TEMP_DIR"

# Split genome into chunks (max 10Mbp with 5kbp overlaps)
echo "Splitting genome into chunks..."
python3 -c "
import os
import sys
import re
#from Bio import SeqIO

genome = sys.argv[1]
out_dir = sys.argv[2]
chunk_size = 10000000
overlap = 50000

def parse_genome(file):
    genomes = {}
    seqid = None
    sequences = []
    with open(file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                sequences = ''.join(sequences)
                if seqid is not None:
                    yield seqid, sequences
                    
                seqid = line[1:].split()[0]
                sequences = []
            else:
                sequences.append(line)
    
    sequences = ''.join(sequences)
    return seqid, sequences

def sanitize(name):
    return re.sub(r'[^\w]', '_', name)
    
sequence_order = []
for record, sequence in parse_genome(genome):
    sequence_order.append(record)
    seq_len = len(sequence)
    start = 0
    chunk_idx = 0
    while start < seq_len:
        end = min(start + chunk_size, seq_len)
        chunk_id = f'{record}_{start}_{end}'
        with open(f'{out_dir}/{chunk_id}.fasta', 'w') as out:
            out.write(f'>{chunk_id}\n{str(sequence[start:end])}\n')    
        if end == seq_len:
            break
        start = end - overlap
        chunk_idx += 1
		
#with open(f'{out_dir}/seq_order.txt', 'w') as out:
#    for s in sequence_order:
#        print(s, file = out)

" "$GENOME" "$TEMP_DIR"

# Calculate memory per chunk (reserve 5GB for system)
CHUNK_MEM=$(( (MEMGB - 5) / (CPU * 2) ))g
echo "Using $CHUNK_MEM per chunk with $CPU parallel jobs"


# Use Python for multiprocessing
echo "Processing chunks using Python multiprocessing..."
python3 -c "
import os
import sys
import subprocess
from multiprocessing import Pool, cpu_count

def process_chunk(args):
    try:
        chunk = args[0]
        HSDIR = args[1]
        
        prefix = chunk.replace('.fasta', '')
        
        # Direct orientation processing
        cmd1 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar scanHead -lcv_filepath {HSDIR}/TrainingSet/head.lcvs -g {chunk} -buffer_size 0 -output {prefix}.head'
        cmd2 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar scanTail -lcv_filepath {HSDIR}/TrainingSet/tail.lcvs -g {chunk} -buffer_size 0 -output {prefix}.tail'
        cmd3 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar pairends -head_score {prefix}.head -tail_score {prefix}.tail -output {prefix}.pairends'
        cmd4 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar draw -pscore {prefix}.pairends -g {chunk} -output {prefix}.draw -pure_helitron'
        
        # Reverse complement processing
        cmd5 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar scanHead -lcv_filepath {HSDIR}/TrainingSet/head.lcvs -g {chunk} -buffer_size 0 --rc -output {prefix}.rc.head'
        cmd6 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar scanTail -lcv_filepath {HSDIR}/TrainingSet/tail.lcvs -g {chunk} -buffer_size 0 --rc -output {prefix}.rc.tail'
        cmd7 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar pairends -head_score {prefix}.rc.head -tail_score {prefix}.rc.tail --rc -output {prefix}.rc.pairends'
        cmd8 = f'java -Xmx2g -jar {HSDIR}/HelitronScanner.jar draw -pscore {prefix}.rc.pairends -g {chunk} -output {prefix}.draw.rc -pure_helitron'

        cmds = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8]
        cmds = [c.split() for c in cmds]
        
        for c in cmds:
            subprocess.run(c, check=True)
            
        return (chunk, True, '')
        
    except subprocess.CalledProcessError as e:
        return (chunk, False, str(e))

if __name__ == '__main__':
    # Get parameters from environment
    TEMPDIR = sys.argv[1]
    HSDIR = sys.argv[2]
    CPU = int(sys.argv[3])
    
    # Get chunks from command line arguments
    chunks = [os.path.join(TEMPDIR, f) for f in os.listdir(TEMPDIR) if f.endswith('.fasta')]
    if not chunks:
        sys.exit('No chunk files provided')
    
    print(f'Processing {len(chunks)} chunks with {CPU} workers')
    
    # Create process pool
    with Pool(processes=CPU) as pool:
        results = []
        # Prepare arguments for each task
        task_args = [(chunk, HSDIR,) for chunk in chunks]
        
        # Process chunks in parallel
        for result in pool.imap_unordered(process_chunk, task_args):
            chunk, success, error = result
            if success:
                print(f'Successfully processed {chunk}')
            else:
                print(f'ERROR processing {chunk}: {error}', file=sys.stderr)
                pool.terminate()  # Stop immediately on failure
                sys.exit(1)
" "$TEMP_DIR" "$HSDIR" "$CPU"

# Check Python exit status
if [ $? -ne 0 ]; then
    echo "Error: Python multiprocessing failed" >&2
    exit 1
fi

# Combine and adjust coordinates
echo "Combining results and adjusting coordinates..."
#> ${GENOME}.HelitronScanner.draw.hel.fa
#> ${GENOME}.HelitronScanner.draw.rc.hel.fa

python3 -c "
import os
import re
import sys

#Adjust start, end coordinates by chunk
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
        
#Write to output                
def build_output(outfile, input_files):
    matcher = re.compile(r'_(\d+)_(\d+)_#SUB_(\d+)-(\d+) \[.+Multi_5\'_ends:(.+)?')
    
    headers = {}
    print_sequence = True
    with open(outfile, 'w') as out:
        for f in input_files:
            for line, isheader in adjust_file(f, matcher):
                if isheader:
                    simpleheader = line.split(';')[0] #just in case the ends are different; report only one
                    if not simpleheader in headers:
                        headers[simpleheader] = 0
                        print_sequence = True
                        out.write(line)
                    else:
                        print_sequence = False #This helitron was in an overlap
                else:
                    if print_sequence:
                        out.write(line)

def collect_inputs(tmpdir):
    forward = [os.path.join(tmpdir, f) for f in os.listdir(tempdir) if f.endswith('.draw.hel.fa')]
    reverse = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith('.draw.rc.hel.fa')]

    return forward, reverse

tempdir = sys.argv[1]
forward_out = sys.argv[2]
reverse_out = sys.argv[3]

fin, rin = collect_inputs(tempdir)
build_output(forward_out, fin)
build_output(reverse_out, rin)
" "$TEMP_DIR" "$GENOME.HelitronScanner.draw.hel.fa" "$GENOME.HelitronScanner.draw.rc.hel.fa"

# Run the final formatting step
echo "Running final formatting..."
perl $path/format_helitronscanner_out.pl -genome $GENOME -sitefilter 1 -minscore 12 -keepshorter 1 -extout 0

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "HelitronScanner parallel processing completed successfully!"
