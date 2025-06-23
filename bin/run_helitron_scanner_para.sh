#!/bin/bash -login
### This script was modified from https://github.com/mcstitzer/maize_v4_TE_annotation/blob/master/helitron/run_helitron_scanner.sh
### Original author: Michelle Stitzer, Apr 11, 2018
### Modifier: Shujun Ou (shujun.ou.1@gmail.com), May 1, 2019
### Revised: Tianyu Lu (tianyu@lu.fm), July 26, 2024
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
MEMGB=150 #Gb
### allow user to specify CPU number to run HelitronScanner
if [ ! -z "$2" ];
then CPU=$2
fi

# Create temporary directory
TEMP_DIR=$(mktemp -d)
echo "Created temp directory: $TEMP_DIR"

# Split genome into chunks (max 10Mbp with 5kbp overlaps)
echo "Splitting genome into chunks..."
python3 -c "
import os
import sys
import re
from Bio import SeqIO

genome = sys.argv[1]
out_dir = sys.argv[2]
chunk_size = 10000000
overlap = 5000

def sanitize(name):
    return re.sub(r'[^\w]', '_', name)

for record in SeqIO.parse(genome, 'fasta'):
    seq_len = len(record)
    start = 0
    chunk_idx = 0
    while start < seq_len:
        end = min(start + chunk_size, seq_len)
        chunk_id = f'{sanitize(record.id)}_{start}_{end}'
        with open(f'{out_dir}/{chunk_id}.fasta', 'w') as out:
            out.write(f'>{chunk_id}\n{str(record.seq[start:end])}\n')
        start = end - overlap
        chunk_idx += 1
" "$GENOME" "$TEMP_DIR"

# Function to process a single chunk
process_chunk() {
    chunk=$1
    HSDIR=$2
    CHUNK_MEM=$3
    
    prefix="${chunk%.*}"
    
    # Direct orientation processing
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar scanHead -lcv_filepath ${HSDIR}/TrainingSet/head.lcvs -g $chunk -buffer_size 0 -output ${prefix}.head
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar scanTail -lcv_filepath ${HSDIR}/TrainingSet/tail.lcvs -g $chunk -buffer_size 0 -output ${prefix}.tail
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar pairends -head_score ${prefix}.head -tail_score ${prefix}.tail -output ${prefix}.pairends
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar draw -pscore ${prefix}.pairends -g $chunk -output ${prefix}.draw -pure_helitron
    
    # Reverse complement processing
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar scanHead -lcv_filepath ${HSDIR}/TrainingSet/head.lcvs -g $chunk -buffer_size 0 --rc -output ${prefix}.rc.head
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar scanTail -lcv_filepath ${HSDIR}/TrainingSet/tail.lcvs -g $chunk -buffer_size 0 --rc -output ${prefix}.rc.tail
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar pairends -head_score ${prefix}.rc.head -tail_score ${prefix}.rc.tail --rc -output ${prefix}.rc.pairends
    java -Xmx${CHUNK_MEM} -jar ${HSDIR}/HelitronScanner.jar draw -pscore ${prefix}.rc.pairends -g $chunk -output ${prefix}.draw.rc -pure_helitron
}
export -f process_chunk

# Calculate memory per chunk (reserve 5GB for system)
CHUNK_MEM=$(( (MEMGB - 5) / (CPU * 2) ))g
echo "Using $CHUNK_MEM per chunk with $CPU parallel jobs"

# Process chunks in parallel
find "$TEMP_DIR" -name "*.fasta" | parallel -j $CPU --halt now,fail=1 "process_chunk {} $HSDIR $CHUNK_MEM"

# Combine and adjust coordinates
echo "Combining results and adjusting coordinates..."
> ${GENOME}.HelitronScanner.draw.hel.fa
> ${GENOME}.HelitronScanner.draw.rc.hel.fa

for chunk in "$TEMP_DIR"/*.fasta; do
    prefix="${chunk%.*}"
    contig=$(basename "$prefix")
    
    # Extract coordinates from filename (contig_start_end)
    IFS='_' read -ra ADDR <<< "$contig"
    start_pos="${ADDR[-2]}"
    
    # Adjust direct orientation coordinates
    awk -v start_pos="$start_pos" '{
        if (/^>/) {
            split($0, parts, ":");
            new_start = parts[3] + start_pos;
            new_end = parts[4] + start_pos;
            printf ">%s_%s_%s:%s:%d:%d:%s:%s:%s\n", parts[1], parts[2], new_start, new_end, parts[5], parts[6], parts[7], parts[8], parts[9];
        } else {
            print $0;
        }
    }' "${prefix}.draw.hel.fa" >> ${GENOME}.HelitronScanner.draw.hel.fa
    
    # Adjust reverse complement coordinates
    awk -v start_pos="$start_pos" '{
        if (/^>/) {
            split($0, parts, ":");
            new_start = parts[3] + start_pos;
            new_end = parts[4] + start_pos;
            printf ">%s_%s_%s:%s:%d:%d:%s:%s:%s\n", parts[1], parts[2], new_start, new_end, parts[5], parts[6], parts[7], parts[8], parts[9];
        } else {
            print $0;
        }
    }' "${prefix}.draw.rc.hel.fa" >> ${GENOME}.HelitronScanner.draw.rc.hel.fa
done

# Run the final formatting step
echo "Running final formatting..."
perl $path/format_helitronscanner_out.pl -genome $GENOME -sitefilter 1 -minscore 12 -keepshorter 1 -extout 0

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "HelitronScanner parallel processing completed successfully!"
