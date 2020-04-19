#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=05:20:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL



#SBATCH --job-name='GT'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


module use /opt/rit/modules
module load genometools
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin
module load python/3
genomeName="Rice"
genome_file="/work/LAS/thomasp-lab/weijia/research/Rice/Genome/Rice.fa"
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.16"
species="Rice"
t=1
dir=$(pwd)
#echo "Start indexing"
#gt suffixerator -db $name -indexname $name -tis -suf -lcp -des -ssp -sds -dna -mirrored
#echo "Finish indexing"
##/usr/bin/time -v gt tirvish -index "Rice.fa" -mintirlen 10 -mintirdist 50 -similar 80.00 -seqids "yes"
#/usr/bin/time -v gt tirvish -index "Rice.fa" -seed 20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist 20000 -similar 80 -mintsd 2 -maxtsd 11 -vic 13 -seqids 'yes' > Rice_TIRvish.gff
#python3 ProcessTIRvish.py -g $genome_file -name $name -p $path -t 1 -d $path 
#python3 GetAllSeq_1.py
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir -s $species
cp $genomeName/*80 temp/ 
