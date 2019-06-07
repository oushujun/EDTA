#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=01:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing


#SBATCH --job-name='Pre'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


#module use /opt/rit/modules
#export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

#module load python/3

genomeFile=$1 #"/work/LAS/thomasp-lab/weijia/research/test1.9/Genome/LNNJ01.1.fsa_nt"
genomeName=$2 #"LNNJ"
#path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.9"
path=$(dirname "$0") #$(pwd)

cp $path/*.sh .

python3 $path/pre.py -g $genomeFile -name $genomeName

for i in Module1 Module2 Module3; do mkdir $i;done
for i in Module1 Module2 Module3; do cp $path/$i/*.sh $i/;done

