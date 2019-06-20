#!/bin/bash

## Modify this job script accordingly and submit with the command:
##    sbatch HPC.sbatch

#SBATCH --time=50:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --partition=freecomputing
#SBATCH --job-name='TIR-Learner'

## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

module load python/3
module load cd-hit
module load gcc

genome=/work/LAS/thomasp-lab/weijia/research/Rice/Genome/Rice.fa
threads=16

sh TIR-Learner.sh $genome $threads


