#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=02:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing


#SBATCH --job-name='Module1'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#module use /opt/rit/modules
#export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

#module load python/3

genomeFile=$1 # "/work/LAS/thomasp-lab/weijia/research/test1.9/Genome/LNNJ01.1.fsa_nt"
genomeName=$2 #"LNNJ"
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.9"
t=16
dir="/work/LAS/thomasp-lab/weijia/research/test1.9/Module1"

mkdir $genomeName
mkdir temp
echo "Module 1, Step 1: Blast Genome against Reference Library"

python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
cp $genomeName/*blast* temp/

echo "Module 1, Step 2: Select 100% coverage entries from Blast results"
python3 $path/Module1/Fullcov.py  -name $genomeName -p $path -d $dir

echo "Module1, Step 3: Making blastDB and get candidate sequences"

python3 $path/Module1/GetSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module1, Step 4: Check TIR and TSD"

python3 $path/Module1/CheckTIRTSD.py -name $genomeName -p $path -t $t -d $dir

echo "Module1, Step 5: Write to Gff3"

python3 $path/Module1/WriteToGff_M1.py -name $genomeName -p $path -t $t -d $dir

echo "Module1, Step 6: Check Low Complexity"

python3 $path/Module1/Lowcomp_M1.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
