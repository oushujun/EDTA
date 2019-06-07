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


#SBATCH --job-name='post'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

module load python/3

genomeFile="/work/LAS/thomasp-lab/weijia/research/test1.9/Genome/LNNJ01.1.fsa_nt"
genomeName="LNNJ"
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.9"
dir="/work/LAS/thomasp-lab/weijia/research/test1.9"
t=16
mkdir $genomeName


echo "Get Final GFF" 
python3 $path/CombineAll.py -name $genomeName -p $path -t $t -d $dir  

mv *.gff3 $genomeName
rm *Low

echo "Get fasta file"
python3 $path/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir


mkdir TIR-Learner-Result
mv $genomeName/*FinalAnn.gff3 TIR-Learner-Result/
mv $genomeName/*FinalAnn.fa TIR-Learner-Result/
rm -r $genomeName"_combine"
rm -r $genomeName
for i in Module1 Module2 Module3; do rm -r $i/$genomeName ; done
for i in Module1 Module2 Module3; do rm -r $i/temp ; done

