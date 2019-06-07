#!/bin/bash
# Contributors: Shujun Ou (shujun.ou.1@gmail.com); Weijia Su (weijia@iastate.edu)
# 06/07/2019

####################
### Installation ###
####################

# conda create -n TIR-Learner
# conda activate TIR-Learner
# conda install -y -c anaconda scikit-learn=0.19.0
# conda install -y -c anaconda biopython pandas glob2
# conda install -y -c conda-forge multiprocess regex
# conda install -y -c bioconda blast


#############
### Usage ###
#############

# sh TIR-Learner.sh genome.fa Index_name CPU_num


## pre-set parameters
genomeFile=$1 #the genome file
genomeName=$2 #base name for all outputs
path=$(dirname "$0") #program path
dir=$(pwd) #current work directory
grfp="~/las/bin/GenericRepeatFinder/bin/"
t=16 #CPUs

## allow user to specify CPU number
if [ ! -z "$3" ];
        then t=$3
fi


####################
### Preparations ###
####################

#mkdir temp
mkdir $genomeName
cd $genomeName

if [ ! 0 ]; then #debug
echo test
fi

## prepare blast databases
ln -s ../$genomeFile $genomeFile
for i in A C H M T; do
	cp $path/RefLib/Rice_DT${i}_RefLib Rice_DT${i}_RefLib
	makeblastdb -in Rice_DT${i}_RefLib -dbtype nucl -parse_seqids
	done

## Check sequence file
python3 $path/pre.py -g $genomeFile -name $genomeName


################
### Module 1 ###
################

echo "Module 1, Step 1: Blast Genome against Reference Library"
python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
#cp $genomeName/*blast* temp/

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


################
### Module 2 ###
################

echo "Module 2, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir -grfp $grfp

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module2/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 3 : GRF result blast reference sequences"
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 4: Get sequences from 80% similarity"
python3 $path/Module2/GetSeq_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 5: Check TIR and TSD"
python3 $path/Module2/CheckTIRTSD_M2.py -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 6: Write to gff3"
python3 $path/Module2/WriteToGff_M2.py -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 7: Remove Low complexity"
python3 $path/Module2/Lowcomp_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir


################
### Module 3 ###
################

echo "Module 3, Step 1: Get dataset"
python3 $path/Module3/getDataset.py -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 2: ML prediction"
python3 $path/Module3/ML_Ensemble.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 3: Check TIR/TSD"
python3 $path/Module3/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 4: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 5: Remove Low complexity"
python3 $path/Module3/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir


######################
### Postprocessing ###
######################

echo "Get Final GFF" 
python3 $path/CombineAll.py -name $genomeName -p $path -t $t -d $dir

echo "Get fasta file"
python3 $path/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

cd ..
mkdir TIR-Learner-Result
mv $genomeName/${genomeName}_combine/*FinalAnn.gff3 TIR-Learner-Result/
mv $genomeName/${genomeName}_combine/*FinalAnn.fa TIR-Learner-Result/

