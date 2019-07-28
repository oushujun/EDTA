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



#SBATCH --job-name='TIRLearner'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin
module load python/3

## pre-set parameters
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.17" #program path
dir=$(pwd) #current work directory
rawFile="/work/LAS/thomasp-lab/weijia/research/testData/ChrXXI.fasta" #the genome file with relative path
species="other" # One of the following "Maize", "Rice" or "others"
len=5000 # Maximum of element length used in GRF
t=16 # CPUs

## read parameters and help doc
helpFunction()
{
   echo ""
   echo "Usage: sh $0 -g genome.fa -s [Rice|Maize|others] -p [path] -l [int] -t [int]"
   echo -e "\t-g	The genome file in multi-FASTA format"
   echo -e "\t-s	Specify species: Rice, Maize, or others. Default: others"
   echo -e "\t-p	Specify the path to the GRF program"
   echo -e "\t-l	Maximum spacer length between the terminal repeat (used in GRF). Default: 5000 (bp)"
   echo -e "\t-t	Maximum thread number used in this program. Default: 16"
   echo ""
   exit 1 # Exit script after printing help
}

while getopts "g:s:p:l:t:" opt
do
   case ${opt} in
      g ) rawFile="$OPTARG" ;;
      s ) species="$OPTARG" ;;
      p ) grfp="$OPTARG" ;;
      l ) len="$OPTARG" ;;
      t ) t="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# read file name and path
genomeFile=$rawFile #the genome file with real path
genomeName="XII"


####################
### Preparations ###
####################

echo "############################################################ Pre-Processing ###########################################################"

# Check sequence file
python3 $path/pre.py -g $genomeFile -name $genomeName

## mkdir temp
for i in Module1 Module2 Module3; do [ -d $i ] || mkdir $i; done


################
### Module 1 ###
################

echo "############################################################ Module 1 Begin ###########################################################"

cd $dir"/Module1"
[ -d $genomeName ] || mkdir $genomeName
[ -d temp ] || mkdir temp

echo "Module 1, Step 1: Blast Genome against Reference Library"
python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1" -s $species
cp $genomeName/*blast* temp/

cd $genomeName
file=tem_blastRestul

if [ -s $file ]

then

cd $dir"/Module1"
echo "Module 1, Step 2: Select 100% coverage entries from Blast results"
python3 $path/Module1/Fullcov.py  -name $genomeName -p $path -d $dir"/Module1" -s $species

echo "Module 1, Step 3: Making blastDB and get candidate sequences"
python3 $path/Module1/GetSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module 1, Step 4: Check TIR and TSD"
python3 $path/Module1/CheckTIRTSD.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module 1, Step 5: Write to Gff3"
python3 $path/Module1/WriteToGff_M1.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module 1, Step 6: Check Low Complexity"
python3 $path/Module1/Lowcomp_M1.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"


echo "############################################################ Module 1 Finished ########################################################"


################
### Module 2 ###
################

echo "############################################################ Module 2 Begin ###########################################################"

cd $dir"/Module2/"
[ -d $genomeName ] || mkdir $genomeName
[ -d temp ] || mkdir temp

echo "Module 2, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
python3 $path/Module2/RunGT.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2" -l $len
cp -r $genomeName/$genomeName* temp/

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module2/ProcessTIRvish.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 3: GRF result blast reference sequences"
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir"/Module2" -s $species
cp $genomeName/*80 temp/

echo "Module 2, Step 4: Get sequences from 80% similarity"
python3 $path/Module2/GetSeq_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 5: Check TIR and TSD"
python3 $path/Module2/CheckTIRTSD_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 6: Write to gff3"
python3 $path/Module2/WriteToGff_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 7: Remove Low complexity"
python3 $path/Module2/Lowcomp_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"


echo "############################################################ Module 2 Finished ########################################################"


################
### Module 3 ###
################
#
echo "############################################################ Module 3 Begin ###########################################################"
#
cd $dir"/Module3/"
[ -d $genomeName ] || mkdir $genomeName
[ -d temp ] || mkdir temp
#
mv ../Module2/$genomeName/*_nonHomo.fa $genomeName/

echo "Module 3, Step 1: Get dataset"
python3 $path/Module3/getDataset.py -name $genomeName -p $path -t $t -d $dir"/Module3"
cp $genomeName/*.csv temp/
cp $genomeName/*nonHomo.fa temp/

echo "Module 3, Step 2: ML prediction"
python3 $path/Module3/ML_Ensemble.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3" -s $species
#
echo "Module 3, Step 3: Check TIR/TSD"
python3 $path/Module3/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"
#
echo "Module 3, Step 4: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"
#
echo "Module 3, Step 5: Remove Low complexity"
python3 $path/Module3/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3"
#
echo "############################################################ Module 3 Finished ########################################################"
#
#
######################
### Postprocessing ###
######################
#
echo "############################################################ Post Processing  #########################################################"

cd $dir
[ -d $genomeName ] || mkdir $genomeName

echo "Get Final GFF"
python3 $path/CombineAll.py -name $genomeName -p $path -t $t -d $dir

mv *.gff3 $genomeName
rm *Low

echo "Get fasta file"
python3 $path/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

[ -d "TIR-Learner-Result" ] || mkdir TIR-Learner-Result
mv $genomeName/*FinalAnn.gff3 TIR-Learner-Result/
mv $genomeName/*FinalAnn.fa TIR-Learner-Result/
rm -r $genomeName"_combine"
rm $genomeFile-+-db.nin $genomeFile-+-db.nsd $genomeFile-+-db.nsq $genomeFile-+-db.nhr $genomeFile-+-db.nog $genomeFile-+-db.nsi #shujun

echo "############################################################ TIR-Learner is finished! #################################################"


else


#echo "No Blast Result, Jump to Module3"

[ -d "Module3_New" ] || mkdir Module3_New
cd $dir"/Module3_New/"
[ -d $genomeName ] || mkdir $genomeName
[ -d "temp" ] || mkdir temp

echo "############################################################ Module 3 Begin ###########################################################"

echo "Module 3, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New" -grfp $grfp -l $len	
cp -r $genomeName/$genomeName* temp/

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module3_New/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"
cp -r $genomeName/*-p temp/

echo "Module 3, Step 3: Get dataset"
python3 $path/Module3_New/getDataset.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"


echo "Module 3, Step 5: Check TIR/TSD"
python3 $path/Module3_New/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"

echo "Module 3, Step 6: Write to Gff"
python3 $path/Module3_New/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"

echo "Module 3, Step 7: Remove Low complexity"
python3 $path/Module3_New/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"

echo "############################################################ Module 3 Finished ########################################################"

echo "Get Final GFF"
python3 $path/Module3_New/CombineAll.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"

mv *.gff3 $genomeName
rm *Low

echo "Get fasta file"
python3  $path/Module3_New/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"

mkdir TIR-Learner-Result
mv $genomeName/*FinalAnn.gff3 TIR-Learner-Result/
mv $genomeName/*FinalAnn.fa TIR-Learner-Result/
rm -r $genomeName"_combine"
rm $genomeFile-+-db.nin $genomeFile-+-db.nsd $genomeFile-+-db.nsq $genomeFile-+-db.nhr $genomeFile-+-db.nog $genomeFile-+-db.nsi #shujun

cp -r TIR-Learner-Result $dir


echo "############################################################ TIR-Learner is finished! #################################################"

fi

