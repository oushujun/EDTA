#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=05:20:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=3   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing

#SBATCH --job-name='testNew'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin
module load python/3
module load genometools
## pre-set parameters
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.16" #program path
dir=$(pwd) #current work directory
rawFile="/work/LAS/thomasp-lab/weijia/research/Genomes/MaizeB73/MaizeB73.fa" #the genome file with relative path
species="Maize" # One of the following "Maize", "Rice" or "others"
len=20000 # Maximum of element length used in GRF
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
genomeName="TIR-Learner"




echo "############################################################ Module 3 Begin ###########################################################"
#gt suffixerator -db $genomeFile -indexname $genomeName -tis -suf -lcp -des -ssp -sds -dna -mirrored
#gt tirvish -index $genomeName -seed 20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist $len -similar 80 -mintsd 2 -maxtsd 11 -vic 13 -seqids 'yes' > $genomeName"_TIRvish.gff"
#cd $dir"/Module3_New/"
#cp -r $genomeName/$genomeName* temp/

#echo "Module 2, Step 2: Process GRF results"
#python3 $path/Module2/ProcessTIRvish.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"
#cp -r $genomeName/*-p temp/
python3 GetSeq.py 
#python3 CheckTIR.py
python3 CheckTIRTSD_M2.py
echo "Module 3, Step 3: Get dataset"
#python3 $path/Module3_New/getDataset.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"

#echo "Module 3, Step 4: ML prediction"
#python3 $path/Module3_New/ML_Ensemble.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New" -s $species
#
#echo "Module 3, Step 5: Check TIR/TSD"
#python3 $path/Module3_New/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"
#
#echo "Module 3, Step 6: Write to Gff"
#python3 $path/Module3_New/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"
#
#echo "Module 3, Step 7: Remove Low complexity"
#python3 $path/Module3_New/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"
#
#echo "############################################################ Module 3 Finished ########################################################"
#
#echo "Get Final GFF"
#python3 $path/Module3_New/CombineAll.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"
#
#mv *.gff3 $genomeName
#rm *Low
#
#echo "Get fasta file"
#python3  $path/Module3_New/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"
#
#mkdir TIR-Learner-Result
#mv $genomeName/*FinalAnn.gff3 TIR-Learner-Result/
#mv $genomeName/*FinalAnn.fa TIR-Learner-Result/
#rm -r $genomeName"_combine"
#rm $genomeFile-+-db.nin $genomeFile-+-db.nsd $genomeFile-+-db.nsq $genomeFile-+-db.nhr $genomeFile-+-db.nog $genomeFile-+-db.nsi #shujun
#
#cp -r TIR-Learner-Result $dir
#

echo "############################################################ TIR-Learner is finished! #################################################"



