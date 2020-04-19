#/bin/bash
# Contributors: Weijia Su (weijia@iastate.edu); Shujun Ou (shujun.ou.1@gmail.com)
# 07/28/2019


####################
### Installation ###
####################

# conda create -n TIR-Learner
# conda activate TIR-Learner
# conda install -y -c anaconda scikit-learn=0.19.0
# conda install -y -c anaconda biopython pandas glob2 python=3.6
# conda install -y -c bioconda blast=2.5.0
# conda install -y -c conda-forge multiprocess regex tensorflow=1.14.0 keras=2.2.4

#############
### Usage ###
#############

## pre-set parameters
version="2.5"
path=$(dirname "$0") #program path
dir=$(pwd) #current work directory
rawFile="" #the genome file with relative path
grfp=$(dirname `which grf-main`)
#grfp="$path/../GenericRepeatFinder/bin/"
#grfp="~/las/bin/GenericRepeatFinder/bin/"
species="others" # One of the following "Maize", "Rice" or "others"
len=5000 # Maximum of element length used in GRF
t=4 # CPUs

## read parameters and help doc
helpFunction()
{
   echo ""
   echo "TIR-Learner Version $version"
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

# check genome
if [ ! -s "$rawFile" ]; then
   echo "ERROR: The genome file is not found or is empty"
   helpFunction
fi

# check GRF
if [ ! -s "${grfp}/grf-main" ]; then
   echo "The GRF program is not found or executable in the folder: $grfp"
   helpFunction
fi

# read file name and path
resolve_link() { #ref: https://github.com/basherpm/basher/issues/49#issuecomment-459985008
  if type -p realpath >/dev/null; then
    realpath "$1"
  elif type -p greadlink >/dev/null; then
    greadlink -f "$1"
  else
    readlink -f "$1"
  fi
}
genomeFile=`resolve_link $rawFile` #the genome file with real path
genomeName="TIR-Learner"


#####################
#### Preparations ###
#####################

echo "############################################################ Pre-Processing ###########################################################"

## Check sequence file
python3 $path/pre.py -g $genomeFile -name $genomeName

### mkdir temp
if [ $species = "Maize" ] || [ $species = "Rice" ]; then

for i in Module1 Module2 Module3; do [ -d $i ] || mkdir $i; done


################
### Module 1 ###
################

echo "############################################################ Module 1 Begin ###########################################################"

echo "Module 1, Step 1: Blast Genome against Reference Library"
cd $dir"/Module1"
[ -d $genomeName ] || mkdir $genomeName
[ -d temp ] || mkdir temp
python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1" -s $species
for i in $genomeName/*blast*; do cp $i temp/; done

echo "Module 1, Step 2: Select 100% coverage entries from Blast results"
cd $dir"/Module1"
python3 $path/Module1/Fullcov.py  -name $genomeName -p $path -d $dir"/Module1" -s $species

# test if Module 1 Step 2 produces blast results
for i in $dir/Module1/$genomeName/*-select.csv; do cat $i; done > test.select.csv
if [ ! -s test.select.csv ]; then
	echo "
	ERROR: No sequence is found similar to the TIR database of $species!
	You may have specified the wrong species. Please double check or set species=others and rerun TIR-Learner
	"
	rm test.select.csv
	exit
fi
rm test.select.csv

echo "Module 1, Step 3: Making blastDB and get candidate sequences"
python3 $path/Module1/GetFastaSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module 1, Step 4: Check TIR and TSD"
python3 $path/Module1/CheckTIRTSD_M1.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module 1, Step 5: Write to Gff3"
python3 $path/Module1/WriteToGff_M1.py -name $genomeName -p $path -t $t -d $dir"/Module1"


echo "############################################################ Module 1 Finished ########################################################"


##################
##### Module 2 ###
##################

echo "############################################################ Module 2 Begin ###########################################################"

cd $dir"/Module2/"
[ -d $genomeName ] || mkdir $genomeName
[ -d temp ] || mkdir temp

echo "Module 2, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2" -grfp $grfp -l $len

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module2/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 3: GRF result blast reference sequences"
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir"/Module2" -s $species
for i in $genomeName/*80; do cp $i temp/; done

echo "Module 2, Step 4: Get sequences from 80% similarity"
python3 $path/Module2/GetFastaSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 5: Check TIR and TSD"
python3 $path/Module2/CheckTIRTSD_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2, Step 6: Write to gff3"
python3 $path/Module2/WriteToGff_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

# clean up intermediate files
rm TIR-Learner/TIR-Learner*fasta > /dev/null 2>&1
rm TIR-Learner/TIR-Learner*-+-200* > /dev/null 2>&1

echo "############################################################ Module 2 Finished ########################################################"


################
### Module 3 ###
################

echo "############################################################ Module 3 Begin ###########################################################"

cd $dir"/Module3/"
[ -d $genomeName ] || mkdir $genomeName
[ -d temp ] || mkdir temp

mv ../Module2/$genomeName/*_nonHomo.fa $genomeName/

echo "Module 3, Step 1: Get dataset and ML prediction"
python3 $path/Module3/getDataset.py -name $genomeName -p $path -t $t -d $dir"/Module3" -g $genomeFile
for i in $genomeName/*nonHomo.fa; do cp $i temp/; done

echo "Module 2, Step 3: Check TIR/TSD"
python3 $path/Module3/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "Module 3, Step 4: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"

# clean up intermediate files
rm TIR-Learner/TIR-Learner*fasta > /dev/null 2>&1
rm TIR-Learner/TIR-Learner*-+-200* > /dev/null 2>&1

echo "############################################################ Module 3 Finished ########################################################"


######################
### Postprocessing ###
######################

echo "############################################################ Post Processing  #########################################################"

cd $dir
[ -d $genomeName ] || mkdir $genomeName

echo "Get Final GFF"
python3 $path/Module3/CombineAll.py -name $genomeName -p $path -t $t -d $dir
cp *.gff3 $genomeName

echo "Get fasta file"
python3 $path/Module3/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

[ -d "TIR-Learner-Result" ] || mkdir TIR-Learner-Result
cp $genomeName/*FinalAnn*.gff3 TIR-Learner-Result/
cp $genomeName/*FinalAnn*.fa TIR-Learner-Result/

echo "############################################################ TIR-Learner is finished! #################################################"


else


####################
### Module 3 new ###
####################

echo "Jump to Module3"

[ -d "$dir/Module3_New" ] || mkdir $dir/Module3_New
cd "$dir/Module3_New/"
[ -d $genomeName ] || mkdir $genomeName
[ -d "temp" ] || mkdir temp


echo "############################################################ Module 3 Begin ###########################################################"

echo "Module 3, Step 1: Split Genome and Run GRF program to find Inverted Repeats"
#rm ./$genomeName/*
python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New" -grfp $grfp -l $len
for i in $genomeName/$genomeName*; do cp $i temp/; done

echo "Module 3, Step 2: Process GRF results"
python3 $path/Module2/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"
for i in $genomeName/*-p; do cp $i temp/; done

echo "Module 3, Step 3: Get dataset"
export OMP_NUM_THREADS=1
python3 $path/Module3_New/getDataset.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"

echo "Module 3, Step 4: Check TIR/TSD"
python3 $path/Module3_New/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"

echo "Module 3, Step 5: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"

rm temp/*fasta > /dev/null 2>&1

echo "############################################################ Module 3 Finished ########################################################"

echo "Get Final GFF"
python3 $path/Module3_New/CombineAll.py -name $genomeName -p $path -t $t -d $dir"/Module3_New"
mv *.gff3 $genomeName

echo "Get fasta file"
python3  $path/Module3/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3_New"

[ -d TIR-Learner-Result ] || mkdir TIR-Learner-Result
mv $genomeName/*FinalAnn*.gff3 TIR-Learner-Result/
mv $genomeName/*FinalAnn*.fa TIR-Learner-Result/

cp -r TIR-Learner-Result $dir

# clean up intermediate files
rm TIR-Learner/TIR-Learner*fasta > /dev/null 2>&1
rm TIR-Learner/TIR-Learner*-+-200* > /dev/null 2>&1

echo "############################################################ TIR-Learner is finished! #################################################"

fi

