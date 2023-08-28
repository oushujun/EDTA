#!/usr/bin/env bash
# Contributors: Weijia Su (weijia@iastate.edu); Tianyu Lu (tlu83@wisc.edu); Shujun Ou (shujun.ou.1@gmail.com)
# 2023-08-09

####################
### Installation ###
####################

# conda create -n TIR-Learner
# conda activate TIR-Learner
# conda install -y -c anaconda scikit-learn
# conda install -y -c anaconda -c bioconda blast biopython pandas glob2 python
# conda install -y -c conda-forge multiprocess regex tensorflow keras

#############
### Usage ###
#############

## pre-set parameters
version="3.0"
path=$(dirname "$0") #program path
rawFile=""           #the genome file with relative path
grfp=$(dirname "$(which grf-main)")
species="Others" # One of the following "Maize", "Rice" or "Others"
len=5000         # Maximum of element length used in GRF
t=16             # Number of processor
outputDir="None" # Output directory
debug=""         # Debug mode argument statement

## read parameters and help doc
helpFunction() {
  echo ""
  echo "TIR-Learner Version $version"
  echo "Usage: sh $0 -g genome.fa -s [Rice|Maize|Others] -p [path] -l [int] -t [int] -o []"
  echo -e "\t-g	Genome file in multi-FASTA format"
  echo -e "\t-s	Specify species: Rice, Maize, or Others. Default: Others"
  echo -e "\t-p	Specify path to the GRF program"
  echo -e "\t-l	Maximum spacer length between the terminal repeat (used in GRF). Default: 5000 (bp)"
  echo -e "\t-t	Maximum processors used in this program. Default: 16"
  echo -e "\t-o	Output directory. Default: directory of the genome file"
  echo -e "\t-d Enable debug mode, which will output each module's result in csv file (Optional)"
  echo ""
  exit 1 # Exit script after printing help
}

while getopts "g:s:p:l:t:o:d" opt; do
  case ${opt} in
  g) rawFile="$OPTARG" ;;
  s) species="$OPTARG" ;;
  p) grfp="$OPTARG" ;;
  l) len="$OPTARG" ;;
  t) t="$OPTARG" ;;
  o) outputDir="$OPTARG" ;;
  d) debug="-d";;
  ?) helpFunction ;; # Print helpFunction in case parameter is non-existent
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
genomeFile=$(resolve_link $rawFile) #the genome file with real path
genomeName="TIR-Learner"

#####################
#### Preparations ###
#####################

# echo "Current Work Directory: "$dir

echo "############################################################ Pre-Processing ###########################################################"

## Check sequence file
python3 $path/pre.py -g $genomeFile -name $genomeName

#################
### Execution ###
#################

python3 $path/bin/main.py -f $genomeFile -n $genomeName -s $species -c $path -g $grfp -l $len -t $t -o $outputDir $debug

echo "############################################################ TIR-Learner is finished! #################################################"
