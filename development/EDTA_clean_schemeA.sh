

############################################################
##### Perform advanced EDTA filterings on TE candidates ####
##### Shujun Ou (shujun.ou.1@gmail.com, 05/21/2019)     ####
############################################################


## Input (Filtered candidates and high-quality candidates)
#       $genome.LTR.fa.stg0, $genome.LTR.fa.stg0.HQ
#       $genome.TIR.fa.stg0, $genome.TIR.fa.stg0.HQ
#       $genome.Helitron.fa.stg0, $genome.Helitron.fa.stg0.HQ

## Output:
#       $genome.LTR.fa.stg1
#       $genome.TIR.fa.stg1
#       $genome.Helitron.fa.stg1


# predefined
#!/bin/bash -login
genome="Rice_MSU7.fasta"
threads=36

# user input
genome=$1
LTRstg0=$2
LTRstg0HQ=$3

TIRstg0=$4
TIRstg0HQ=$5

Helitronstg0=$6
Helitronstg0HQ=$7

### allow user to specify CPU number
if [ ! -z "$8" ];
        then threads=$8
fi

# Make links to raw results
ln -s $LTRstg0 $genome.LTR.fa.stg0
ln -s $LTRstg0HQ $genome.LTR.fa.stg0.HQ
ln -s $TIRstg0 $genome.TIR.fa.stg0
ln -s $TIRstg0HQ $genome.TIR.fa.stg0.HQ
ln -s $Helitronstg0 $genome.Helitron.fa.stg0
ln -s $Helitronstg0HQ $genome.Helitron.fa.stg0.HQ


#############################
###### rm contaminants ######
#############################

# copy raw and HQ LTR TIR Helitron files for EDTA processing
mkdir $genome.combine.EDTA_process
cp $genome.LTR.fa.stg0 $genome.TIR.fa.stg0 $genome.Helitron.fa.stg0 $genome.LTR.fa.stg0.HQ $genome.TIR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ $genome.combine.EDTA_process
cd $genome.combine.EDTA_process

# remove mite and helitron in LTR candidates
#cat $genome.TIR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ > $genome.TIR.Helitron.fa.stg0.HQ
#RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR.Helitron.fa.stg0.HQ $genome.LTR.fa.stg0
#perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.stg0.masked > $genome.LTR.fa.stg1

# remove LTR and helitron in TIR candidates
#cat $genome.LTR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ > $genome.LTR.Helitron.fa.stg0.HQ
#RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.Helitron.fa.stg0.HQ $genome.TIR.fa.stg0
#perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg1

# remove LTR in TIR candidates
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.stg0.HQ $genome.TIR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg1

# remove TIR in LTR candidates
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR.fa.stg1 $genome.LTR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.stg0.masked > $genome.LTR.fa.stg1

# remove LTR and TIR in Helitron candidates
#cat $genome.LTR.fa.stg0.HQ $genome.TIR.fa.stg0.HQ > $genome.LTR.TIR.fa.stg0.HQ
#RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.fa.stg0.HQ $genome.Helitron.fa.stg0
cat $genome.LTR.fa.stg1 $genome.TIR.fa.stg1 > $genome.LTR.TIR.fa.stg1
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.fa.stg1 $genome.Helitron.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.masked > $genome.Helitron.fa.stg1



# aggregate clean sublibraries and cluster
cat $genome.LTR.fa.stg1 $genome.TIR.fa.stg1 $genome.Helitron.fa.stg1 > $genome.LTR.TIR.Helitron.fa.stg1

