

############################################################
##### Perform standard EDTA filterings on TE candidates ####
##### Shujun Ou (shujun.ou.1@gmail.com, 05/21/2019)     ####
############################################################

## Input:
#	$genome.LTR.raw.fa
#	$genome.TIR.raw.fa
#	$genome.MITE.raw.fa
#	$genome.Helitron.raw.fa

## Output:
#	$genome.LTR.fa.stg0, $genome.LTR.fa.stg0.HQ
#	$genome.TIR.fa.stg0, $genome.TIR.fa.stg0.HQ
#	$genome.MITE.fa.stg0, $genome.MITE.fa.stg0.HQ
#	$genome.Helitron.fa.stg0, $genome.Helitron.fa.stg0.HQ

# predefined
#!/bin/bash -login
genome="Rice_MSU7.fasta"
threads=36

# user input
genome=$1
LTRraw=$2
TIRraw=$3
MITEraw=$4
Helitronraw=$5

### allow user to specify CPU number
if [ ! -z "$6" ];
        then threads=$6
fi

# Make links to raw results
ln -s $LTRraw $genome.LTR.raw.fa
ln -s $TIRraw $genome.TIR.raw.fa
ln -s $MITEraw $genome.MITE.raw.fa
ln -s $Helitronraw $genome.Helitron.raw.fa

if [ ! 0 ]; then

###########################
######  Process LTR  ######
###########################

# copy raw LTR to a folder for EDTA processing
mkdir $genome.LTR.EDTA_process
cp $genome.LTR.raw.fa $genome.LTR.EDTA_process
cd $genome.LTR.EDTA_process

# clean up tandem repeats and short seq with cleanup_tandem.pl
perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.LTR.raw.fa > $genome.LTR.raw.fa.renamed
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTR.raw.fa.renamed > $genome.LTR.fa.stg0

# identify mite contaminants with MITE-Hunter
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c $threads -i $genome.LTR.fa.stg0
cat *_Step8_* > $genome.LTR.fa.stg0.mite

# identify Helitron contaminants with HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome.LTR.fa.stg0 $threads
cat $genome.LTR.fa.stg0.HelitronScanner.draw.hel.fa $genome.LTR.fa.stg0.HelitronScanner.draw.rc.hel.fa ~/las/git_bin/TElib_benchmark/database/HelitronScanner.training.set.fa > $genome.LTR.fa.stg0.helitron

# remove potential mite and helitron contaminants
cat $genome.LTR.fa.stg0.mite $genome.LTR.fa.stg0.helitron > $genome.LTR.fa.stg0.mite.helitron
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.stg0.mite.helitron $genome.LTR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.stg0.masked > $genome.LTR.fa.stg0.cln

# extract LTR regions from stg0.cln as HQ
grep "_LTR" $genome.LTR.fa.stg0.cln > $genome.LTR.fa.stg0.cln.list
perl ~/las/git_bin/TElib_benchmark/util/output_by_list.pl 1 $genome.LTR.fa.stg0.cln 1 $genome.LTR.fa.stg0.cln.list -FA > $genome.LTR.fa.stg0.HQ

# return to the root folder
cp $genome.LTR.fa.stg0 $genome.LTR.fa.stg0.HQ ../
cd ..


###########################
######  Process TIR  ######
###########################

# make a TIR folder for EDTA processing
mkdir $genome.TIR.EDTA_process

# convert TIR-Learner names into RepeatMasker readible names, seperate MITE (<600bp) and TIR elements
perl ~/las/git_bin/TElib_benchmark/util/rename_tirlearner.pl $genome.TIR.raw.fa | perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl - > $genome.TIR.EDTA_process/$genome.TIR.raw.fa.renamed

# clean up tandem repeats and short seq with cleanup_tandem.pl
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.EDTA_process/$genome.TIR.raw.fa.renamed > $genome.TIR.EDTA_process/$genome.TIR_1.fa.stg0


###########################
######  Process MITE ######
###########################

# Enter the EDTA processing folder
cp $genome.MITE.raw.fa $genome.TIR.EDTA_process
cd $genome.TIR.EDTA_process

# convert name to RM readible
perl -i -nle 's/MITEhunter//; print $_ and next unless /^>/; my $id = (split)[0]; print "${id}#MITE/unknown"' $genome.MITE.raw.fa
perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.MITE.raw.fa > $genome.MITE.raw.fa.renamed

# remove MITEs existed in TIR-Learner results, clean up tandem repeats and short seq with cleanup_tandem.pl
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR_1.fa.stg0 $genome.MITE.raw.fa.renamed
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.MITE.raw.fa.renamed.masked > $genome.MITE.fa.stg0

# aggregate TIR-Learner and MITE-Hunter results together
cat $genome.TIR_1.fa.stg0 $genome.MITE.fa.stg0 | perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl - > $genome.TIR.fa.stg0

# identify LTR contaminants with LTRharvest
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.TIR.fa.stg0 -out $genome.TIR.fa.stg0.LTR
perl -i -nle 's/#.*\[/_/; s/\]//; s/,/_/g; print $_' $genome.TIR.fa.stg0.LTR

# identify Helitron contaminants with HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome.TIR.fa.stg0 $threads
cat $genome.TIR.fa.stg0.HelitronScanner.draw.hel.fa $genome.TIR.fa.stg0.HelitronScanner.draw.rc.hel.fa ~/las/git_bin/TElib_benchmark/database/HelitronScanner.training.set.fa > $genome.TIR.fa.stg0.helitron

# remove potential LTR and helitron contaminants
cat $genome.TIR.fa.stg0.LTR $genome.TIR.fa.stg0.helitron > $genome.TIR.fa.stg0.LTR.helitron
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR.fa.stg0.LTR.helitron $genome.TIR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg0.HQ

# return to the root folder
cp $genome.TIR.fa.stg0 $genome.TIR.fa.stg0.HQ ../
cd ..


##############################
###### Process Helitron ######
##############################

# make a Helitron folder for EDTA processing
mkdir $genome.Helitron.EDTA_process
cp $genome.Helitron.raw.fa $genome.Helitron.EDTA_process
cd $genome.Helitron.EDTA_process

# format raw candidates
perl -nle 'print $_ and next unless /^>/; my $line=(split)[0]; $line=~s/#SUB_//; print "$line#DNA/Helitron"' $genome.Helitron.raw.fa > $genome.Helitron.raw.fa.renamed

# clean up DNA TE and LINE coding sequence, and plant protein coding sequence
#perl ~/las/git_bin/TElib_benchmark/util/cleanup_proteins.pl $genome.Helitron.raw.fa.renamed
#perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.Helitron.raw.fa.renamed.clean.clean > $genome.Helitron.fa.stg0

# remove long non-repetitive sequences in candidates
mkdir 
Red -gnm ./genome/ -msk ./output &
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nonTE.pl

# remove single copy candidates
perl ~/las/git_bin/TElib_benchmark/util/filter_copy_number.pl $genome.Helitron.raw.fa.renamed > $genome.Helitron.raw.fa.renamed.2cp
fi

# clean up tandem repeats and short seq with cleanup_tandem.pl
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.Helitron.raw.fa.renamed.2cp > $genome.Helitron.fa.stg0

# put rm long here

# identify LTR contaminants with LTRharvest
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.Helitron.fa.stg0 -indexname $genome.Helitron.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.Helitron.fa.stg0 -out $genome.Helitron.fa.stg0.LTR
perl -i -nle 's/#.*\[/_/; s/\]//; s/,/_/g; print $_' $genome.Helitron.fa.stg0.LTR

# identify mite contaminants with MITE-Hunter
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome.Helitron.fa.stg0
cat *_Step8_* > $genome.Helitron.fa.stg0.mite

# remove potential LTR and TIR contaminants
cat $genome.Helitron.fa.stg0.LTR $genome.Helitron.fa.stg0.mite > $genome.Helitron.fa.stg0.LTR.mite
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.Helitron.fa.stg0.LTR.mite $genome.Helitron.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.masked > $genome.Helitron.fa.stg0.HQ

# return to the root folder
cp $genome.Helitron.fa.stg0 $genome.Helitron.fa.stg0.HQ ../
cd ..

