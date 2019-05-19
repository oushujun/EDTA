

genome="Rice_MSU7.fasta"

if [ ! 0 ]; then
###########################
###### LTR_retriever ######
###########################

# run LTRharvest and LTR_FINDER
#~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
#~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.TIR.fa.stg0 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $name.harvest.scn

# run LTR_retriever

perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.LTRlib.fa > $genome.LTR.fa

# copy raw LTRlib to a folder for EDTA processing
mkdir $genome.LTR.EDTA_process
cp $genome.LTR.fa $genome.LTR.EDTA_process
cd $genome.LTR.EDTA_process

# identify mite contaminants with MITE-Hunter
#mkdir $genome.LTRlib.fa_MITE
#cp $genome.LTRlib.fa $genome.LTRlib.fa_MITE
#cd $genome.LTRlib.fa_MITE
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome.LTR.fa
cat $genome.LTR.fa.mites_Step8_* > $genome.LTR.fa.mite
#cp $genome.LTR.fa.mite ../
#cd ..

# identify Helitron contaminants with HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome.LTR.fa
cat $genome.LTR.fa.HelitronScanner.draw.hel.fa $genome.LTR.fa.HelitronScanner.draw.rc.hel.fa > $genome.LTR.fa.helitron

# remove potential mite and helitron contaminants
cat $genome.LTR.fa.mite $genome.LTR.fa.helitron > $genome.LTR.fa.mite.helitron
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.mite.helitron -cutoff 225 $genome.LTR.fa
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.masked > $genome.LTR.fa.HQ

# copy results to 1 folder up
cp $genome.LTR.fa.HQ ../
cd ..


###########################
######  TIR-Learner  ######
###########################

# run TIR-Learner

# convert TIR-Learner names into RepeatMasker readible names, seperate MITE (<600bp) and TIR elements
perl ~/las/git_bin/TElib_benchmark/util/rename_tirlearner.pl $genome.TIR-Learner.fa > $genome.TIR-Learner.fa.renamed

# clean up tandem repeats and short seq with cleanup_tandem.pl
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR-Learner.fa.renamed > $genome.TIR-Learner.fa.renamed.stg0


###########################
######  MITE-Hunter  ######
###########################

# run MITE-Hunter
#mkdir ${genome}_MITE
#cp $genome.LTRlib.fa $genome.LTRlib.fa_MITE
#cd ${genome}_MITE
#ln -s ../$genome $genome
#perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome
#cat $genome.mites_Step8_* > $genome.MITE-Hunter.fa
#cp $genome.MITE-Hunter.fa ../
#cd ..

# convert name to RM readible
perl -i -nle 's/MITEhunter//; print $_ and next unless /^>/; my $id = (split)[0]; print "${id}#MITE/unknown"' $genome.MITE-Hunter.fa
perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.MITE-Hunter.fa > $genome.MITE-Hunter.fa.renamed

# remove MITEs existed in TIR-Learner results, clean up tandem repeats and short seq with cleanup_tandem.pl
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.TIR-Learner.fa.renamed.stg0 $genome.MITE-Hunter.fa.renamed
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.MITE-Hunter.fa.renamed.masked > $genome.MITE-Hunter.fa.stg0

# aggregate TIR-Learner and MITE-Hunter results together
cat $genome.TIR-Learner.fa.renamed.stg0 $genome.MITE-Hunter.fa.stg0 | perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl - > $genome.TIR.fa.stg0

# identify LTR contaminants with LTRharvest
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.TIR.fa.stg0 -out $genome.TIR.fa.stg0.LTR
perl -i -nle 's/#.*\[/_/; s/\]//; s/,/_/g; print $_' $genome.TIR.fa.stg0.LTR

# identify Helitron contaminants with HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome.TIR.fa.stg0
cat $genome.TIR.fa.stg0.HelitronScanner.draw.hel.fa $genome.TIR.fa.stg0.HelitronScanner.draw.rc.hel.fa > $genome.TIR.fa.stg0.helitron

# remove potential LTR and helitron contaminants
cat $genome.TIR.fa.stg0.LTR $genome.TIR.fa.stg0.helitron > $genome.TIR.fa.stg0.LTR.helitron
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.TIR.fa.stg0.LTR.helitron $genome.TIR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.HQ


#############################
###### HelitronScanner ######
#############################

# run HelitronScanner
#sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome

# filtre out low-quality Helitron candidates
#perl ~/las/git_bin/TElib_benchmark/util/format_helitronscanner_out.pl $genome
#perl -i -nle 'print $_ and next unless /^>/; my $line=(split)[0]; $line=~s/#SUB_//; print "$line#DNA/Helitron"' $genome.HelitronScanner.filtered.fa

# clean up DNA TE and LINE coding sequence, and plant protein coding sequence
perl ~/las/git_bin/TElib_benchmark/util/cleanup_proteins.pl $genome.HelitronScanner.filtered.fa
perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.HelitronScanner.filtered.fa.clean.clean > $genome.Helitron.fa.stg0

# identify LTR contaminants with LTRharvest
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.Helitron.fa.stg0 -indexname $genome.Helitron.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.Helitron.fa.stg0 -out $genome.Helitron.fa.stg0.LTR
perl -i -nle 's/#.*\[/_/; s/\]//; s/,/_/g; print $_' $genome.Helitron.fa.stg0.LTR

# identify mite contaminants with MITE-Hunter
mkdir $genome.Helitron.fa.stg0_MITE
cp $genome.Helitron.fa.stg0 $genome.Helitron.fa.stg0_MITE/
cd $genome.Helitron.fa.stg0_MITE/
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome.Helitron.fa.stg0
cat *_Step8_* > $genome.Helitron.fa.stg0.mite
cp $genome.Helitron.fa.stg0.mite ../
cd ..

# remove potential LTR and helitron contaminants
cat $genome.Helitron.fa.stg0.LTR $genome.Helitron.fa.stg0.mite > $genome.Helitron.fa.stg0.LTR.mite
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.Helitron.fa.stg0.LTR.mite $genome.Helitron.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.masked > $genome.Helitron.fa.stg0.HQ
fi

#############################
###### rm contaminants ######
#############################

# Filtered candidates
#$genome.LTRlib.fa
#$genome.TIR.fa.stg0
#$genome.Helitron.fa.stg0

# High quality candidates
#$genome.LTRlib.fa.HQ
#$genome.TIR.fa.HQ
#$genome.Helitron.fa.stg0.HQ

# remove mite and helitron in LTR candidates
cat $genome.TIR.fa.HQ $genome.Helitron.fa.stg0.HQ > $genome.TIR.Helitron.fa.HQ
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.TIR.Helitron.fa.HQ $genome.LTRlib.fa
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTRlib.fa.masked > $genome.LTRlib.fa.stg1

# remove LTR and helitron in TIR candidates
cat $genome.LTRlib.fa.HQ $genome.Helitron.fa.stg0.HQ > $genome.LTR.Helitron.fa.HQ
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.LTR.Helitron.fa.HQ $genome.TIR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg1

# remove LTR and TIR in Helitron candidates
cat $genome.LTRlib.fa.HQ $genome.TIR.fa.HQ > $genome.LTR.TIR.fa.HQ
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.fa.HQ $genome.Helitron.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.masked > $genome.Helitron.fa.stg1

# aggregate clean sublibraries and cluster
cat $genome.LTRlib.fa.stg1 $genome.TIR.fa.stg1 $genome.Helitron.fa.stg1 > $genome.LTR.TIR.Helitron.fa
exit
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in $genome.LTR.TIR.Helitron.fa -minlen 80 -threads 30 > $genome.LTR.TIR.Helitron.fa.cln
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in $genome.LTR.TIR.Helitron.fa.cln -minlen 80 -threads 30 > $genome.LTR.TIR.Helitron.fa.cln2
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in $genome.LTR.TIR.Helitron.fa.cln2 -minlen 80 -threads 30 > $genome.LTR.TIR.Helitron.fa.cln3
mv $genome.LTR.TIR.Helitron.fa.cln3 $genome.LTR.TIR.Helitron.lib.fa
