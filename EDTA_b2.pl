#!/usr/bin/perl -w
use strict;

my $genome = "Rice_MSU7.fasta";

#task
#1. double check cleamup_tandem.pl para

###########################
###### LTR_retriever ######
###########################

# run LTRharvest and LTR_FINDER
#~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
#~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.TIR.fa.stg0 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $name.harvest.scn

# run LTR_retriever

# identify mite contaminants with MITE-Hunter
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome.LTRlib.fa
cat MITE-Hunter_edu/*_Step8_* > $genome.LTRlib.fa.mite

# identify Helitron contaminants with HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh Rice_MSU7.fasta.LTRlib.fa
cat $genome.LTRlib.fa.HelitronScanner.draw.hel.fa $genome.LTRlib.fa.HelitronScanner.draw.rc.hel.fa > $genome.LTRlib.fa.helitron

# remove potential mite and helitron contaminants
cat $genome.LTRlib.fa.mite $genome.LTRlib.fa.helitron > $genome.LTRlib.fa.mite.helitron
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.LTRlib.fa.mite.helitron -cutoff 225 $genome.LTRlib.fa
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTRlib.fa.masked > $genome.LTRlib.fa.HQ


###########################
######  TIR-Learner  ######
###########################

# run TIR-Learner

# convert TIR-Learner names into RepeatMasker readible names, seperate MITE (<600bp) and TIR elements
perl ~/las/git_bin/TElib_benchmark/util/rename_tirlearner.pl $genome.TIR-Learner.fa > $genome.TIR-Learner.fa.renamed

# clean up tandem repeats and short seq with cleanup_tandem.pl
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR-Learner.fa.renamed > $genome.TIR-Learner.fa.renamed.stg0



# clean up LTR seq with LTR-retriever-derived library, clean up tandem repeats and short seq with cleanup_tandem.pl
#RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib Rice_MSU7.fasta.LTRlib.fa.masked.cln -cutoff 225 TIR-Learner_Rice_0425.fa.renamed &
#perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -minlen 80 -cleanN 1 -cleanT 1 -trf 1 -f TIR-Learner_Rice_0425.fa.renamed > TIR-Learner_Rice_0425.fa.renamed.masked.cln

# make TIR and LTR library
#cat TIR-Learner_Rice_0425.fa.renamed.masked.cln $genome.LTRlib.fa.masked.cln > $genome.LTR_TIR.lib.fa



###########################
######  MITE-Hunter  ######
###########################

# run MITE-Hunter
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome

# get all candidate sequences
cat $genome.mite/$genome.mites_Step8_* > $genome.MITE-Hunter.fa

# convert name to RM readible
perl -i -nle 's/MITEhunter//; print $_ and next unless /^>/; my $id = (split)[0]; print "${id}#MITE/unknown"' $genome.MITE-Hunter.fa

# remove MITEs existed in TIR-Learner results, clean up tandem repeats and short seq with cleanup_tandem.pl
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.TIR-Learner.fa.renamed.tg0 $genome.MITE-Hunter.fa
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.MITE-Hunter.fa > $genome.MITE-Hunter.fa.stg0

# aggregate TIR-Learner and MITE-Hunter results together
cat $genome.TIR-Learner.fa.renamed.stg0 $genome.MITE-Hunter.fa.stg0 > $genome.TIR.fa.stg0

# identify LTR contaminants with LTRharvest
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.TIR.fa.stg0 -seqids yes > $genome.TIR.fa.stg0.harvest.scn
perl ~/las/git_bin/TElib_benchmark/util/get_range.pl 1 $genome.TIR.fa.stg0 $genome.TIR.fa.stg0.harvest.scn -i -g
perl ~/las/git_bin/TElib_benchmark/util/call_seq_by_list.pl $genome.TIR.fa.stg0.harvest.scn.list -C $genome.TIR.fa.stg0 > $genome.TIR.fa.stg0.LTR
perl -i -nle 's/\|.*//; print $_' $genome.TIR.fa.stg0.LTR

# identify Helitron contaminants with HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome.TIR.fa.stg0
cat $genome.TIR.fa.stg0.HelitronScanner.draw.hel.fa $genome.TIR.fa.stg0.HelitronScanner.draw.rc.hel.fa > $genome.TIR.fa.stg0.helitron

# remove potential LTR and helitron contaminants
cat $genome.TIR.fa.stg0.LTR $genome.TIR.fa.stg0.helitron > $genome.TIR.fa.stg0.LTR.helitron
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $genome.TIR.fa.stg0.LTR.helitron $genome.TIR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.HQ




# clean up MITEs with TIR-Learner results, clean up tandem repeats and short seq with cleanup_tandem.pl
nohup RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib Rice_MSU7.fasta.LTR_TIR.lib.fa -cutoff 225 MITE-Hunter2_0421_edu.fa.mod
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f MITE-Hunter2_0421_edu.fa.mod.masked > MITE-Hunter2_0421_edu.fa.mod.masked.cln

# make LTR-TIR-MITE library
cat Rice_MSU7.fasta.LTR_TIR.lib.fa MITE-Hunter2_0421_edu.fa.mod.masked.cln > Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa

# Cluster sequences and remove nested insertions
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa -minlen 80 -threads 30 > Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.cln

# run HelitronScanner to find Helitron contaminations
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.cln
perl ~/las/git_bin/TElib_benchmark/util/format_helitronscanner_out.pl Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.cln

# Add curated Helitron head and tail sequences to the Helitron mask lib
cat ~/las/git_bin/TElib_benchmark/database/HelitronScanner.training.set.fa >> Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.HelitronScanner.filtered.fa

# remove Helitron contaminations
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.HelitronScanner.filtered.fa -cutoff 225 Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.masked > Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa.masked.cln


#############################
###### HelitronScanner ######
#############################

# run HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome

# filtre out low-quality Helitron candidates
perl ~/las/git_bin/TElib_benchmark/util/format_helitronscanner_out.pl $genome
perl -nle 'print $_ and next unless /^>/; my $line=(split)[0]; $line=~s/#SUB_//; print "$line#DNA/Helitron"' $genome.HelitronScanner.filtered.fa > $genome.HelitronScanner.filtered.fa.mod

#clean up non-Helitron with LTR-TIR-MITE library, clean up tandem repeats and short seq with cleanup_tandem.pl
nohup RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa -cutoff 225 NIP_TIGR7.fasta.HelitronScanner.filtered.fa.mod
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f NIP_TIGR7.fasta.HelitronScanner.filtered.fa.mod.masked > NIP_TIGR7.fasta.HelitronScanner.filtered.fa.mod.masked.cln

#clean up DNA TE and LINE coding sequence, and plant protein coding sequence
perl ~/las/git_bin/TElib_benchmark/util/cleanup_proteins.pl NIP_TIGR7.fasta.HelitronScanner.filtered.fa.mod.masked.cln



## make structural library, reduce redundancy, and mask the genome
cat Rice_MSU7.fasta.LTR_TIR_MITE.lib.fa NIP_TIGR7.fasta.HelitronScanner.filtered.fa.mod.masked.cln.clean.clean > Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa
cd-hit-est -i Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa -o Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa.clust -c 0.8 -G 0.8 -s 0.9 -T 20 -aL 0.9 -aS 0.95 -T 30 &
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa -minlen 80 -threads 30 > Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa.cln &
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa.cln -minlen 80 -threads 30 > Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa.cln2 &
perl ~/las/git_bin/TElib_benchmark/util/cleanup_nested.pl -in Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa.cln2 -minlen 80 -threads 30 > Rice_MSU7.fasta.LTR_TIR_MITE_Helitron.lib.fa.cln3 &



