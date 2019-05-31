

###########################
###### LTR_retriever ######
###########################

# run LTRharvest and LTR_FINDER
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna
~/las/git_bin/TElib_benchmark/bin/genometools-1.5.10/bin/gt ltrharvest -index $genome.TIR.fa.stg0 -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $name.harvest.scn
cp $genome.LTRlib.fa ../$genome.LTR.raw.fa
cd ../



###########################
######  TIR-Learner  ######
###########################

# run TIR-Learner



###########################
######  MITE-Hunter  ######
###########################

# run MITE-Hunter
mkdir ${genome}_MITE
cp $genome.LTRlib.fa $genome.LTRlib.fa_MITE
cd ${genome}_MITE
ln -s ../$genome $genome
perl ~/las/git_bin/TElib_benchmark/bin/MITE-Hunter2/MITE_Hunter_manager.pl -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c 16 -i $genome
cat *_Step8_* > $genome.MITE-Hunter.raw.fa
cp $genome.MITE-Hunter.raw.fa ../
cd ..

#############################
###### HelitronScanner ######
#############################

# run HelitronScanner
sh ~/las/git_bin/Plant_TE_annotation/helitron/run_helitron_scanner.sh $genome $threads

# filtre out low-quality Helitron candidates
perl ~/las/git_bin/TElib_benchmark/util/format_helitronscanner_out.pl $genome
perl -i -nle 'print $_ and next unless /^>/; my $line=(split)[0]; $line=~s/#SUB_//; print "$line#DNA/Helitron"' $genome.HelitronScanner.filtered.fa

# filter out low-copy Helitron candidates

cp $genome.HelitronScanner.filtered.fa $genome.Helitron.raw.fa

# clean up DNA TE and LINE coding sequence, and plant protein coding sequence
perl ~/las/git_bin/TElib_benchmark/util/cleanup_proteins.pl $genome.Helitron.raw.fa.renamed
perl ~/las/git_bin/TElib_benchmark/util/rename_TE.pl $genome.Helitron.raw.fa.renamed.clean.clean > $genome.Helitron.fa.stg0


