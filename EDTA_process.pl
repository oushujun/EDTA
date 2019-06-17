#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

#######################################################################
##### Perform EDTA basic and advcanced filterings on TE candidates ####
##### Shujun Ou (shujun.ou.1@gmail.com, 05/21/2019)                ####
#######################################################################

## Input:
#	$genome.LTR.raw.fa
#	$genome.TIR.raw.fa
#	$genome.MITE.raw.fa
#	$genome.Helitron.raw.fa

## Output:
#	$genome.LTR.TIR.Helitron.fa.stg1

my $usage = "\nPerform EDTA basic and advcanced filterings for raw TE candidates and generate the stage 1 library
	perl EDTA_process.pl [options]
		-genome	[File]	The genome FASTA
		-ltr	[File]	The raw LTR library FASTA
		-tir	[File]	The raw TIR library FASTA
		-mite	[File]	The raw MITE library FASTA
		-helitron	[File]	The raw Helitron library FASTA
		-repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
		-blast [path]	The directory containing Blastn (default: read from ENV)
		-protlib [File] Protein-coding aa sequences to be removed from TE candidates. (default lib: alluniRefprexp082813 (plant))
					You may use uniprot_sprot database available from here: 
					ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
		-threads	[int]	Number of theads to run this script
		-help|-h	Display this help info
\n";

# user input
my $genome = '';
my $LTRraw = '';
my $TIRraw = '';
my $MITEraw = '';
my $Helitronraw = '';

# pre-defined
my $threads = 4;
my $script_path = $FindBin::Bin;
my $rename_TE = "$script_path/util/rename_TE.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $MITE_Hunter = "$script_path/bin/MITE-Hunter2/MITE_Hunter_manager.pl";
my $HelitronScanner = "$script_path/util/run_helitron_scanner.sh";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $rename_tirlearner = "$script_path/util/rename_tirlearner.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $protlib = "$script_path/database/alluniRefprexp082813";
my $genometools = "$script_path/bin/genometools-1.5.10/bin/gt";
my $repeatmasker = " ";
my $blast = " ";

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$LTRraw = $ARGV[$k+1] if /^-ltr$/i and $ARGV[$k+1] !~ /^-/;
	$TIRraw = $ARGV[$k+1] if /^-tir/i and $ARGV[$k+1] !~ /^-/;
	$MITEraw = $ARGV[$k+1] if /^-mite/i and $ARGV[$k+1] !~ /^-/;
	$Helitronraw = $ARGV[$k+1] if /^-helitron/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blast = $ARGV[$k+1] if /^-blast/i and $ARGV[$k+1] !~ /^-/;
	$protlib = $ARGV[$k+1] if /^-protlib/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
        }

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "LTR raw library file $LTRraw not exists!\n$usage" unless -e $LTRraw;
die "TIR raw library file $TIRraw not exists!\n$usage" unless -e $TIRraw;
die "MITE raw library file $MITEraw not exists!\n$usage" unless -e $MITEraw;
die "Helitron raw library file $Helitronraw not exists!\n$usage" unless -e $Helitronraw;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The MITE_Hunter is not found in $MITE_Hunter!\n" unless -s $MITE_Hunter;
die "The HelitronScanner is not found in $HelitronScanner!\n" unless -s $HelitronScanner;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The protein-coding sequence library is not found in $protlib!\n" unless -s $protlib;
die "The GenomeTools is not found in $genometools!\n" unless -s $genometools;

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

# Make working directories
`mkdir $genome.EDTA.LTR` unless -e "$genome.EDTA.LTR" && -d "$genome.EDTA.LTR";
`mkdir $genome.EDTA.TIR` unless -e "$genome.EDTA.TIR" && -d "$genome.EDTA.TIR";
`mkdir $genome.EDTA.Helitron` unless -e "$genome.EDTA.Helitron" && -d "$genome.EDTA.Helitron";
`mkdir $genome.EDTA.combine` unless -e "$genome.EDTA.combine" && -d "$genome.EDTA.combine";


###########################
######  Process LTR  ######
###########################

# enter the LTR folder for EDTA processing
chdir "$genome.EDTA.LTR";
`ln -s ../$LTRraw $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $rename_TE $genome.LTR.raw.fa > $genome.LTR.raw.fa.renamed`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTR.raw.fa.renamed > $genome.LTR.fa.stg0`;

# identify mite contaminants with MITE-Hunter
`rm genome* 2>/dev/null`;
`perl $MITE_Hunter -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c $threads -i $genome.LTR.fa.stg0 2>/dev/null`;
`cat *_Step8_* > $genome.LTR.fa.stg0.mite`;

# identify Helitron contaminants with HelitronScanner
`sh $HelitronScanner $genome.LTR.fa.stg0 $threads`;
`cat $genome.LTR.fa.stg0.HelitronScanner.draw.hel.fa $genome.LTR.fa.stg0.HelitronScanner.draw.rc.hel.fa | perl -nle \'print \$_ and next unless /^>/; my \$line=(split)[0]; \$line=~s/\#SUB_//; \$line=~s/[\#\\/]/_/g; print \"\$line\#DNA\/Helitron\"\' > $genome.LTR.fa.stg0.helitron`;
`cat $script_path/database/HelitronScanner.training.set.fa >> $genome.LTR.fa.stg0.helitron`;

# remove potential mite and helitron contaminants
`cat $genome.LTR.fa.stg0.mite $genome.LTR.fa.stg0.helitron > $genome.LTR.fa.stg0.mite.helitron`;
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.stg0.mite.helitron $genome.LTR.fa.stg0 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.stg0.masked > $genome.LTR.fa.stg0.cln`;

# extract LTR regions from stg0.cln as HQ
`grep \"_LTR\" $genome.LTR.fa.stg0.cln > $genome.LTR.fa.stg0.cln.list`;
`perl $output_by_list 1 $genome.LTR.fa.stg0.cln 1 $genome.LTR.fa.stg0.cln.list -FA > $genome.LTR.fa.stg0.HQ`;

# copy results to the combine folder
`cp $genome.LTR.fa.stg0 $genome.LTR.fa.stg0.HQ ../$genome.EDTA.combine`;
chdir '..';


###########################
######  Process TIR  ######
###########################

# enter the TIR folder for EDTA processing
chdir "$genome.EDTA.TIR";
`ln -s ../$TIRraw $genome.TIR.raw.fa` unless -s "$genome.TIR.raw.fa";

# convert TIR-Learner names into RepeatMasker readible names, seperate MITE (<600bp) and TIR elements
`perl $rename_tirlearner $genome.TIR.raw.fa | perl $rename_TE - > $genome.TIR.raw.fa.renamed`;

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.raw.fa.renamed > $genome.TIR_1.fa.stg0`;


###########################
######  Process MITE  #####
###########################

# Link raw MITEs to the EDTA processing folder
`ln -s ../$MITEraw $genome.MITE.raw.fa` unless -s "$genome.MITE.raw.fa";

# convert name to RM readible
`perl -i -nle \'s/MITEhunter//; print \$_ and next unless /^>/; my \$id = (split)[0]; print \"\${id}#MITE/unknown\"\' $genome.MITE.raw.fa`;
`perl $rename_TE $genome.MITE.raw.fa > $genome.MITE.raw.fa.renamed`;

# remove MITEs existed in TIR-Learner results, clean up tandem repeats and short seq with cleanup_tandem.pl
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR_1.fa.stg0 $genome.MITE.raw.fa.renamed 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.MITE.raw.fa.renamed.masked > $genome.MITE.fa.stg0`;

# aggregate TIR-Learner and MITE-Hunter results together
`cat $genome.TIR_1.fa.stg0 $genome.MITE.fa.stg0 | perl $rename_TE - > $genome.TIR.fa.stg0`;

# identify LTR contaminants with LTRharvest
`$genometools suffixerator -db $genome.TIR.fa.stg0 -indexname $genome.TIR.fa.stg0 -tis -suf -lcp -des -ssp -sds -dna`;
`$genometools ltrharvest -index $genome.TIR.fa.stg0 -out $genome.TIR.fa.stg0.LTR`;
`rm $genome.TIR.fa.stg0.esq $genome.TIR.fa.stg0.lcp $genome.TIR.fa.stg0.llv $genome.TIR.fa.stg0.md5 $genome.TIR.fa.stg0.prj $genome.TIR.fa.stg0.sds $genome.TIR.fa.stg0.suf`;
`perl -i -nle \'s/#.*\\[/_/; s/\\]//; s/,/_/g; print \$_\' $genome.TIR.fa.stg0.LTR`;

# identify Helitron contaminants with HelitronScanner
`sh $HelitronScanner $genome.TIR.fa.stg0 $threads`;
`cat $genome.TIR.fa.stg0.HelitronScanner.draw.hel.fa $genome.TIR.fa.stg0.HelitronScanner.draw.rc.hel.fa | perl -nle \'print \$_ and next unless /^>/; my \$line=(split)[0]; \$line=~s/\#SUB_//; \$line=~s/[\#\\/]/_/g; print \"\$line\#DNA\/Helitron\"\' > $genome.TIR.fa.stg0.helitron`;
`cat $script_path/database/HelitronScanner.training.set.fa >> $genome.TIR.fa.stg0.helitron`;

# remove potential LTR and helitron contaminants
`cat $genome.TIR.fa.stg0.LTR $genome.TIR.fa.stg0.helitron > $genome.TIR.fa.stg0.LTR.helitron`;
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR.fa.stg0.LTR.helitron $genome.TIR.fa.stg0 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg0.HQ`;

# copy results to the combine folder
`cp $genome.TIR.fa.stg0 $genome.TIR.fa.stg0.HQ ../$genome.EDTA.combine`;
chdir '..';


##############################
###### Process Helitron ######
##############################

# enter the Helitron folder for EDTA processing
chdir "$genome.EDTA.Helitron";
`ln -s ../$Helitronraw $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";

# format raw candidates
`perl -nle \'print \$_ and next unless /^>/; my \$line=(split)[0]; \$line=~s/\#SUB_//; print \"\$line\#DNA\/Helitron\"\' $genome.Helitron.raw.fa | perl $rename_TE - > $genome.Helitron.raw.fa.renamed`;

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.Helitron.raw.fa.renamed > $genome.Helitron.fa.stg0`;

# copy results to the combine folder
`cp $genome.Helitron.fa.stg0 ../$genome.EDTA.combine`;
chdir '..';


#################################
###### Advanced filterings ######
#################################

# enter the combine folder for EDTA processing
chdir "$genome.EDTA.combine";

# copy LTR.stg0 as LTR.stg1
`cp $genome.LTR.fa.stg0 $genome.LTR.fa.stg1`;

# remove LTR in TIR candidates
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.stg0.HQ $genome.TIR.fa.stg0 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg1`;

# remove LTR and TIR in Helitron candidates
`cat $genome.LTR.fa.stg0.HQ $genome.TIR.fa.stg0.HQ > $genome.LTR.TIR.fa.stg0.HQ`;
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.fa.stg0.HQ $genome.Helitron.fa.stg0 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.masked > $genome.Helitron.fa.stg1`;

# aggregate clean sublibraries and cluster
`cat $genome.LTR.fa.stg1 $genome.TIR.fa.stg1 $genome.Helitron.fa.stg1 > $genome.LTR.TIR.Helitron.fa.stg1.raw`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw -threads $threads -minlen 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln2`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw.cln2 -threads $threads -minlen 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln3`;

# remove protein-coding sequences
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.fa.stg1.raw.cln3 -rmdnate 0 -rmline 1 -rmprot 1 -protlib $protlib -blast $blast -threads $threads`;
`mv $genome.LTR.TIR.Helitron.fa.stg1.raw.cln3.clean $genome.LTR.TIR.Helitron.fa.stg1`;

chdir '..';

