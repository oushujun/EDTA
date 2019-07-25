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
#my $mindiff = 1; #minimum richness difference between $TE1 and $TE2 for a sequence to be considered as real to $TE1
my $threads = 4;
my $script_path = $FindBin::Bin;
my $TE_purifier = "$script_path/util/TE_purifier.pl";
my $rename_TE = "$script_path/util/rename_TE.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
#my $MITE_Hunter = "$script_path/bin/MITE-Hunter2/MITE_Hunter_manager.pl";
my $HelitronScanner = "$script_path/util/run_helitron_scanner.sh";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $rename_tirlearner = "$script_path/util/rename_tirlearner.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $cleanup_gff = "$script_path/util/filter_gff_TIR.pl";
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
#	$mindiff = $ARGV[$k+1] if /^-mindiff/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blast = $ARGV[$k+1] if /^-blast/i and $ARGV[$k+1] !~ /^-/;
	$protlib = $ARGV[$k+1] if /^-protlib/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
        }

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "LTR raw library file $LTRraw not exists!\n$usage" unless -e $LTRraw;
die "TIR raw library file $TIRraw not exists!\n$usage" unless -e $TIRraw;
die "MITE raw library file $MITEraw not exists!\n$usage" unless -e $MITEraw;
die "Helitron raw library file $Helitronraw not exists!\n$usage" unless -e $Helitronraw;
die "The script TE_purifier.pl is not found in $TE_purifier!\n" unless -s $TE_purifier;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;
#die "The MITE_Hunter is not found in $MITE_Hunter!\n" unless -s $MITE_Hunter;
die "The HelitronScanner is not found in $HelitronScanner!\n" unless -s $HelitronScanner;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script filter_gff_TIR.pl is not found in $cleanup_gff!\n" unless -s $cleanup_gff;
die "The protein-coding sequence library is not found in $protlib!\n" unless -s $protlib;
die "The GenomeTools is not found in $genometools!\n" unless -s $genometools;

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

# Make working directories
`mkdir $genome.EDTA.combine` unless -e "$genome.EDTA.combine" && -d "$genome.EDTA.combine";

# enter the combine folder for EDTA processing
chdir "$genome.EDTA.combine";
`ln -s ../$LTRraw $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";
`ln -s ../$TIRraw $genome.TIR.raw.fa` unless -s "$genome.TIR.raw.fa";
`ln -s ../$MITEraw $genome.MITE.raw.fa` unless -s "$genome.MITE.raw.fa";
`ln -s ../$Helitronraw $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";


if (0){

##### Make stg0 files

###########################
######  Process LTR  ######
###########################

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $rename_TE $genome.LTR.raw.fa > $genome.LTR.raw.fa.renamed`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTR.raw.fa.renamed > $genome.LTR.fa.stg0`;


###################################
######  Process TIR and MITE ######
###################################

# convert names into RepeatMasker readible names, seperate MITE (<600bp) and TIR elements
`perl $rename_tirlearner $genome.TIR.raw.fa | perl $rename_TE - > $genome.TIR.raw.fa.renamed`;
`perl -nle \'s/MITEhunter//; print \$_ and next unless /^>/; my \$id = (split)[0]; print \"\${id}#MITE/unknown\"\' $genome.MITE.raw.fa | perl $rename_TE - > $genome.MITE.raw.fa.renamed`;

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.raw.fa.renamed > $genome.TIR_1.fa.stg0`;

# remove MITEs existed in TIR-Learner results, clean up tandem repeats and short seq with cleanup_tandem.pl
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR_1.fa.stg0 $genome.MITE.raw.fa.renamed 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.MITE.raw.fa.renamed.masked > $genome.MITE.fa.stg0`;

# aggregate TIR-Learner and MITE-Hunter results together
`cat $genome.TIR_1.fa.stg0 $genome.MITE.fa.stg0 | perl $rename_tirlearner - > $genome.TIR.fa.stg0`;


##############################
###### Process Helitron ######
##############################

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl -nle \'print \$_ and next unless /^>/; my \$line=(split)[0]; \$line=~s/\#SUB_//; print \"\$line\#DNA\/Helitron\"\' $genome.Helitron.raw.fa | perl $rename_TE - > $genome.Helitron.raw.fa.renamed`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.Helitron.raw.fa.renamed > $genome.Helitron.fa.stg0`;

}


#################################
###### Advanced filterings ######
#################################

# Purify TE1 with difference in target richness between TE1 and TE2
sub Purifier() {
	my ($TE1, $TE2, $mindiff) = ($_[0], $_[1], $_[2]);
	`perl $TE_purifier -TE1 $TE1 -TE2 $TE2 -t $threads -mindiff $mindiff`;
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 1 -minrm 1 -trf 0 -f $TE1-$TE2.fa > $TE1.HQ`;
	}

# purify LTR
&Purifier("$genome.LTR.fa.stg0", "$genome.TIR.fa.stg0", 5);
&Purifier("$genome.LTR.fa.stg0.HQ", "$genome.Helitron.fa.stg0", 5);
`mv $genome.LTR.fa.stg0.HQ.HQ $genome.LTR.fa.HQ1`;
# $genome.LTR.fa.stg0 =  $genome.LTR.fa.stg1

# purify TIR
&Purifier("$genome.TIR.fa.stg0", "$genome.LTR.fa.HQ1", 1);
&Purifier("$genome.TIR.fa.stg0.HQ", "$genome.Helitron.fa.stg0", 1);
`mv $genome.TIR.fa.stg0.HQ.HQ $genome.TIR.fa.HQ1`;
# $genome.TIR.fa.HQ1 = $genome.TIR.fa.stg1

# purify Helitron
&Purifier("$genome.Helitron.fa.stg0", "$genome.LTR.fa.HQ1", 1);
&Purifier("$genome.Helitron.fa.stg0.HQ", "$genome.TIR.fa.HQ1", 1);
`mv $genome.Helitron.fa.stg0.HQ.HQ $genome.Helitron.fa.HQ1`;
# $genome.Helitron.fa.HQ1 = $genome.Helitron.fa.stg1

#clean up database files
`for i in $genome.LTR.fa.stg0 $genome.TIR.fa.stg0 $genome.Helitron.fa.stg0; do rm \$i.nhr \$i.nin \$i.nsq; done 2>/dev/null`;
#exit;

# aggregate clean sublibraries and cluster
`cat $genome.LTR.fa.stg0 $genome.TIR.fa.HQ1 $genome.Helitron.fa.HQ1 | perl -nle 's/>/\\n>/g unless /^>/; print \$_' > $genome.LTR.TIR.Helitron.fa.stg1.raw`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw -threads $threads -minlen 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln2`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw.cln2 -threads $threads -minlen 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln3`;

# remove protein-coding sequences
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.fa.stg1.raw.cln3 -rmdnate 0 -rmline 1 -rmprot 1 -protlib $protlib -blast $blast -threads $threads`;
`perl $rename_TE $genome.LTR.TIR.Helitron.fa.stg1.raw.cln3.clean > $genome.LTR.TIR.Helitron.fa.stg1`;

chdir '..';
