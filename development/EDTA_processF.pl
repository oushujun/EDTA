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
#	$genome.Helitron.raw.fa

## Output:
#	$genome.LTR.TIR.Helitron.fa.stg1

my $usage = "\nPerform EDTA basic and advcanced filterings for raw TE candidates and generate the stage 1 library
	perl EDTA_processF.pl [options]
		-genome	[File]	The genome FASTA
		-ltr	[File]	The raw LTR library FASTA
		-tir	[File]	The raw TIR library FASTA
		-helitron	[File]	The raw Helitron library FASTA
		-mindiff	[float]	The minimum fold difference in richness between main TE and contaminants
		-repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
		-blast [path]	The directory containing Blastn (default: read from ENV)
		-protlib [File] Protein-coding aa sequences to be removed from TE candidates. (default lib: alluniRefprexp082813 (plant))
					You may use uniprot_sprot database available from here:
					ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
		-threads|-t	[int]	Number of theads to run this script
		-help|-h	Display this help info
\n";

# user input
my $genome = '';
my $LTRraw = '';
my $TIRraw = '';
my $Helitronraw = '';

# pre-defined
my $mindiff = 1; #minimum richness difference between $TE1 and $TE2 for a sequence to be considered as real to $TE1
my $maxit = 1; # rounds of iterated purge
my $threads = 4;
my $script_path = $FindBin::Bin;
my $TE_purifier = "$script_path/util/TE_purifier.pl";
my $rename_TE = "$script_path/util/rename_TE.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $rename_tirlearner = "$script_path/util/rename_tirlearner.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $cleanup_gff = "$script_path/util/filter_gff_TIR.pl";
my $protlib = "$script_path/database/alluniRefprexp082813";
my $repeatmasker = " ";
my $blast = " ";

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$LTRraw = $ARGV[$k+1] if /^-ltr$/i and $ARGV[$k+1] !~ /^-/;
	$TIRraw = $ARGV[$k+1] if /^-tir/i and $ARGV[$k+1] !~ /^-/;
	$Helitronraw = $ARGV[$k+1] if /^-helitron/i and $ARGV[$k+1] !~ /^-/;
	$mindiff = $ARGV[$k+1] if /^-mindiff/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blast = $ARGV[$k+1] if /^-blast/i and $ARGV[$k+1] !~ /^-/;
	$protlib = $ARGV[$k+1] if /^-protlib/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
        }

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "LTR raw library file $LTRraw not exists!\n$usage" unless -s $LTRraw;
die "TIR raw library file $TIRraw not exists!\n$usage" unless -s $TIRraw;
die "Helitron raw library file $Helitronraw not exists!\n$usage" unless -s $Helitronraw;
die "The script TE_purifier.pl is not found in $TE_purifier!\n" unless -s $TE_purifier;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script filter_gff_TIR.pl is not found in $cleanup_gff!\n" unless -s $cleanup_gff;
die "The protein-coding sequence library is not found in $protlib!\n" unless -s $protlib;

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
`ln -s ../$Helitronraw $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";


##################################
######  define subroutines  ######
##################################

# identify TIR contaminants with TIRvish
sub Purifier() {
	my ($TE1, $TE2, $mindiff, $k) = ($_[0], $_[1], $_[2], $_[3]);
	`perl $TE_purifier -TE1 $TE1 -TE2 $TE2 -t $threads -mindiff $mindiff`;
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 1 -minrm 1 -trf 0 -f $TE1-$TE2.fa > $TE1.HQ$k`;
	}

# purge sequence from library
sub RMclean() {
	my ($lib, $file, $minlen, $trf) = ($_[0], $_[1], $_[2], $_[3]);
	`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $lib $file 2>/dev/null`;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen $minlen -minscore 3000 -trf $trf -cleanN 1 -cleanT 1 -f $file.masked > $file.cln`;
	}


##### Make stg0 and HQ0 files

###########################
######  Process LTR  ######
###########################

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $rename_TE $genome.LTR.raw.fa > $genome.LTR.raw.fa.renamed`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.LTR.raw.fa.renamed > $genome.LTR.fa.stg0`;


###########################
######  Process TIR  ######
###########################

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl $rename_tirlearner $genome.TIR.raw.fa | perl $rename_TE - > $genome.TIR.raw.fa.renamed`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.raw.fa.renamed > $genome.TIR.fa.stg0`;


##############################
###### Process Helitron ######
##############################

# clean up tandem repeats and short seq with cleanup_tandem.pl
`perl -nle \'print \$_ and next unless /^>/; my \$line=(split)[0]; \$line=~s/\#SUB_//; print \"\$line\#DNA\/Helitron\"\' $genome.Helitron.raw.fa | perl $rename_TE - > $genome.Helitron.raw.fa.renamed`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.Helitron.raw.fa.renamed > $genome.Helitron.fa.stg0`;


#################################
###### Advanced filterings ######
#################################

# predefine variables
my $LTR = "$genome.LTR.fa.stg0";
my $TIR = "$genome.TIR.fa.stg0";
my $HEL = "$genome.Helitron.fa.stg0";

# use iterations to purge contaminants
for (my $i=1; $i<=$maxit; $i++){
	# purify LTR
	&Purifier("$LTR", "$TIR", $mindiff, $i);
	&Purifier("$LTR.HQ$i", "$HEL", $mindiff, $i);
	`mv $LTR.HQ$i.HQ$i $LTR.HQ$i`;

	# purify Helitron
	&Purifier("$HEL", "$TIR", $mindiff, $i);
	&Purifier("$HEL.HQ$i", "$LTR", $mindiff, $i);
	`mv $HEL.HQ$i.HQ$i $HEL.HQ$i`;

	# purify TIR
	&Purifier("$TIR", "$LTR", $mindiff, $i);
	&Purifier("$TIR.HQ$i", "$HEL", $mindiff, $i);
	`mv $TIR.HQ$i.HQ$i $TIR.HQ$i`;

	# rename to iterate
	$LTR = "$LTR.HQ$i";
	$TIR = "$TIR.HQ$i";
	$HEL = "$HEL.HQ$i";
	}

# aggregate clean sublibraries and cluster
`cat $genome.LTR.fa.stg0 $genome.TIR.fa.stg0.*HQ$maxit $genome.Helitron.fa.stg0.*HQ$maxit | perl -nle 's/>/\\n>/g unless /^>/; print \$_' > $genome.LTR.TIR.Helitron.fa.stg1.raw`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.fa.stg1.raw -threads $threads -minlen 80 -cov 0.95 -iter 3 -blastplus $blast`;

# remove protein-coding sequences
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.fa.stg1.raw.cln -rmdnate 0 -rmline 1 -rmprot 1 -protlib $protlib -blast $blast -threads $threads`;
`perl $rename_TE $genome.LTR.TIR.Helitron.fa.stg1.raw.cln.clean > $genome.LTR.TIR.Helitron.fa.stg1`;

chdir '..';

