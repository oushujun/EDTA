#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

#######################################################################
##### Perform EDTA basic and advcanced filterings on TE candidates ####
##### Shujun Ou (shujun.ou.1@gmail.com, 11/04/2019)                ####
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
		-mindiff_ltr	[float]	The minimum fold difference in richness between LTRs and contaminants (default: 1)
		-mindiff_tir	[float]	The minimum fold difference in richness between TIRs and contaminants (default: 1)
		-mindiff_hel	[float]	The minimum fold difference in richness between Helitrons and contaminants (default: 4)
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
my $mindiff_LTR = 1; #minimum richness difference between $TE1 and $TE2 for a sequence to be considered as real to $TE1
my $mindiff_TIR = 1;
my $mindiff_HEL = 2;

my $threads = 4;
my $script_path = $FindBin::Bin;
my $TE_purifier = "$script_path/util/TE_purifier.pl";
my $rename_TE = "$script_path/util/rename_TE.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
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
	$mindiff_LTR = $ARGV[$k+1] if /^-mindiff_ltr/i and $ARGV[$k+1] !~ /^-/;
	$mindiff_TIR = $ARGV[$k+1] if /^-mindiff_tir/i and $ARGV[$k+1] !~ /^-/;
	$mindiff_HEL = $ARGV[$k+1] if /^-mindiff_hel/i and $ARGV[$k+1] !~ /^-/;
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
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The protein-coding sequence library is not found in $protlib!\n" unless -s $protlib;

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

# Make working directories
#`mkdir $genome.EDTA.combine` unless -e "$genome.EDTA.combine" && -d "$genome.EDTA.combine";
`mkdir $genome.EDTA.combineI2` unless -e "$genome.EDTA.combineI2" && -d "$genome.EDTA.combineI2";

# enter the combine folder for EDTA processing
#chdir "$genome.EDTA.combine";
chdir "$genome.EDTA.combineI2";
`ln -s ../$LTRraw $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";
`ln -s ../$TIRraw $genome.TIR.raw.fa` unless -s "$genome.TIR.raw.fa";
`ln -s ../$Helitronraw $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";


##################################
######  define subroutines  ######
##################################

sub Purifier() {
	my ($TE1, $TE2, $mindiff) = ($_[0], $_[1], $_[2]);
	`perl $TE_purifier -TE1 $TE1 -TE2 $TE2 -t $threads -mindiff $mindiff`;
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 1 -minrm 1 -trf 0 -f $TE1-$TE2.fa > $TE1.HQ`;
	}


#################################
###### Advanced filterings ######
#################################

## predefined variables
my $LTR = "$genome.LTR.raw.fa";
my $TIR = "$genome.TIR.raw.fa";
my $HEL = "$genome.Helitron.raw.fa";

## Purge contaminants
# purify LTR
&Purifier("$LTR", "$TIR", $mindiff_LTR);
&Purifier("$LTR.HQ", "$HEL", $mindiff_LTR);
`mv $LTR.HQ.HQ $LTR.HQ`;

# purify Helitron
#&Purifier("$HEL", "$TIR", $mindiff_HEL);
#&Purifier("$HEL.HQ", "$LTR", $mindiff_HEL);
&Purifier("$HEL", "$TIR", 1);
&Purifier("$HEL.HQ", "$LTR", 4);
`mv $HEL.HQ.HQ $HEL.HQ`;

# purify TIR
&Purifier("$TIR", "$LTR", $mindiff_TIR);
&Purifier("$TIR.HQ", "$HEL", $mindiff_TIR);
`mv $TIR.HQ.HQ $TIR.HQ`;

## clean LTRs in TIRs and Helitrons
`cat $TIR.HQ $HEL.HQ | perl -nle 's/>/\\n>/g unless /^>/; print \$_' > $genome.TIR.Helitron.fa.stg1.raw`;
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $LTR.HQ $genome.TIR.Helitron.fa.stg1.raw 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.Helitron.fa.stg1.raw.masked > $genome.TIR.Helitron.fa.stg1.raw.cln`;

## cluster TIRs and Helitrons and make stg1 raw library
`perl $cleanup_nested -in $genome.TIR.Helitron.fa.stg1.raw.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast`;
`cat $LTR $genome.TIR.Helitron.fa.stg1.raw.cln.cln > $genome.LTR.TIR.Helitron.fa.stg1.raw.cln`;

## remove protein-coding sequences
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.fa.stg1.raw.cln -rmdnate 0 -rmline 1 -rmprot 1 -protlib $protlib -blast $blast -threads $threads`;
`perl $rename_TE $genome.LTR.TIR.Helitron.fa.stg1.raw.cln.clean > $genome.LTR.TIR.Helitron.fa.stg1`;

chdir '..';

