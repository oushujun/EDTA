#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;

#####################################################################
##### Perform EDTA basic and advance filtering on TE candidates #####
##### Shujun Ou (shujun.ou.1@gmail.com, 11/04/2019)             #####
#####################################################################

## Input:
#	$genome.LTR.raw.fa
#	$genome.nonLTR.raw.fa
#	$genome.TIR.raw.fa
#	$genome.Helitron.raw.fa

## Output:
#	$genome.EDTA.fa.stg1

my $usage = "\nPerform EDTA basic and advance filtering for raw TE candidates and generate the stage 1 library
	perl EDTA_processF.pl [options]
		-genome	[File]	The genome FASTA
		-ltr	[File]	The raw LTR library FASTA
		-nonltr	[File]	The raw nonLTR library FASTA
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
my $nonLTRraw = '';
my $TIRraw = '';
my $Helitronraw = '';
my $err = '';

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
	$nonLTRraw = $ARGV[$k+1] if /^-nonltr$/i and $ARGV[$k+1] !~ /^-/;
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

# define RepeatMasker -pa parameter
my $rm_threads = int($threads/4);

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "LTR raw library file $LTRraw not exists!\n$usage" unless -s $LTRraw;
die "nonLTR raw library file $nonLTRraw not exists!\n$usage" unless -e $nonLTRraw; # allow empty file
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
`mkdir $genome.EDTA.combine` unless -e "$genome.EDTA.combine" && -d "$genome.EDTA.combine";

# enter the combine folder for EDTA processing
chdir "$genome.EDTA.combine";
`cp ../$LTRraw $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";
`cp ../$nonLTRraw $genome.nonLTR.raw.fa` unless -s "$genome.nonLTR.raw.fa";
`cp ../$TIRraw $genome.TIR.raw.fa` unless -s "$genome.TIR.raw.fa";
`cp ../$Helitronraw $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";


##################################
######  define subroutines  ######
##################################

# purify $TE2 contaminants in $TE1
sub Purifier() {
	my ($TE1, $TE2, $mindiff) = ($_[0], $_[1], $_[2]);
	`perl $TE_purifier -TE1 $TE1 -TE2 $TE2 -t $threads -mindiff $mindiff`;
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 1 -minrm 1 -trf 0 -f $TE1-$TE2.fa > $TE1.HQ`;
	}


#################################
###### Advance filtering ######
#################################

## predefined variables
my $LTR = "$genome.LTR.raw.fa";
my $nonLTR = "$genome.nonLTR.raw.fa";
my $TIR = "$genome.TIR.raw.fa";
my $HEL = "$genome.Helitron.raw.fa";
my $HEL_TM = "TE_Trimmer_consensus_merged.fasta_flex5";

## Purge contaminants
# purify LTR
&Purifier("$LTR", "$TIR", $mindiff_LTR);
&Purifier("$LTR.HQ", "$HEL", $mindiff_LTR);
`mv $LTR.HQ.HQ $LTR.HQ`;

#&Purifier("$HEL_TM", "$LTR.HQ", $mindiff_HEL);
#&Purifier("$HEL", "$HEL_TM.HQ", $mindiff_HEL);
#`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 1 -minrm 1 -trf 0 -f $HEL.masked > $HEL.HQ`;
#`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -minrm 1 -trf 0 -f $HEL.masked > $HEL.HQ`;
#&Purifier("$HEL.HQ", "$TIR", 1);
#`cat "$HEL.HQ.HQ" "$HEL_TM.HQ" > "$HEL.HQ"`;

#exit;
#&Purifier("$HEL.HQ", "$TIR", $mindiff_HEL);

# purify Helitron
&Purifier("$HEL", "$HEL_TM", $mindiff_HEL);
`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -minrm 1 -trf 0 -f $HEL.masked > $HEL.raw`;
`cat $HEL_TM >> $HEL.raw`;
&Purifier("$HEL.raw", "$TIR", $mindiff_HEL);
&Purifier("$HEL.raw.HQ", "$LTR", $mindiff_HEL);
`mv $HEL.raw.HQ.HQ $HEL.HQ`;

# purify TIR
&Purifier("$TIR", "$LTR", $mindiff_TIR);
&Purifier("$TIR.HQ", "$HEL.HQ", $mindiff_TIR);
#&Purifier("$TIR.HQ", "$HEL", $mindiff_TIR);
`mv $TIR.HQ.HQ $TIR.HQ`;

if (0){
# clean LTRs in nonLTRs
$err = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $LTR.HQ $genome.nonLTR.raw.fa 2>&1`;
if ($err !~ /done/) {
	`cp $genome.nonLTR.raw.fa $genome.nonLTR.raw.fa.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: $1/s;
	print STDERR "\n$err\n";
	}
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.nonLTR.raw.fa.masked > $genome.nonLTR.raw.fa.cln`;
}

# clean nonLTRs in LTRs
if (-s "$genome.nonLTR.raw.fa"){
	$err = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.nonLTR.raw.fa $LTR 2>&1`;
	if ($err !~ /done/) {
        	`cp $LTR $LTR.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: $1/s;
	        print STDERR "\n$err\n";
        	}
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $LTR.masked > $LTR.cln`;
	} else {
		`cp $LTR $LTR.cln`;
	}


## clean LTRs, Helitrons, and nonLTRs in TIRs
#`cat $TIR.HQ $HEL.HQ | perl -nle 's/>/\\n>/g unless /^>/; print \$_' > $genome.TIR.Helitron.fa.stg1.raw`;
#`cat $LTR.HQ $genome.nonLTR.raw.fa.cln > $genome.LTR.nonLTR.fa`;
#`cat $LTR.HQ $genome.nonLTR.raw.fa > $genome.LTR.nonLTR.fa`;
#`cat $LTR.HQ $genome.nonLTR.raw.fa $HEL_TM.HQ > $genome.LTR.nonLTR.Hel.fa`;
`cat $LTR.HQ $genome.nonLTR.raw.fa $HEL.HQ > $genome.LTR.nonLTR.Hel.fa`;
#$err = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.nonLTR.fa $genome.TIR.Helitron.fa.stg1.raw 2>&1`;
$err = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.nonLTR.Hel.fa $TIR.HQ 2>&1`;
if ($err !~ /done/) {
	#`cp $genome.TIR.Helitron.fa.stg1.raw $genome.TIR.Helitron.fa.stg1.raw.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: $1/s;
	`cp $TIR.HQ $TIR.HQ.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: $1/s;
	print STDERR "\n$err\n";
	}
	#`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.Helitron.fa.stg1.raw.masked > $genome.TIR.Helitron.fa.stg1.raw.cln`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $TIR.HQ.masked > $TIR.HQ.cln`;
`cat $HEL.HQ $TIR.HQ.cln > $genome.TIR.Helitron.fa.stg1.raw.cln`;

## cluster TIRs and Helitrons and make stg1 raw library
`perl $cleanup_nested -in $genome.TIR.Helitron.fa.stg1.raw.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast`;
#`perl $cleanup_nested -in $TIR.HQ.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast`;
`cat $LTR $genome.nonLTR.raw.fa $genome.TIR.Helitron.fa.stg1.raw.cln.cln > $genome.EDTA.stg1.raw.cln`;
#`cat $LTR $genome.nonLTR.raw.fa $genome.TIR.Helitron.fa.stg1.raw.cln.cln > $genome.EDTA.stg1.raw.cln`;
#`cat $LTR.cln $genome.nonLTR.raw.fa $genome.TIR.Helitron.fa.stg1.raw.cln.cln > $genome.EDTA.stg1.raw.cln`;
#`cat $LTR.cln $genome.nonLTR.raw.fa $HEL.HQ $TIR.HQ.cln.cln > $genome.EDTA.fa.stg1`;

## clean up the folder
#`rm *.masked *.ori.out *.nhr *.nin *.nsq 2>/dev/null`;

chdir '..';
