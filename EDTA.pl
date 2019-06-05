#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

########################################################
##### Extensive de-novo TE Annotator (EDTA) v1.0    ####
##### Shujun Ou (shujun.ou.1@gmail.com, 05/31/2019) ####
########################################################

## Input: $genome
## Output: $genome.TElib.fa

my $usage = '';

# pre-defined
my $genome = '';
my $threads = 4;
my $script_path = $FindBin::Bin;
my $EDTA_raw = "$script_path/EDTA_raw.pl";
my $EDTA_process = "$script_path/EDTA_process.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $repeatmodeler = " ";
my $repeatmasker = " ";
my $blast = " ";

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$repeatmodeler = $ARGV[$k+1] if /^-repeatmodeler/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blast = $ARGV[$k+1] if /^-blast/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
	}

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "The script EDTA_raw.pl is not found in $EDTA_raw!\n" unless -s $EDTA_raw;
die "The script EDTA_process.pl is not found in $EDTA_process!\n" unless -s $EDTA_process;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

if (0){
# Get raw TE candidates
`perl $EDTA_raw -genome $genome -threads $threads`;

# Filter raw TE candidates and the make stage 1 library
`perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -tir $genome.EDTA.raw/$genome.TIR.raw.fa -mite $genome.EDTA.raw/$genome.MITE.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.raw.fa -repeatmasker $repeatmasker -blast $blast -threads $threads`;

# Make the final working directory
`mkdir $genome.EDTA.final` unless -e "$genome.EDTA.final" && -d "$genome.EDTA.final";
chdir "$genome.EDTA.final";

# clean up LINE retrotransposases in the stage 1 library
`cp ../$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1 $genome.LTR.TIR.Helitron.fa.stg1`;
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.fa.stg1 -rmdnate 0 -rmline 1 -rmprot 0 -blast $blast -threads $threads`;

# RepeatMask the genome with the cleanned stage 1 library
`ln -s ../$genome $genome` unless -e $genome;
`${repeatmasker}RepeatMasker -pa $threads -qq -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1.clean $genome`;
}

chdir "$genome.EDTA.final"; #tst

#####################################
###### Final TE/SINE/LINE scan ######
#####################################

`${repeatmodeler}BuildDatabase -name $genome.masked -engine ncbi $genome.masked`;
`${repeatmodeler}RepeatModeler -engine ncbi -pa $threads -database $genome.masked`;

exit;
# rename RepeatModeler candidates and make stage 2 library
#`cat TBD > $genome.RepeatModeler.raw.fa`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.RepeatModeler.raw.fa > $genome.RepeatModeler.raw.fa.stg0`;
`cat $genome.RepeatModeler.raw.fa.stg0 $genome.LTR.TIR.Helitron.fa.stg1.clean > $genome.LTR.TIR.Helitron.others.fa.stg2`;

# clean up coding sequences in the stage 2 library
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.others.fa.stg2 -rmdnate 0 -rmline 0 -rmprot 1 -blast $blast -threads $threads`;

# final rounds of redundancy removal and make final library
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.others.fa.stg2.clean -threads $threads -minlene 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.others.fa.stg2.cln`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.others.fa.stg2.cln -threads $threads -minlene 80 -cov 0.95 -blastplus $blast > $genome.TElib.fa`;

# copy the library one folder up
`cp $genome.TElib.fa ../`;



