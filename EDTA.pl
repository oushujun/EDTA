#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

my $version = "v1.1";
#v1.0 05/31/2019
#v1.1 06/05/2019

print "
########################################################
##### Extensive de-novo TE Annotator (EDTA) $version    ####
##### Shujun Ou (shujun.ou.1\@gmail.com)             ####
########################################################
\n\n\n";

## Input: $genome
## Output: $genome.TElib.fa

my $usage = "\nGenerates a structure-based high-quality TE library
	perl EDTA.pl [options]
		-genome	[File]	The genome FASTA
		-step	[all|filter|final]	Specify which steps you want to run EDTA.
															all: run the entire pipeline (default)
															filter: start from raw TEs to the end.
															final: start from filtered TEs to finalizing the run.
		-repeatmodeler [path]	The directory containing RepeatModeler (default: read from ENV)
		-repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
		-blast [path]	The directory containing BLASTx and BLASTn (default: read from ENV)
		-trf [path]	The directory containing TRF (default: included in this package)
		-threads	[int]	Number of theads to run this script (default: 4)
		-help|-h	Display this help info
\n";

# pre-defined
my $genome = '';
my $step = "ALL";
my $threads = 4;
my $script_path = $FindBin::Bin;
my $EDTA_raw = "$script_path/EDTA_raw.pl";
my $EDTA_process = "$script_path/EDTA_process.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $repeatmodeler = "";
my $repeatmasker = "";
my $blast = "";
my $trf = "";

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$step = uc $ARGV[$k+1] if /^-step$/i and $ARGV[$k+1] !~ /^-/;
	$repeatmodeler = $ARGV[$k+1] if /^-repeatmodeler$/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker$/i and $ARGV[$k+1] !~ /^-/;
	$blast = $ARGV[$k+1] if /^-blast$/i and $ARGV[$k+1] !~ /^-/;
	$trf = $ARGV[$k+1] if /^-trf$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
	}

my $date=`date`;
chomp ($date);
print "$date\tDependency checking: ";

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "The script EDTA_raw.pl is not found in $EDTA_raw!\n" unless -s $EDTA_raw;
die "The script EDTA_process.pl is not found in $EDTA_process!\n" unless -s $EDTA_process;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;

#makeblastdb
$blast=`which makeblastdb 2>/dev/null` if $blast eq '';
$blast=~s/makeblastdb\n//;
die "makeblastdb is not exist in the BLAST+ path $blast!\n" unless -X "${blast}makeblastdb";
#blastn
$blast=`which blastn 2>/dev/null` if $blast eq '';
$blast=~s/blastn\n//;
die "blastn is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastn";
#blastx
$blast=`which blastx 2>/dev/null` if $blast eq '';
$blast=~s/blastx\n//;
die "blastx is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastx";
#RepeatMasker
my $rand=int(rand(1000000));
$repeatmasker=`which RepeatMasker 2>/dev/null` if $repeatmasker eq '';
$repeatmasker=~s/RepeatMasker\n//;
die "RepeatMasker is not exist in the RepeatMasker path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";
`cp $script_path/database/dummy060817.fa ./dummy060817.fa.$rand`;
my $RM_test=`${repeatmasker}RepeatMasker -e ncbi -q -pa 1 -no_is -norna -nolow dummy060817.fa.$rand -lib dummy060817.fa.$rand 2>/dev/null`;
die "The RMblast engine is not installed in RepeatMasker!\n" unless $RM_test=~s/done//gi;
`rm dummy060817.fa.$rand*`;
#trf
$trf="$script_path/bin/TRF/trf409.legacylinux64" if $trf eq ''; #path to the trf program
`$trf 2>/dev/null`;
$trf="$script_path/bin/TRF/trf409.macosx" if $?==32256;
`$trf 2>/dev/null`;
die "Error: No Tandem Repeat Finder is working on the current system.
						Both trf409.macosx and trf409.legacylinux64 were tested, and failed.
						Please report it to https://github.com/oushujun/LTR_retriever/issues" if $?==32256;
#cd-hit-est
#$cdhitpath=`which cd-hit-est 2>/dev/null` if $cdhitpath eq '';
#$cdhitpath=~s/cd-hit-est\n//;
#die "cd-hit-est is not exist in the CDHIT path $cdhitpath!\n" if (!(-X "${cdhitpath}cd-hit-est") and $cdhit);
#die "The path of CDHIT is not specified!\n" unless -X "${cdhitpath}cd-hit-est";

print "All passed!\n";

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

goto $step;

if (0){

ALL:
# Get raw TE candidates
$date=`date`;
chomp ($date);
print "$date\tObtain raw TE libraries using various structure-based programs: ";
`perl $EDTA_raw -genome $genome -threads $threads`;
die "ERROR: Raw LTR results not found in $genome.EDTA.raw/$genome.LTR.raw.fa" unless -s "$genome.EDTA.raw/$genome.LTR.raw.fa";
die "ERROR: Raw TIR results not found in $genome.EDTA.raw/$genome.TIR.raw.fa" unless -s "$genome.EDTA.raw/$genome.TIR.raw.fa";
die "ERROR: Raw MITE results not found in $genome.EDTA.raw/$genome.MITE.raw.fa" unless -s "$genome.EDTA.raw/$genome.MITE.raw.fa";
die "ERROR: Raw Helitron results not found in $genome.EDTA.raw/$genome.Helitron.raw.fa" unless -s "$genome.EDTA.raw/$genome.Helitron.raw.fa";

FILTER:
# Filter raw TE candidates and the make stage 1 library
$date=`date`;
chomp ($date);
print "$date\tPerform EDTA basic and advcanced filterings for raw TE candidates and generate the stage 1 library: ";
`perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -tir $genome.EDTA.raw/$genome.TIR.raw.fa -mite $genome.EDTA.raw/$genome.MITE.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.raw.fa -repeatmasker $repeatmasker -blast $blast -threads $threads`;
die "ERROR: Stage 1 library not found in $genome.LTR.TIR.Helitron.fa.stg1" unless -s "$genome.LTR.TIR.Helitron.fa.stg1";

FINAL:
my $date=`date`;
chomp ($date);
print "$date\tPerform EDTA final steps to generate a non-redundant comprehensive TE library: ";

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
