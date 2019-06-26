#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

my $version = "v1.2";
#v1.0 05/31/2019
#v1.1 06/05/2019
#v1.2 06/16/2019

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
		-step	[all|filter|final] Specify which steps you want to run EDTA.
						all: run the entire pipeline (default)
						filter: start from raw TEs to the end.
						final: start from filtered TEs to finalizing the run.
		-protlib [File] Protein-coding aa sequences to be removed from TE candidates. (default lib: alluniRefprexp082813 (plant))
					You may use uniprot_sprot database available from here:
					ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
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
my $protlib = "$script_path/database/alluniRefprexp082813";
my $rename_TE = "$script_path/util/rename_TE.pl";
my $GRF = "";
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
	$protlib = $ARGV[$k+1] if /^-protlib/i and $ARGV[$k+1] !~ /^-/;
	$trf = $ARGV[$k+1] if /^-trf$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
	}

my $date=`date`;
chomp ($date);
print "$date\tDependency checking:\n";

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "The script EDTA_raw.pl is not found in $EDTA_raw!\n" unless -s $EDTA_raw;
die "The script EDTA_process.pl is not found in $EDTA_process!\n" unless -s $EDTA_process;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The protein-coding sequence library is not found in $protlib!\n" unless -s $protlib;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;

# makeblastdb
$blast=`which makeblastdb 2>/dev/null` if $blast eq '';
$blast=~s/makeblastdb\n//;
die "makeblastdb is not exist in the BLAST+ path $blast!\n" unless -X "${blast}makeblastdb";
# blastn
$blast=`which blastn 2>/dev/null` if $blast eq '';
$blast=~s/blastn\n//;
die "blastn is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastn";
# blastx
$blast=`which blastx 2>/dev/null` if $blast eq '';
$blast=~s/blastx\n//;
die "blastx is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastx";
# RepeatMasker
my $rand=int(rand(1000000));
$repeatmasker=`which RepeatMasker 2>/dev/null` if $repeatmasker eq '';
$repeatmasker=~s/RepeatMasker\n//;
die "RepeatMasker is not exist in the RepeatMasker path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";
`cp $script_path/database/dummy060817.fa ./dummy060817.fa.$rand`;
my $RM_test=`${repeatmasker}RepeatMasker -e ncbi -q -pa 1 -no_is -norna -nolow dummy060817.fa.$rand -lib dummy060817.fa.$rand 2>/dev/null`;
die "The RMblast engine is not installed in RepeatMasker!\n" unless $RM_test=~s/done//gi;
`rm dummy060817.fa.$rand*`;
# trf
$trf="$script_path/bin/TRF/trf409.legacylinux64" if $trf eq ''; #default path to the trf program
`$trf 2>/dev/null`;
$trf="$script_path/bin/TRF/trf409.macosx" if $?==32256;
`$trf 2>/dev/null`;
die "Error: No Tandem Repeat Finder is working on the current system.
	Both trf409.macosx and trf409.legacylinux64 were tested, and failed.
	Please report it to https://github.com/oushujun/EDTA/issues" if $?==32256;
# GRF
$GRF = "$script_path/bin/GenericRepeatFinder/bin/grf-main" if $GRF eq ''; #default path to the GRF program 
`$GRF 2>/dev/null`;
die "Error: The Generic Repeat Finder (GRF) is not working on the current system.
	Please reinstall it in $GRF following instructions in https://github.com/bioinfolabmu/GenericRepeatFinder.
	If you continus to encounter this issue, please report it to https://github.com/oushujun/EDTA/issues" if $?==32256;

#cd-hit-est
#$cdhitpath=`which cd-hit-est 2>/dev/null` if $cdhitpath eq '';
#$cdhitpath=~s/cd-hit-est\n//;
#die "cd-hit-est is not exist in the CDHIT path $cdhitpath!\n" if (!(-X "${cdhitpath}cd-hit-est") and $cdhit);
#die "The path of CDHIT is not specified!\n" unless -X "${cdhitpath}cd-hit-est";

print "\t\tAll passed!\n";

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

goto $step;


ALL:
# Get raw TE candidates
$date=`date`;
chomp ($date);
print "$date\tObtain raw TE libraries using various structure-based programs: \n";
`perl $EDTA_raw -genome $genome -threads $threads`;
die "ERROR: Raw LTR results not found in $genome.EDTA.raw/$genome.LTR.raw.fa" unless -e "$genome.EDTA.raw/$genome.LTR.raw.fa";
die "ERROR: Raw TIR results not found in $genome.EDTA.raw/$genome.TIR.raw.fa" unless -e "$genome.EDTA.raw/$genome.TIR.raw.fa";
die "ERROR: Raw MITE results not found in $genome.EDTA.raw/$genome.MITE.raw.fa" unless -e "$genome.EDTA.raw/$genome.MITE.raw.fa";
die "ERROR: Raw Helitron results not found in $genome.EDTA.raw/$genome.Helitron.raw.fa" unless -e "$genome.EDTA.raw/$genome.Helitron.raw.fa";
$date=`date`;
chomp ($date);
print "$date\tObtain raw TE libraries finished.\n";


FILTER:
# Filter raw TE candidates and the make stage 1 library
$date=`date`;
chomp ($date);
print "$date\tPerform EDTA basic and advcanced filterings for raw TE candidates and generate the stage 1 library: \n";
`perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -tir $genome.EDTA.raw/$genome.TIR.raw.fa -mite $genome.EDTA.raw/$genome.MITE.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.raw.fa -repeatmasker $repeatmasker -blast $blast -threads $threads -protlib $protlib`;
die "ERROR: Stage 1 library not found in $genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1" unless -s "$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1";
$date=`date`;
chomp ($date);
print "$date\tEDTA basic and advcanced filters finished.\n";


#####################################
###### Final TE/SINE/LINE scan ######
#####################################

FINAL:
$date=`date`;
chomp ($date);
print "$date\tPerform EDTA final steps to generate a non-redundant comprehensive TE library:\n";

# Make the final working directory
`mkdir $genome.EDTA.final` unless -e "$genome.EDTA.final" && -d "$genome.EDTA.final";
chdir "$genome.EDTA.final";
if (1){
`rm -rf *`;

# RepeatMask the genome with the cleanned stage 1 library
`cp ../$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1 ./`;
`ln -s ../$genome $genome` unless -e $genome;
`${repeatmasker}RepeatMasker -pa $threads -qq -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1 $genome 2>/dev/null`;

# Scan the repeatmasked genome with RepeatModeler for any remaining TEs
`${repeatmodeler}BuildDatabase -name $genome.masked -engine ncbi $genome.masked`;
`${repeatmodeler}RepeatModeler -engine ncbi -pa $threads -database $genome.masked 2>/dev/null`;
}

# rename RepeatModeler candidates and make stage 2 library
`cat RM_*/round-*/family-*fa | perl -nle \'print \$_ and next unless /^>/; my \$name=(split)[2]; print \">\$name\"\' > $genome.RepeatModeler.raw.fa`;
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1 $genome.RepeatModeler.raw.fa 2>/dev/null`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.RepeatModeler.raw.fa.masked > $genome.RepeatModeler.fa.stg1`;
`cat $genome.RepeatModeler.fa.stg1 $genome.LTR.TIR.Helitron.fa.stg1 > $genome.LTR.TIR.Helitron.others.fa.stg2`;

# clean up coding sequences in the stage 2 library
`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.others.fa.stg2 -rmdnate 0 -rmline 0 -rmprot 1 -protlib $protlib -blast $blast -threads $threads`;

# final rounds of redundancy removal and make final library
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.others.fa.stg2.clean -threads $threads -minlene 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.others.fa.stg2.cln`;
`perl $cleanup_nested -in $genome.LTR.TIR.Helitron.others.fa.stg2.cln -threads $threads -minlene 80 -cov 0.95 -blastplus $blast > $genome.LTR.TIR.Helitron.others.fa.stg2.cln2`;

# rename all TEs in the final library
`perl $rename_TE $genome.LTR.TIR.Helitron.others.fa.stg2.cln2 > $genome.TElib.fa`;

# check results
die "ERROR: Final TE library not found in $genome.TElib.fa" unless -s "$genome.TElib.fa";
`cp $genome.TElib.fa ../`;
$date=`date`;
chomp ($date);
print "$date\tEDTA final stage finished! Check out the final TE library: $genome.TElib.fa\n\n";


