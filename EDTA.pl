#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $version = "v1.7.3";
#v1.0 05/31/2019
#v1.1 06/05/2019
#v1.2 06/16/2019
#v1.3 07/20/2019
#v1.4 08/07/2019
#v1.5 08/14/2019
#v1.6 11/09/2019
#v1.7 12/25/2019

print "
########################################################
##### Extensive de-novo TE Annotator (EDTA) $version  ####
##### Shujun Ou (shujun.ou.1\@gmail.com)             ####
########################################################
\n\n\n";

## Input: $genome
## Output: $genome.EDTA.TElib.fa

my $usage = "\nThis is the Extensive de-novo TE Annotator that generates a high-quality structure-based TE library. Usage:
	perl EDTA.pl [options]
		--genome	[File]	The genome FASTA
		--species [Rice|Maize|others]	Specify the species for identification of TIR candidates. Default: others
		--step	[all|filter|final|anno] Specify which steps you want to run EDTA.
						all: run the entire pipeline (default)
						filter: start from raw TEs to the end.
						final: start from filtered TEs to finalizing the run.
						anno: perform whole-genome annotation/analysis after TE library construction.
		--overwrite	[0|1]	If previous raw TE results are found, decide to overwrite (1, rerun) or not (0, default).
		--cds	[File]	Provide a FASTA file containing the coding sequence (no introns, UTRs, nor TEs) of this genome or its close relative.
		--curatedlib	[File]	Provided a curated library to keep consistant naming and classification for known TEs.
					TEs in this file will be trusted 100%, so please ONLY provide MANUALLY CURATED ones.
					This option is not mandatory. It's totally OK if no file is provided (default).
		--sensitive	[0|1]	Use RepeatModeler to identify remaining TEs (1) or not (0, default).
					This step is very slow and MAY help to recover some TEs.
		--anno	[0|1]	Perform (1) or not perform (0, default) whole-genome TE annotation after TE library construction.
		--evaluate [0|1]	Evaluate (1) classification consistency of the TE annotation. (-anno 1 required). Default: 0.
					This step is slow and does not affect the annotation result.
		--exclude	[File]	Exclude bed format regions from TE annotation. Default: undef. (-anno 1 required).
		--force	[0|1]	When no confident TE candidates are found: 0, interrupt and exit (default); 1, use rice TEs to continue.
		--repeatmodeler [path]	The directory containing RepeatModeler (default: read from ENV)
		--repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
		--blast [path]	The directory containing BLASTx and BLASTn (default: read from ENV)
		--check_dependencies Check if dependencies are fullfiled and quit
		--trf [path]	The directory containing TRF (default: read from ENV)
		--threads|-t	[int]	Number of theads to run this script (default: 4)
		--help|-h	Display this help info
\n";

# pre-defined
my $genome = '';
my $check_dependencies = undef;
my $species = "others";
my $step = "ALL";
my $overwrite = 0; #0, no rerun. 1, rerun even old results exist.
my $HQlib = '';
my $cds = ''; #a fasta file containing cds of this genome.
my $sensitive = 0; #0, will not run RepeatModeler to get remaining TEs (default). 1, run RepeatModeler
my $anno = 0; #0, will not annotate whole-genome TE (default). 1, annotate with RepeatMasker
my $evaluate = 0; #1 will evaluate the consistancy of the TE annotation
my $exclude = ''; #a bed file exclude from TE annotation
my $force = 0; #if there is no confident TE found in EDTA_raw, 1 will use rice TEs as raw lib, 0 will error and interrupt.
my $threads = 4;
my $script_path = $FindBin::Bin;
my $EDTA_raw = "$script_path/EDTA_raw.pl";
my $EDTA_process = "$script_path/EDTA_processI.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $cleanup_TE = "$script_path/util/cleanup_TE.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $count_nested = "$script_path/util/count_nested.pl";
my $count_base = "$script_path/util/count_base.pl";
my $make_masked = "$script_path/util/make_masked.pl";
my $make_gff3 = "$script_path/util/make_gff3_with_RMout.pl";
my $protlib = "$script_path/database/alluniRefprexp082813";
my $rice_LTR = "$script_path/database/rice6.9.5.liban.LTR";
#my $rice_nonLTR = "$script_path/database/rice6.9.5.liban.nonLTR";
my $rice_TIR = "$script_path/database/rice6.9.5.liban.TIR";
my $rice_helitron = "$script_path/database/rice6.9.5.liban.Helitron";
my $rename_TE = "$script_path/util/rename_TE.pl";
my $rename_RM = "$script_path/util/rename_RM_TE.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
my $buildSummary = "$script_path/util/buildSummary.pl"; #modified from RepeatMasker. Robert M. Hubley (rhubley@systemsbiology.org)
my $filter_gff = "$script_path/util/filter_gff.pl";
my $bed2gff = "$script_path/util/bed2gff.pl";
my $gff2bed = "$script_path/util/gff2bed.pl";
my $get_frag = "$script_path/util/get_frag.pl";
my $keep_nest = "$script_path/util/keep_nest.pl";
my $combine_overlap = "$script_path/util/combine_overlap.pl";
my $reclassify = "$script_path/util/classify_by_lib_RM.pl";
my $rename_by_list = "$script_path/util/rename_by_list.pl";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $TEsorter = "";
my $mdust = "";
my $GRF = "";
my $repeatmodeler = "";
my $repeatmasker = "";
my $blast = "";
my $trf = "";
my $beta2 = 0; #0, beta2 is not ready. 1, try it out.
my $help = undef;

# read parameters
if ( !GetOptions( 'genome=s'            => \$genome,
                  'species=s'           => \$species,
                  'step=s'              => \$step,
                  'overwrite=i'         => \$overwrite,
                  'curatedlib=s'        => \$HQlib,
                  'cds=s'                => \$cds,
                  'sensitive=s'          => \$sensitive,
		  'anno=i'               => \$anno,
		  'evaluate=i'           => \$evaluate,
		  'exclude=s'            => \$exclude,
		  'force=i'              => \$force,
		  'tesorter=s'           => \$TEsorter,
		  'repeatmodeler=s'      => \$repeatmodeler,
		  'repeatmasker=s'       => \$repeatmasker,
		  'blast=s'              => \$blast,
		  'protlib=s'            => \$protlib,
  		  'trf=s'                => \$trf,
		  'threads|t=i'          => \$threads,
		  'check_dependencies!'  => \$check_dependencies,
		  'help|h!'              => \$help ) )

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($help) {
    pod2usage( { -verbose => 0,
                 -exitval => 0,
                 -message => "$usage\n" } );
}

if ( ! $genome and (! $check_dependencies) ){
    pod2usage( {
           -message => "At least 1 parameter mandatory:\n1) Input fasta file: --genome\n".
           "$usage\n\n",
           -verbose => 0,
           -exitval => 2 } );
}


# check bolean
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n";}
if ($sensitive != 0 and $sensitive != 1){ die "The expected value for the sensitive parameter is 0 or 1!\n";}
if ($anno != 0 and $anno != 1){ die "The expected value for the anno parameter is 0 or 1!\n";}
if ($evaluate != 0 and $evaluate != 1){ die "The expected value for the evaluate parameter is 0 or 1!\n";}
if ($force != 0 and $force != 1){ die "The expected value for the force parameter is 0 or 1!\n";}

my $date=`date`;
chomp ($date);
print "$date\tDependency checking:\n";

# check files and dependencies
die "The script EDTA_raw.pl is not found in $EDTA_raw!\n" unless -s $EDTA_raw;
die "The script EDTA_processF.pl is not found in $EDTA_process!\n" unless -s $EDTA_process;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script cleanup_TE.pl is not found in $cleanup_TE!\n" unless -s $cleanup_TE;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script count_nested.pl is not found in $count_nested!\n" unless -s $count_nested;
die "The script count_base.pl is not found in $count_base!\n" unless -s $count_base;
die "The script make_masked.pl is not found in $make_masked!\n" unless -s $make_masked;
die "The script make_gff3_with_RMout.pl is not found in $make_gff3!\n" unless -s $make_gff3;
die "The protein-coding sequence library is not found in $protlib!\n" unless -s $protlib;
die "The rice LTR sequence library is not found in $rice_LTR!\n" unless -s $rice_LTR;
#die "The rice nonLTR sequence library is not found in $rice_nonLTR!\n" unless -s $rice_nonLTR;
die "The rice TIR sequence library is not found in $rice_TIR!\n" unless -s $rice_TIR;
die "The rice Helitron sequence library is not found in $rice_helitron!\n" unless -s $rice_helitron;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script buildSummary.pl is not found in $buildSummary!\n" unless -s $buildSummary;
die "The script filter_gff.pl is not found in $filter_gff!\n" unless -s $filter_gff;
die "The script bed2gff.pl is not found in $bed2gff!\n" unless -s $bed2gff;
die "The script gff2bed.pl is not found in $gff2bed!\n" unless -s $gff2bed;
die "The script get_frag.pl is not found in $get_frag!\n" unless -s $get_frag;
die "The script keep_nest.pl is not found in $keep_nest!\n" unless -s $keep_nest;
die "The script combine_overlap.pl is not found in $combine_overlap!\n" unless -s $combine_overlap;
die "The script classify_by_lib_RM.pl is not found in $reclassify!\n" unless -s $reclassify;
die "The script rename_by_list.pl is not found in $rename_by_list!\n" unless -s $rename_by_list;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
print "All dpendencies satisfied\n";
if ($check_dependencies){
	print "exit\n";
	exit;
}

# makeblastdb, blastn, blastx
$blast=`which makeblastdb 2>/dev/null` if $blast eq '';
$blast=~s/makeblastdb\n//;
$blast="$blast/" if $blast ne '' and $blast !~ /\/$/;
die "makeblastdb/blastn/blastx is not exist in the BLAST+ path $blast!\n" unless -X "${blast}makeblastdb";
die "blastn is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastn";
die "blastx is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastx";
# TEsorter
$TEsorter=`which TEsorter 2>/dev/null` if $TEsorter eq '';
$TEsorter=~s/TEsorter\n//;
$TEsorter="$TEsorter/" if $TEsorter ne '' and $TEsorter !~ /\/$/;
die "TEsorter is not exist in the TEsorter path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
# RepeatMasker
my $rand=int(rand(1000000));
$repeatmasker=`which RepeatMasker 2>/dev/null` if $repeatmasker eq '';
$repeatmasker=~s/RepeatMasker\n//;
$repeatmasker="$repeatmasker/" if $repeatmasker ne '' and $repeatmasker !~ /\/$/;
die "RepeatMasker is not exist in the RepeatMasker path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";
`cp $script_path/database/dummy060817.fa ./dummy060817.fa.$rand`;
my $RM_test=`${repeatmasker}RepeatMasker -e ncbi -q -pa 1 -no_is -norna -nolow dummy060817.fa.$rand -lib dummy060817.fa.$rand 2>/dev/null`;
die "The RMblast engine is not installed in RepeatMasker!\n" unless $RM_test=~s/done//gi;
`rm dummy060817.fa.$rand*`;
# RepeatModeler
$repeatmodeler=`which RepeatModeler 2>/dev/null` if $repeatmodeler eq '';
$repeatmodeler=~s/RepeatModeler\n//;
$repeatmodeler="$repeatmodeler/" if $repeatmodeler ne '' and $repeatmodeler !~ /\/$/;
die "RepeatModeler is not exist in the RepeatModeler path $repeatmodeler!\n" unless -X "${repeatmodeler}RepeatModeler";
# trf
$trf=`which trf 2>/dev/null` if $trf eq '';
$trf=~s/\n$//;
`$trf 2>/dev/null`;
die "Error: No Tandem Repeat Finder is working on the current system.
	Please report it to https://github.com/oushujun/EDTA/issues" if $?==32256;
die "\n\tTandem Repeat Finder not found!\n\n" unless $trf ne '';
# GRF
$GRF = "$script_path/bin/GenericRepeatFinder/bin/grf-main" if $GRF eq ''; #default path to the GRF program
`$GRF 2>/dev/null`;
die "Error: The Generic Repeat Finder (GRF) is not working on the current system.
	Please reinstall it in $GRF following instructions in https://github.com/bioinfolabmu/GenericRepeatFinder.
	If you continus to encounter this issue, please report it to https://github.com/oushujun/EDTA/issues\n" if $?==32256;
# mdust
$mdust=`which mdust 2>/dev/null` if $mdust eq '';
$mdust=~s/mdust\n//;
$mdust="$mdust/" if $mdust ne '' and $mdust !~ /\/$/;
die "mdust is not working on the current system. Please reinstall it in this folder $mdust.
	If you continus to encounter this issue, please report it to https://github.com/oushujun/EDTA/issues\n" unless -X "${mdust}mdust";

print "\t\t\t\tAll passed!\n";

# make a softlink to the user-provided files
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

# check $HQlib
if ($HQlib ne ''){
	if (-s $HQlib){
		print "\n\tCustom library $HQlib is provided via -curatedlib. Please make sure this is a manually curated library but not machine generated.\n\n";
		my $HQlib_file = basename($HQlib);
		`ls -s $HQlib $HQlib_file` unless -e $HQlib_file;
		$HQlib = $HQlib_file;
		} else {
		die "\n\tERROR: The custom library $HQlib you specified is not found!\n\n";
		}
	}

if ($cds ne ''){
	if (-s $cds){
		print "\n\tA CDS file is provided via -cds. Please make sure there is no TE-related sequences in this file.\n\n";
		my $cds_file = basename($cds);
		`ls -s $cds $cds_file` unless -e $cds_file;
		$cds = $cds_file;
		} else {
		die "\n\tERROR: The CDS file $cds you specified is not found!\n\n";
		}
	}

if ($exclude ne ''){
	if (-s $exclude){
		print "\n\tA BED file is provided via -exclude. Regions specified by this file will be excluded from TE annotation and masking.\n\n";
		my $exclude_file = basename($exclude);
		`ls -s $exclude $exclude_file ` unless -e $exclude_file;
		$exclude = $exclude_file;
		} else {
		die "\n\tERROR: The exclusion BED file $exclude you specified is not found!\n\n";
		}
	}

goto $step;


##################################################
####### Get raw LTR/TIR/Helitron candidates ######
##################################################

ALL:

# report status
$date=`date`;
chomp ($date);
print "$date\tObtain raw TE libraries using various structure-based programs: \n";

# Get raw TE candidates
`perl $EDTA_raw -genome $genome -overwrite $overwrite -species $species -threads $threads -mdust $mdust -blastplus $blast -tesorter $TEsorter`;

chdir "$genome.EDTA.raw";

# Force to use rice TEs when raw.fa is empty
if ($force eq 1){
	`cp $rice_LTR $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";
#	`cp $rice_LTR $genome.nonLTR.raw.fa` unless -s "$genome.nonLTR.raw.fa";
	`cp $rice_TIR $genome.TIR.raw.fa` unless -s "$genome.TIR.raw.fa";
	`cp $rice_helitron $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";
	}

# check results and report status
die "ERROR: Raw LTR results not found in $genome.EDTA.raw/$genome.LTR.raw.fa" unless -s "$genome.LTR.raw.fa";
die "ERROR: Raw TIR results not found in $genome.EDTA.raw/$genome.TIR.raw.fa" unless -s "$genome.TIR.raw.fa";
die "ERROR: Raw Helitron results not found in $genome.EDTA.raw/$genome.Helitron.raw.fa" unless -s "$genome.Helitron.raw.fa";

# combine intact TEs
`cat $genome.LTR.intact.fa $genome.TIR.intact.fa $genome.Helitron.intact.fa > $genome.EDTA.intact.fa`;
`cat $genome.LTR.intact.fa.gff3 $genome.TIR.intact.fa.gff $genome.Helitron.intact.fa.gff | perl $gff2bed - structural > $genome.EDTA.intact.bed`;
`perl $bed2gff $genome.EDTA.intact.bed`;
`mv $genome.EDTA.intact.bed.gff $genome.EDTA.intact.gff`;
`cp $genome.EDTA.intact.gff ../`;

$date=`date`;
chomp ($date);
print "$date\tObtain raw TE libraries finished.
\t\t\t\tAll intact TEs found by EDTA: $genome.EDTA.intact.fa\t$genome.EDTA.intact.gff\n\n";
chdir "..";


##################################################
####### Filter LTR/TIR/Helitron candidates #######
##################################################

FILTER:

# report status
$date=`date`;
chomp ($date);
print "$date\tPerform EDTA advcance filtering for raw TE candidates and generate the stage 1 library: \n\n";

# Filter raw TE candidates and the make stage 1 library
`perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -tir $genome.EDTA.raw/$genome.TIR.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.raw.fa -repeatmasker $repeatmasker -blast $blast -threads $threads -protlib $protlib`;

# check results and report status
die "ERROR: Stage 1 library not found in $genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1" unless -s "$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1";
$date=`date`;
chomp ($date);
print "$date\tEDTA advcance filtering finished.\n\n";


#####################################
###### Final TE/SINE/LINE scan ######
#####################################

FINAL:

# report status
$date=`date`;
chomp ($date);
print "$date\tPerform EDTA final steps to generate a non-redundant comprehensive TE library:\n\n";

# Make the final working directory
`mkdir $genome.EDTA.final` unless -e "$genome.EDTA.final" && -d "$genome.EDTA.final";
chdir "$genome.EDTA.final";
`rm -rf $genome.* 2>/dev/null`;
`cp ../$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1 ./`;
`cp ../$cds ./` if $cds ne '';
`cp ../$HQlib ./` if $HQlib ne '';
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.fa ./`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.bed ./`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.gff ./`;

# identify remaining TEs in the genome
if ($sensitive == 1){
	print "\t\t\t\tUse RepeatModeler to identify any remaining TEs that are missed by structure-based methods.\n\n";
	# RepeatMask the genome with the cleanned stage 1 library
	`ln -s ../$genome $genome` unless -e $genome;
	`${repeatmasker}RepeatMasker -pa $threads -qq -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1 $genome 2>/dev/null`;

	# Scan the repeatmasked genome with RepeatModeler for any remaining TEs
	`${repeatmodeler}BuildDatabase -name $genome.masked -engine ncbi $genome.masked`;
	`${repeatmodeler}RepeatModeler -engine ncbi -pa $threads -database $genome.masked 2>/dev/null`;
	`rm $genome.masked.nhr $genome.masked.nin $genome.masked.nnd $genome.masked.nni $genome.masked.nog $genome.masked.nsq`;

	# rename RepeatModeler candidates and make stage 2 library
	`perl $rename_RM RM_*/consensi.fa.classified > $genome.RepeatModeler.raw.fa`;
	if (-s "$genome.RepeatModeler.raw.fa"){
		`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1 $genome.RepeatModeler.raw.fa 2>/dev/null`;
		`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.RepeatModeler.raw.fa.masked > $genome.RepeatModeler.fa.stg1`;
		`cat $genome.RepeatModeler.fa.stg1 $genome.LTR.TIR.Helitron.fa.stg1 > $genome.LTR.TIR.Helitron.others.fa.stg2`;

		# clean up coding sequences in the stage 2 library
		`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.others.fa.stg2 -rmdnate 0 -rmline 0 -rmprot 1 -protlib $protlib -blast $blast -threads $threads`;
		} else {
		print "\t\t\t\tRepeatModeler is finished, but no consensi.fa.classified files found.\n\n";
		`cp $genome.LTR.TIR.Helitron.fa.stg1 $genome.LTR.TIR.Helitron.others.fa.stg2.clean`;
		}
	} else {
	print "\t\t\t\tSkipping the RepeatModeler step (-sensitive 0).\n\t\t\t\tRun EDTA.pl -step final -sensitive 1 if you want to use RepeatModeler.\n\n";
	`cp $genome.LTR.TIR.Helitron.fa.stg1 $genome.LTR.TIR.Helitron.others.fa.stg2.clean`;
	}

# rename file
`cp $genome.LTR.TIR.Helitron.others.fa.stg2.clean $genome.EDTA.raw.fa`;

if ($cds ne ''){
	# report status
	$date=`date`;
	chomp ($date);
	print "$date\tRemove CDS in the EDTA library:\n\n";

	# cleanup CDS with TEsorter
	`perl $cleanup_TE -cds $cds -minlen 300 -tesorter $TEsorter -repeatmasker $repeatmasker -t $threads -rawlib $genome.EDTA.raw.fa`;
	$cds = "$cds.mod.noTE";

	# remove cds in the EDTA library
	if (-s "$cds"){
		`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.raw.fa 2>/dev/null`;
		`${repeatmasker}RepeatMasker -pa $threads -qq -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.intact.fa 2>/dev/null`;
		`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.3 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.EDTA.raw.fa.masked > $genome.EDTA.raw.fa.cln`;
		`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.8 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 1 -f $genome.EDTA.intact.fa.masked > $genome.EDTA.intact.fa`;

		# remove gene seq in intact TEs
		if (-s "$genome.EDTA.intact.fa.masked.cleanup"){
			`grep -v -P "Only|head|tail" $genome.EDTA.intact.fa.masked.cleanup | awk '{if (\$2>=0.8) print \$1}' |sort -u | awk '{print "Parent\\t"\$1"\\nID\\t"\$1}' > $genome.EDTA.intact.fa.masked.cleanup.rmlist`;
			`perl $output_by_list 1 $genome.EDTA.intact.fa 2 $genome.EDTA.intact.fa.masked.cleanup.rmlist -ex -FA > $genome.EDTA.intact.fa.rmTE`;
			`mv $genome.EDTA.intact.fa.rmTE $genome.EDTA.intact.fa`; #update intact.fa

			`perl $filter_gff $genome.EDTA.intact.gff $genome.EDTA.intact.fa.masked.cleanup.rmlist > $genome.EDTA.intact.gff.new`;
			`perl -nle 'my \$id = \$1 if /=(repeat_region[0-9]+);/; print "Parent\t\$id" if defined \$id' $genome.EDTA.intact.gff.removed >> $genome.EDTA.intact.fa.masked.cleanup.rmlist`;
			`perl $filter_gff $genome.EDTA.intact.gff $genome.EDTA.intact.fa.masked.cleanup.rmlist > $genome.EDTA.intact.gff.new`;
			`mv $genome.EDTA.intact.gff.new $genome.EDTA.intact.gff`; #update intact.gff
			`perl $gff2bed $genome.EDTA.intact.gff structural > $genome.EDTA.intact.bed`; #update intact.bed
			}
		} else {
		print STDERR "\t\t\t\tWarning: No CDS left after clean up ($cds.mod.noTE empty). Will not clean CDS in the raw lib.\n\n";
		`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.cln`;
		}

	} else {
	print "\t\t\t\tSkipping the CDS cleaning step (-cds [File]) since no CDS file is provided or it's empty.\n\n";
	`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.cln`;
	}

# Final rounds of redundancy removal and make final EDTA library
`perl $cleanup_nested -in $genome.EDTA.raw.fa.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast 2>/dev/null`;

# rename all TEs in the EDTA library
`perl $rename_TE $genome.EDTA.raw.fa.cln.cln > $genome.EDTA.TElib.fa`;

# check results
die "ERROR: Final TE library not found in $genome.EDTA.TElib.fa" unless -s "$genome.EDTA.TElib.fa";
`cp $genome.EDTA.TElib.fa ../`;

if ($HQlib ne ''){
	# report status
	$date=`date`;
	chomp ($date);
	print "$date\tCombine the high-quality TE library $HQlib with the EDTA library:\n\n";

	# remove known TEs in the EDTA library
	`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib ../$HQlib $genome.EDTA.TElib.fa 2>/dev/null`;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 0 -f $genome.EDTA.TElib.fa.masked > $genome.EDTA.TElib.novel.fa`;
	`mv $genome.EDTA.TElib.fa $genome.EDTA.TElib.ori.fa`;
	`cat $HQlib $genome.EDTA.TElib.novel.fa > $genome.EDTA.TElib.fa`;
	`cp $genome.EDTA.TElib.novel.fa $genome.EDTA.TElib.fa ../`;

	# reclassify intact TEs with known TEs
	`${repeatmasker}RepeatMasker -pa $threads -qq -no_is -norna -nolow -div 40 -lib ../$HQlib $genome.EDTA.intact.fa 2>/dev/null`;
	`perl $reclassify -seq $genome.EDTA.intact.fa -RM $genome.EDTA.intact.fa.out`;
	`perl $rename_by_list $genome.EDTA.intact.bed $genome.EDTA.intact.fa.rename.list 1 > $genome.EDTA.intact.bed.rename`;
	`perl $bed2gff $genome.EDTA.intact.bed.rename`;
	`mv $genome.EDTA.intact.bed.rename.gff $genome.EDTA.intact.gff`; #update intact.gff
	`cp $genome.EDTA.intact.gff ../`; #replace the intact gff that has no lib family info
	}

# report status
$date=`date`;
chomp ($date);
print "$date\tEDTA final stage finished! You may check out:
		\t\tThe final EDTA TE library: $genome.EDTA.TElib.fa\n";
print "		\t\tFamily names of intact TEs have been updated by $HQlib: $genome.EDTA.intact.gff\n" if $HQlib ne '';
print "\tComparing to the curated library you provided, this are the novel TEs EDTA found: $genome.EDTA.TElib.novel.fa
	The high-quality library you provided has been incorporated into the final library: $genome.EDTA.TElib.fa\n\n" if $HQlib ne '';
chdir "..";


#####################################
###### Post-library annotation ######
#####################################

ANNO:
if ($anno == 1){
	# report status
	$date=`date`;
	chomp ($date);
	print "$date\tPerform post-EDTA analysis for whole-genome annotation:\n\n";

	# Make the post-library annotation working directory
	`mkdir $genome.EDTA.anno` unless -e "$genome.EDTA.anno" && -d "$genome.EDTA.anno";
	chdir "$genome.EDTA.anno";
	`rm -rf $genome.* 2>/dev/null`;
	`cp ../$genome.EDTA.final/$genome.EDTA.TElib.fa ./`;
	`cp ../$exclude ./` if $exclude ne '';
	`ln -s ../$genome $genome` unless -e $genome;

	# annotate TEs using RepeatMasker
	`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.EDTA.TElib.fa $genome 2>/dev/null`;
	die "ERROR: RepeatMasker results not found in $genome.out!\n\n" unless -s "$genome.out" or -s "$genome.mod.out";

	# exclude regions from TE annotation and make whole-genome TE annotation
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv 30 -minscore 300 -minlen 80 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
	`mv $genome.out.new $genome.EDTA.RM.out`;
	`perl $make_gff3 $genome.EDTA.RM.out`;
	`mv $genome.EDTA.RM.out.gff $genome.EDTA.RM.gff`;

	# combine homology-based and strutrual-based annotation
	`perl $gff2bed $genome.EDTA.RM.gff homology > $genome.EDTA.RM.bed`;
	`cp ../$genome.EDTA.final/$genome.EDTA.intact.bed ./`;
	`perl $combine_overlap $genome.EDTA.intact.bed $genome.EDTA.intact.bed.cmb 5`;
	`perl $keep_nest $genome.EDTA.intact.bed $genome.EDTA.RM.bed $threads`;
	`perl $keep_nest $genome.EDTA.RM.bed $genome.EDTA.intact.bed.cmb $threads`;
	`sort -suV $genome.EDTA.intact.bed-$genome.EDTA.RM.bed $genome.EDTA.RM.bed-$genome.EDTA.intact.bed.cmb > $genome.EDTA.TEanno.bed`;
	`perl $bed2gff $genome.EDTA.TEanno.bed`;
	`mv $genome.EDTA.TEanno.bed.gff $genome.EDTA.TEanno.gff`;

	# make summary table for the annotation
	my $genome_size = `perl $count_base $genome`;
	$genome_size = (split /\s+/, $genome_size)[1] - (split /\s+/, $genome_size)[2];
	`perl -nle 'my (\$chr, \$s, \$e, undef, \$supfam, undef, \$anno)=(split); next if \$supfam=~/target_site_duplication|long_terminal_repeat/i; \$anno=~s/ID=//; \$anno=~s/;.*//; print "10000 0.001 0.001 0.001 \$chr \$s \$e NA NA \$anno \$supfam"' $genome.EDTA.TEanno.bed > $genome.EDTA.TEanno.out`;
	`perl $buildSummary -maxDiv 40 -genome_size $genome_size $genome.EDTA.TEanno.out > $genome.EDTA.TEanno.sum 2>/dev/null`;
	my $tot_TE = `grep Total $genome.EDTA.TEanno.sum|grep %|awk '{print \$4}'`;
	chomp $tot_TE;

	# make low-threshold masked genome for MAKER
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
	`mv $genome.new.masked $genome.MAKER.masked`;
	my $maker_TE = `perl $count_base $genome.MAKER.masked`;
	$maker_TE = (split /\s+/, $maker_TE)[-1];
	$maker_TE = sprintf("%.2f%%", $maker_TE*100);

	# check results and report status
	die "ERROR: TE annotation results not found in $genome.EDTA.TEanno.gff!\n\n" unless -s "$genome.EDTA.TEanno.gff";
	print "ERROR: The masked genome for MAKER annotation is not found in $genome.MAKER.masked!\n\n" unless -s "$genome.MAKER.masked";
	$date=`date`;
	chomp ($date);
	print "$date\tTE annotation using the EDTA library has finished! Check out:\n";
	print "\t\t\t\tWhole-genome TE annotation (total TE: $tot_TE): $genome.EDTA.TEanno.gff\n";
	print "\t\t\t\tWhole-genome TE annotation summary: $genome.EDTA.TEanno.sum\n";
	print "\t\t\t\tLow-threshold TE masking for MAKER gene annotation (masked: $maker_TE): $genome.MAKER.masked\n\n";
	`cp $genome.MAKER.masked $genome.EDTA.TEanno.gff $genome.EDTA.TEanno.sum ../`;

	# evaluate the annotation consistency
	if ($evaluate == 1){
		# report status
		$date=`date`;
		chomp ($date);
		print "$date\tEvaluate the level of inconsistency for whole-genome TE annotation (slow step):\n\n";

		# extract whole-genome TE and perform all-v-all blast, then summarize the results
		`awk '{if (\$5~/[0-9]+/ && \$1>300 && \$7-\$6>80) print \$11"\t"\$5":"\$6".."\$7}' $genome.EDTA.TEanno.out | perl $call_seq - -C $genome > $genome.EDTA.TE.fa`;
		`perl $cleanup_nested -in $genome.EDTA.TE.fa -threads $threads -minlen 80 -miniden 80 -cov 0.95 -blastplus $blast 2>/dev/null`;
		`for i in nested all redun; do perl $count_nested -in $genome.EDTA.TE.fa.stat -cat \$i > $genome.EDTA.TE.fa.stat.\$i.sum; done`;

		# check results and report status
		die "ERROR: TE annotation stats results not found in $genome.EDTA.TE.fa.stat!\n\n" unless -s "$genome.EDTA.TE.fa.stat";
		$date=`date`;
		chomp ($date);
		print "$date\tEvaluation of TE annotation finished! Check out these files:\n
				Overall: $genome.EDTA.TE.fa.stat.all.sum
				Nested: $genome.EDTA.TE.fa.stat.nested.sum
				Non-nested: $genome.EDTA.TE.fa.stat.redun.sum\n\n";
		}

	}
