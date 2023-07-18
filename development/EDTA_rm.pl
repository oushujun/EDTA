#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $version = "v2.1.4";
#v1.0 05/31/2019
#v1.1 06/05/2019
#v1.2 06/16/2019
#v1.3 07/20/2019
#v1.4 08/07/2019
#v1.5 08/14/2019
#v1.6 11/09/2019
#v1.7 12/25/2019
#v1.8 02/09/2020
#v1.9 07/24/2020
#v2.0 11/25/2021
#v2.1 10/10/2022

print "
########################################################
##### Extensive de-novo TE Annotator (EDTA) $version  ####
##### Shujun Ou (shujun.ou.1\@gmail.com)             ####
########################################################
\n\n\n";

## Input: $genome
## Output: $genome.EDTA.TElib.fa

my $usage = "\nThis is the Extensive de-novo TE Annotator that generates a high-quality
structure-based TE library. Usage:

perl EDTA.pl [options]
	--genome [File]		The genome FASTA file. Required.
	--species [Rice|Maize|others]	Specify the species for identification of TIR
					candidates. Default: others
	--step [all|filter|final|anno]	Specify which steps you want to run EDTA.
					all: run the entire pipeline (default)
					filter: start from raw TEs to the end.
					final: start from filtered TEs to finalizing the run.
					anno: perform whole-genome annotation/analysis after
						TE library construction.
	--overwrite [0|1]	If previous raw TE results are found, decide to overwrite
				(1, rerun) or not (0, default).
	--cds [File]	Provide a FASTA file containing the coding sequence (no introns,
			UTRs, nor TEs) of this genome or its close relative.
	--curatedlib [File]	Provided a curated library to keep consistant naming and
				classification for known TEs. TEs in this file will be
				trusted 100%, so please ONLY provide MANUALLY CURATED ones.
				This option is not mandatory. It's totally OK if no file is
				provided (default).
	--sensitive [0|1]	Use RepeatModeler to identify remaining TEs (1) or not (0,
				default). This step is slow but MAY help to recover some TEs.
	--anno [0|1]	Perform (1) or not perform (0, default) whole-genome TE annotation
			after TE library construction.
	--rmout	[File]	Provide your own homology-based TE annotation instead of using the
			EDTA library for masking. File is in RepeatMasker .out format. This
			file will be merged with the structural-based TE annotation. (--anno 1
			required). Default: use the EDTA library for annotation.
	--evaluate [0|1]	Evaluate (1) classification consistency of the TE annotation.
				(--anno 1 required). Default: 0. This step is slow and does
				not change the annotation result.
	--exclude [File]	Exclude regions (bed format) from TE masking in the MAKER.masked
				output. Default: undef. (--anno 1 required).
	--force	[0|1]	When no confident TE candidates are found: 0, interrupt and exit
			(default); 1, use rice TEs to continue.
	--u [float]	Neutral mutation rate to calculate the age of intact LTR elements.
			Intact LTR age is found in this file: *EDTA_raw/LTR/*.pass.list.
			Default: 1.3e-8 (per bp per year, from rice).
	--repeatmodeler [path]	The directory containing RepeatModeler (default: read from ENV)
	--repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
	--check_dependencies Check if dependencies are fullfiled and quit
	--threads|-t [int]	Number of theads to run this script (default: 4)
	--debug	 [0|1]	Retain intermediate files (default: 0)
	--help|-h 	Display this help info
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
my $rmout = ''; #a RM .out file for custom homology-based annotation.
my $evaluate = 0; #1 will evaluate the consistancy of the TE annotation
my $exclude = ''; #a bed file exclude from TE annotation
my $force = 0; #if there is no confident TE found in EDTA_raw, 1 will use rice TEs as raw lib, 0 will error and interrupt.
my $miu = 1.3e-8; #mutation rate, per bp per year, from rice
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
my $filter_gff = "$script_path/util/filter_gff3.pl";
my $RMout2bed = "$script_path/util/RMout2bed.pl";
my $bed2gff = "$script_path/util/bed2gff.pl";
my $gff2bed = "$script_path/util/gff2bed.pl";
my $get_frag = "$script_path/util/get_frag.pl";
my $keep_nest = "$script_path/util/keep_nest.pl";
my $combine_overlap = "$script_path/util/combine_overlap.pl";
my $split_overlap = "$script_path/util/split_overlap.pl";
my $reclassify = "$script_path/util/classify_by_lib_RM.pl";
my $rename_by_list = "$script_path/util/rename_by_list.pl";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $format_TElib = "$script_path/util/format_TElib.pl";
my $LTR_retriever = "";
my $genometools = "";
my $repeatmodeler = "";
my $repeatmasker = "";
my $TEsorter = "";
my $blastplus = "";
my $mdust = "";
my $trf = "";
my $GRF = "";

my $beta2 = 0; #0, beta2 is not ready. 1, developer mode.
#my $reanno = 0; #0, use existing whole-genome RM results (beta); 1, de novo Repeatmasker using the EDTA library (default)
my $debug = 0;
my $help = undef;

# read parameters
if ( !GetOptions( 'genome=s'            => \$genome,
                  'species=s'           => \$species,
                  'step=s'              => \$step,
                  'overwrite=i'         => \$overwrite,
                  'curatedlib=s'        => \$HQlib,
                  'cds=s'                => \$cds,
                  'protlib=s'            => \$protlib,
                  'sensitive=s'          => \$sensitive,
		  'anno=i'               => \$anno,
		  'rmout=s'              => \$rmout,
		  'evaluate=i'           => \$evaluate,
		  'exclude=s'            => \$exclude,
		  'force=i'              => \$force,
		  'u=s'                  => \$miu,
		  'repeatmodeler=s'      => \$repeatmodeler,
		  'repeatmasker=s'       => \$repeatmasker,
		  'tesorter=s'           => \$TEsorter,
		  'blast=s'              => \$blastplus,
		  'threads|t=i'          => \$threads,
		  'check_dependencies!'  => \$check_dependencies,
                  'debug=i'              => \$debug,
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

if ( (! -s $genome) and (! $check_dependencies) ){
    pod2usage( {
           -message => "At least 1 parameter is required:\n1) Input fasta file: --genome\n".
           "$usage\n\n",
           -verbose => 0,
           -exitval => 2 } );
	}

# check bolean
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"}
if ($sensitive != 0 and $sensitive != 1){ die "The expected value for the sensitive parameter is 0 or 1!\n"}
if ($anno != 0 and $anno != 1){ die "The expected value for the anno parameter is 0 or 1!\n"}
if ($evaluate != 0 and $evaluate != 1){ die "The expected value for the evaluate parameter is 0 or 1!\n"}
if ($force != 0 and $force != 1){ die "The expected value for the force parameter is 0 or 1!\n"}
if ($miu !~ /[0-9\.e\-]+/){ die "The expected value for the u parameter is float value without units!\n"}
if ($debug != 0 and $debug != 1){ die "The expected value for the debug parameter is 0 or 1!\n"}
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"}

# define RepeatMasker -pa parameter
my $rm_threads = int($threads/4);

chomp (my $date = `date`);
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
die "The script filter_gff3.pl is not found in $filter_gff!\n" unless -s $filter_gff;
die "The script RMout2bed.pl is not found in $RMout2bed!\n" unless -s $RMout2bed;
die "The script bed2gff.pl is not found in $bed2gff!\n" unless -s $bed2gff;
die "The script gff2bed.pl is not found in $gff2bed!\n" unless -s $gff2bed;
die "The script get_frag.pl is not found in $get_frag!\n" unless -s $get_frag;
die "The script keep_nest.pl is not found in $keep_nest!\n" unless -s $keep_nest;
die "The script combine_overlap.pl is not found in $combine_overlap!\n" unless -s $combine_overlap;
die "The script split_overlap.pl is not found in $split_overlap!\n" unless -s $split_overlap;
die "The script classify_by_lib_RM.pl is not found in $reclassify!\n" unless -s $reclassify;
die "The script rename_by_list.pl is not found in $rename_by_list!\n" unless -s $rename_by_list;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;

# GenomeTools
chomp ($genometools=`which gt 2>/dev/null`) if $genometools eq '';
$genometools =~ s/\s+$//;
$genometools = dirname($genometools) unless -d $genometools;
$genometools="$genometools/" if $genometools ne '' and $genometools !~ /\/$/;
die "Error: gt is not found in the genometools path $genometools!\n" unless -X "${genometools}gt";
# LTR_retriever
chomp ($LTR_retriever=`which LTR_retriever 2>/dev/null`) if $LTR_retriever eq '';
$LTR_retriever =~ s/\s+$//;
$LTR_retriever = dirname($LTR_retriever) unless -d $LTR_retriever;
$LTR_retriever="$LTR_retriever/" if $LTR_retriever ne '' and $LTR_retriever !~ /\/$/;
die "Error: LTR_retriever is not found in the LTR_retriever path $LTR_retriever!\n" unless -X "${LTR_retriever}LTR_retriever";
# RepeatMasker
my $rand=int(rand(1000000));
chomp ($repeatmasker=`which RepeatMasker 2>/dev/null`) if $repeatmasker eq '';
$repeatmasker =~ s/\s+$//;
$repeatmasker = dirname($repeatmasker) unless -d $repeatmasker;
$repeatmasker="$repeatmasker/" if $repeatmasker ne '' and $repeatmasker !~ /\/$/;
die "Error: RepeatMasker is not found in the RepeatMasker path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";
`cp $script_path/database/dummy060817.fa ./dummy060817.fa.$rand`;
my $RM_test=`${repeatmasker}RepeatMasker -e ncbi -q -pa 1 -no_is -norna -nolow dummy060817.fa.$rand -lib dummy060817.fa.$rand 2>/dev/null`;
die "Error: The RMblast engine is not installed in RepeatMasker!\n" unless $RM_test=~s/done//gi;
`rm dummy060817.fa.$rand* 2>/dev/null`;
# RepeatModeler
chomp ($repeatmodeler=`which RepeatModeler 2>/dev/null`) if $repeatmodeler eq '';
$repeatmodeler =~ s/\s+$//;
$repeatmodeler = dirname($repeatmodeler) unless -d $repeatmodeler;
$repeatmodeler="$repeatmodeler/" if $repeatmodeler ne '' and $repeatmodeler !~ /\/$/;
die "Error: RepeatModeler is not found in the RepeatModeler path $repeatmodeler!\n" unless -X "${repeatmodeler}RepeatModeler";
# makeblastdb, blastn, blastx
chomp ($blastplus=`which makeblastdb 2>/dev/null`) if $blastplus eq '';
$blastplus =~ s/\s+$//;
$blastplus = dirname($blastplus) unless -d $blastplus;
$blastplus="$blastplus/" if $blastplus ne '' and $blastplus !~ /\/$/;
die "Error: makeblastdb is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}makeblastdb";
die "Error: blastn is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";
die "Error: blastx is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastx";
# TEsorter
chomp ($TEsorter=`which TEsorter 2>/dev/null`) if $TEsorter eq '';
$TEsorter =~ s/\s+$//;
$TEsorter = dirname($TEsorter) unless -d $TEsorter;
$TEsorter="$TEsorter/" if $TEsorter ne '' and $TEsorter !~ /\/$/;
die "Error: TEsorter is not found in the TEsorter path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
# mdust
chomp ($mdust=`which mdust 2>/dev/null`) if $mdust eq '';
$mdust =~ s/\s+$//;
$mdust = dirname($mdust) unless -d $mdust;
$mdust = "$mdust/" if $mdust ne '' and $mdust !~ /\/$/;
die "Error: mdust is not found in the mdust path $mdust!\n" unless -X "${mdust}mdust";
# trf
chomp ($trf=`which trf 2>/dev/null`) if $trf eq '';
$trf=~s/\n$//;
`$trf 2>/dev/null`;
die "Error: Tandem Repeat Finder is not found in the TRF path $trf!\n" if $?==32256;
# GRF
chomp ($GRF = `which grf-main 2>/dev/null`) if $GRF eq '';
$GRF =~ s/\n$//;
`$GRF 2>/dev/null`;
die "Error: The Generic Repeat Finder (GRF) is not found in the GRF path: $GRF\n" if $?==32256;

print "\t\t\t\tAll passed!\n\n";
exit if $check_dependencies;

# make a softlink to the user-provided files
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

# check if duplicated sequences found
my $id_mode = 0; #record the mode of id conversion.
my $id_len = `grep \\> $genome|perl -ne 'chomp; s/>//g; my \$len=length \$_; \$max=\$len if \$max<\$len; print "\$max\\n"'`; #find out the longest sequence ID length in the genome
$id_len =~ s/\s+$//;
$id_len = (split /\s+/, $id_len)[-1];
my $raw_id = `grep \\> $genome|wc -l`;
my $old_id = `grep \\> $genome|sort -u|wc -l`;
if ($raw_id > $old_id){
	chomp ($date = `date`);
	die "$date\tERROR: Identical sequence IDs found in the provided genome! Please resolve this issue and try again.\n";
	}

# remove sequence annotations (content after the first space in sequence names) and replace special characters with _
`perl -nle 'my \$info=(split)[0]; \$info=~s/[\\~!@#\\\$%\\^&\\*\\(\\)\\+\\\-\\=\\?\\[\\]\\{\\}\\:;",\\<\\/\\\\\|]+/_/g; \$info=~s/_+/_/g; print \$info' $genome > $genome.mod`;

# try to shortern sequences
my $id_len_max = 13; # allowed longest length of a sequence ID in the input file
if ($id_len > $id_len_max){
	chomp ($date = `date`);
	print "$date\tThe longest sequence ID in the genome contains $id_len characters, which is longer than the limit ($id_len_max)\n";
	print "\t\t\t\tTrying to reformat seq IDs...\n\t\t\t\tAttempt 1...\n";
	`perl -lne 'chomp; if (s/^>+//) {s/^\\s+//; \$_=(split)[0]; s/(.{1,$id_len_max}).*/>\$1/g;} print "\$_"' $genome.mod > $genome.temp`;
	my $new_id = `grep \\> $genome.temp|sort -u|wc -l`;
	chomp ($date = `date`);
	if ($old_id == $new_id){
		$id_mode = 1;
		`mv $genome.temp $genome.mod`;
		print "$date\tSeq ID conversion successful!\n\n";
		} else {
		print "\t\t\t\tAttempt 2...\n";
		`perl -ne 'chomp; if (/^>/) {\$_=">\$1" if /([0-9]+)/;} print "\$_\n"' $genome.mod > $genome.temp`;
		$new_id = `grep \\> $genome.temp|sort -u|wc -l`;
		if ($old_id == $new_id){
			$id_mode = 2;
			`mv $genome.temp $genome.mod`;
			print "$date\tSeq ID conversion successful!\n\n";
			} else {
			`rm $genome.temp 2>/dev/null`;
			die "$date\tERROR: Fail to convert seq IDs to <= $id_len_max characters! Please provide a genome with shorter seq IDs.\n\n";
			}
		}
	}
$genome = "$genome.mod";

# check $HQlib
if ($HQlib ne ''){
	if (-s $HQlib){
		print "\tA custom library $HQlib is provided via --curatedlib. Please make sure this is a manually curated library but not machine generated.\n\n";
		chomp ($HQlib = `realpath $HQlib`);
		my $HQlib_file = basename($HQlib);
		`ln -s $HQlib $HQlib_file` unless -e $HQlib_file;
		$HQlib = $HQlib_file;
		} else {
		die "\tERROR: The custom library $HQlib you specified is not found!\n\n";
		}
	}

if ($cds ne ''){
	if (-s $cds){
		print "\tA CDS file $cds is provided via --cds. Please make sure this is the DNA sequence of coding regions only.\n\n";
		chomp ($cds = `realpath $cds`);
		my $cds_file = basename($cds);
		`ln -s $cds $cds_file` unless -e $cds_file;
		$cds = $cds_file;
		} else {
		die "\tERROR: The CDS file $cds you specified is not found!\n\n";
		}
	}

if ($rmout ne ''){
	if (-s $rmout){
		print "\tA RepeatMasker .out file $rmout is provided via --rmout.\n\n";
		chomp ($rmout = `realpath $rmout`);
		} else {
		die "\tERROR: The RepeatMasker .out file $rmout you specified is not found!\n\n";
		}
	}

if ($exclude ne ''){
	if (-s $exclude){
		print "\tA BED file is provided via --exclude. Regions specified by this file will be excluded from TE annotation and masking.\n\n";
		my $exclude_file = basename($exclude);
		`ln -s $exclude $exclude_file ` unless -e $exclude_file;
		$exclude = $exclude_file;
		} else {
		die "\tERROR: The exclusion BED file $exclude you specified is not found!\n\n";
		}
	}

$step = uc $step;
goto $step;


##################################################
####### Get raw LTR/TIR/Helitron candidates ######
##################################################

ALL:

# report status
chomp ($date = `date`);
print "$date\tObtain raw TE libraries using various structure-based programs: \n";

# Get raw TE candidates
`perl $EDTA_raw --genome $genome --overwrite $overwrite --species $species --u $miu --threads $threads --genometools $genometools --ltrretriever $LTR_retriever --blastplus $blastplus --tesorter $TEsorter --GRF $GRF --trf_path $trf --repeatmasker $repeatmasker --convert_seq_name 0`;

chdir "$genome.EDTA.raw";

# Force to use rice TEs when raw.fa is empty
if ($force eq 1){
	`cp $rice_LTR $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";
#	`cp $rice_LTR $genome.nonLTR.raw.fa` unless -s "$genome.nonLTR.raw.fa";
	`cp $rice_TIR $genome.TIR.raw.fa` unless -s "$genome.TIR.raw.fa";
	`cp $rice_helitron $genome.Helitron.raw.fa` unless -s "$genome.Helitron.raw.fa";
	}

# check results and report status
die "ERROR: Raw LTR results not found in $genome.EDTA.raw/$genome.LTR.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of intact LTRs in your genome. Consider to use the --force 1 parameter to overwrite this check\n" unless -s "$genome.LTR.raw.fa";
die "ERROR: Raw TIR results not found in $genome.EDTA.raw/$genome.TIR.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of intact TIRs in your genome. Consider to use the --force 1 parameter to overwrite this check\n" unless -s "$genome.TIR.raw.fa";
die "ERROR: Raw Helitron results not found in $genome.EDTA.raw/$genome.Helitron.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of intact Helitrons in your genome. Consider to use the --force 1 parameter to overwrite this check\n" unless -s "$genome.Helitron.raw.fa";

# combine intact TEs
`cat $genome.LTR.intact.fa $genome.TIR.intact.fa $genome.Helitron.intact.fa > $genome.EDTA.intact.fa`;
`cat $genome.TIR.intact.bed $genome.Helitron.intact.bed | perl $bed2gff - TE_struc > $genome.EDTA.intact.gff3.raw`;
`cat $genome.LTR.intact.gff3 >> $genome.EDTA.intact.gff3.raw`;
`sort -sV -k1,1 -k4,4 $genome.EDTA.intact.gff3.raw | grep -v '^#' > $genome.EDTA.intact.gff3; rm $genome.EDTA.intact.gff3.raw`;
`cp $genome.EDTA.intact.gff3 ../`;

chomp ($date = `date`);
print "$date\tObtain raw TE libraries finished.
\t\t\t\tAll intact TEs found by EDTA: \n\t\t\t\t\t$genome.EDTA.intact.fa\n\t\t\t\t\t$genome.EDTA.intact.gff3\n\n";
chdir "..";


##################################################
####### Filter LTR/TIR/Helitron candidates #######
##################################################

FILTER:

# report status
chomp ($date = `date`);
print "$date\tPerform EDTA advance filtering for raw TE candidates and generate the stage 1 library: \n\n";

# remove existing results
`rm ./$genome.EDTA.combine/* 2>/dev/null` if $overwrite == 1;

# Filter raw TE candidates and the make stage 1 library
`perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -tir $genome.EDTA.raw/$genome.TIR.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.raw.fa -repeatmasker $repeatmasker -blast $blastplus -threads $threads -protlib $protlib`;

# check results, remove intermediate files, and report status
die "ERROR: Stage 1 library not found in $genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1" unless -s "$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1";
chdir "$genome.EDTA.combine";
`rm ./$genome.LTR.raw* ./$genome.TIR.raw* ./$genome.Helitron.raw* ./$genome.TIR.Helitro* ./$genome.LTR.TIR.Helitron.fa.stg1.* 2>/dev/null` unless $debug eq 1;
chdir "..";
chomp ($date = `date`);
print "$date\tEDTA advance filtering finished.\n\n";


#####################################
###### Final TE/SINE/LINE scan ######
#####################################

FINAL:

# report status
chomp ($date = `date`);
print "$date\tPerform EDTA final steps to generate a non-redundant comprehensive TE library:\n\n";

# Make the final working directory
`mkdir $genome.EDTA.final` unless -e "$genome.EDTA.final" && -d "$genome.EDTA.final";
chdir "$genome.EDTA.final";
`rm ./* 2>/dev/null` if $overwrite == 1;
`cp ../$genome.EDTA.combine/$genome.LTR.TIR.Helitron.fa.stg1 ./`;
`cp ../$cds ./` if $cds ne '';
`cp ../$HQlib ./` if $HQlib ne '';
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.fa ./$genome.EDTA.intact.fa.raw`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.gff3 ./`;
`cp ../$exclude ./` if $exclude ne '';

###### developer notes
# $rmout may become useless



# identify remaining TEs in the genome
if ($sensitive == 1){
	print "\t\t\t\tUse RepeatModeler to identify any remaining TEs that are missed by structure-based methods.\n\n";
	# RepeatMask the genome with the cleanned stage 1 library
	`ln -s ../$genome $genome` unless -e $genome;
	#	if ($rmout ne ''){
		#print STDERR "$date\tA RepeatMasker result file $rmout is provided! Will use this file without running RepeatMasker.\n\n";
		#`perl $make_masked -genome $genome -rmout $rmout -maxdiv 40 -minscore 300 -minlen 80 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
		#`mv $genome.new.masked $genome.masked`;
		#		} else {
			#`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -qq -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1 $genome 2>/dev/null`;
		#		}

	chomp ($date = `date`);
	if ($overwrite eq 0 and -s "$genome.RM.consensi.fa"){
		print STDERR "$date\tExisting RepeatModeler result file $genome.RM.consensi.fa found!\n\t\t\t\tWill keep this file without rerunning this module.\n\t\t\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
		} else {
		`rm -rf ./RM_*/consensi.fa 2>/dev/null`;
		# Scan the repeatmasked genome with RepeatModeler for any remaining TEs
		`${repeatmodeler}BuildDatabase -name $genome -engine ncbi $genome`;
#		`${repeatmodeler}BuildDatabase -name $genome.masked -engine ncbi $genome.masked`;
		`${repeatmodeler}RepeatModeler -engine ncbi -pa $rm_threads -database $genome 2>/dev/null`;
#		`rm $genome.masked.nhr $genome.masked.nin $genome.masked.nnd $genome.masked.nni $genome.masked.nog $genome.masked.nsq 2>/dev/null`;
		`rm $genome.nhr $genome.nin $genome.nnd $genome.nni $genome.nog $genome.nsq 2>/dev/null`;
		`awk '{print \$1}' RM_*/consensi.fa.classified > $genome.RM.consensi.fa`;
		}

	# filter and reclassify RepeatModeler candidates with TEsorter and make stage 2 library
	if (-s "$genome.RM.consensi.fa"){
		
		`${TEsorter}TEsorter $genome.RM.consensi.fa -p $threads`;
		`perl $rename_RM $genome.RM.consensi.fa.rexdb.cls.lib > $genome.RepeatModeler.raw.fa`;
		my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.Helitron.fa.stg1 $genome.RepeatModeler.raw.fa 2>/dev/null`;
		`cp $genome.RepeatModeler.raw.fa $genome.RepeatModeler.raw.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
		`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.RepeatModeler.raw.fa.masked > $genome.RepeatModeler.fa.stg1`;
		`cat $genome.LTR.TIR.Helitron.fa.stg1 $genome.RepeatModeler.fa.stg1 > $genome.LTR.TIR.Helitron.others.fa.stg2`;

		# clean up coding sequences in the stage 2 library
		`perl $cleanup_proteins -seq $genome.LTR.TIR.Helitron.others.fa.stg2 -rmdnate 0 -rmline 0 -rmprot 1 -protlib $protlib -blast $blastplus -threads $threads`;
		} else {
		print "\t\t\t\tRepeatModeler is finished, but no consensi.fa files found.\n\n";
		`cp $genome.LTR.TIR.Helitron.fa.stg1 $genome.LTR.TIR.Helitron.others.fa.stg2.clean`;
		}
	} else {
	print "\t\t\t\tSkipping the RepeatModeler step (--sensitive 0).\n\t\t\t\tRun EDTA.pl --step final --sensitive 1 if you want to use RepeatModeler.\n\n";
	`cp $genome.LTR.TIR.Helitron.fa.stg1 $genome.LTR.TIR.Helitron.others.fa.stg2.clean`;
	}

# rename file
`cp $genome.LTR.TIR.Helitron.others.fa.stg2.clean $genome.EDTA.raw.fa`;

if ($cds ne ''){
	# report status
	chomp ($date = `date`);

	# cleanup TE-related sequences in the CDS file with TEsorter
	print "$date\tClean up TE-related sequences in the CDS file with TEsorter:\n\n";
	`perl $cleanup_TE -cds $cds -minlen 300 -tesorter $TEsorter -repeatmasker $repeatmasker -t $threads -rawlib $genome.EDTA.raw.fa`;
	`rm ./$cds ./$cds.code.r* 2>/dev/null` unless $debug eq 1;
	$cds = "$cds.code.noTE";

	# remove cds-related sequences in the EDTA library
	print "\t\t\t\tRemove CDS-related sequences in the EDTA library:\n\n";
	if (-s "$cds"){
		my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.raw.fa 2>/dev/null`;
		`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
		$rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.intact.fa.raw 2>/dev/null`;
		`cp $genome.EDTA.intact.fa.raw $genome.EDTA.intact.fa.raw.masked` if $rm_status =~ /No repetitive sequences were detected/i;
		`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.3 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.EDTA.raw.fa.masked > $genome.EDTA.raw.fa.cln`;
		`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.8 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 0 -f $genome.EDTA.intact.fa.raw.masked > $genome.EDTA.intact.fa.rmCDS`;

		# remove gene seq in intact TEs
		if (-s "$genome.EDTA.intact.fa.raw.masked.cleanup"){
			`grep -v -P "Only|head|tail" $genome.EDTA.intact.fa.raw.masked.cleanup | awk '{if (\$2>=0.8) print \$1}' |sort -u | awk '{print "Name\\t"\$1"\\nParent\\t"\$1"\\nID\\t"\$1}' > $genome.EDTA.intact.fa.raw.masked.cleanup.rmlist`;
			`perl $output_by_list 1 $genome.EDTA.intact.fa.raw 2 $genome.EDTA.intact.fa.raw.masked.cleanup.rmlist -ex -FA > $genome.EDTA.intact.fa`; #update intact.fa
			`perl $filter_gff $genome.EDTA.intact.gff3 $genome.EDTA.intact.fa.raw.masked.cleanup.rmlist > $genome.EDTA.intact.gff3.new`;
			`perl -nle 'my \$id = \$1 if /=(repeat_region[0-9]+);/; print "Parent\t\$id\nName\t\$id" if defined \$id' $genome.EDTA.intact.gff3.removed >> $genome.EDTA.intact.fa.raw.masked.cleanup.rmlist`;
			`perl $filter_gff $genome.EDTA.intact.gff3 $genome.EDTA.intact.fa.raw.masked.cleanup.rmlist > $genome.EDTA.intact.gff3.new`;
			`mv $genome.EDTA.intact.gff3.new $genome.EDTA.intact.gff3`; #update intact.gff
			} else {
			`cp $genome.EDTA.intact.fa.rmCDS $genome.EDTA.intact.fa`;
			}
		} else {
		print STDERR "\t\t\t\tWarning: No CDS left after clean up ($cds.code.noTE empty). Will not clean CDS in the raw lib.\n\n";
		`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.cln`;
		`cp $genome.EDTA.intact.fa.raw $genome.EDTA.intact.fa`;
		}

	} else {
	print "\t\t\t\tSkipping the CDS cleaning step (--cds [File]) since no CDS file is provided or it's empty.\n\n";
	`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.cln`;
	`cp $genome.EDTA.intact.fa.raw $genome.EDTA.intact.fa`;
	}

# Final rounds of redundancy removal and make final EDTA library
`perl $cleanup_nested -in $genome.EDTA.raw.fa.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blastplus 2>/dev/null`;

# rename all TEs in the EDTA library
`perl $rename_TE $genome.EDTA.raw.fa.cln.cln > $genome.EDTA.TElib.fa`;
#`perl $rename_TE $genome.EDTA.raw.fa.cln.cln | perl $format_TElib - > $genome.EDTA.TElib.fa`;

# check results
die "ERROR: Final TE library not found in $genome.EDTA.TElib.fa" unless -s "$genome.EDTA.TElib.fa";
die "ERROR: Intact TE annotation not found in $genome.EDTA.intact.gff3" unless -s "$genome.EDTA.intact.gff3";
`cp $genome.EDTA.TElib.fa ../`;
`cp $genome.EDTA.intact.gff3 ../`;

if ($HQlib ne ''){
	# report status
	chomp ($date = `date`);
	print "$date\tCombine the high-quality TE library $HQlib with the EDTA library:\n\n";

	# remove known TEs in the EDTA library
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $HQlib $genome.EDTA.TElib.fa 2>/dev/null`;
	`cp $genome.EDTA.TElib.fa $genome.EDTA.TElib.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 0 -f $genome.EDTA.TElib.fa.masked > $genome.EDTA.TElib.novel.fa`;
	`mv $genome.EDTA.TElib.fa $genome.EDTA.TElib.ori.fa`;
	`cat $HQlib $genome.EDTA.TElib.novel.fa > $genome.EDTA.TElib.fa`;
	`cp $genome.EDTA.TElib.novel.fa $genome.EDTA.TElib.fa ../`;
	}

# reclassify intact TEs with the TE lib #113
`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.EDTA.TElib.fa $genome.EDTA.intact.fa 2>/dev/null` unless -s "$genome.EDTA.intact.fa.out" and $overwrite == 0;
die "ERROR: The masked file for $genome.EDTA.intact.fa is not found! The RepeatMasker annotation on this file may be failed. Please check the $genome.EDTA.TElib.fa file for sequence naming formats especially when you provide a library via --curatedlib.\n" unless -s "$genome.EDTA.intact.fa.out";
`perl $reclassify -seq $genome.EDTA.intact.fa -RM $genome.EDTA.intact.fa.out`;
`perl $rename_by_list $genome.EDTA.intact.gff3 $genome.EDTA.intact.fa.rename.list 1 > $genome.EDTA.intact.gff3.rename`;
`mv $genome.EDTA.intact.fa.rename $genome.EDTA.intact.fa`;
`mv $genome.EDTA.intact.gff3.rename $genome.EDTA.intact.gff3`;
`cp $genome.EDTA.intact.gff3 ../`; #replace the intact gff that has no lib family info

# remove intermediate files
`rm $genome.EDTA.intact.fa.raw.* $genome.EDTA.raw.fa.* $genome.EDTA.TElib.fa.* $genome.LTR.TIR.Helitron.fa.stg1.* $genome.masked *.cat.gz 2>/dev/null` if $debug eq 0;

# report status
chomp ($date = `date`);
print "$date\tEDTA final stage finished! You may check out:
		\t\tThe final EDTA TE library: $genome.EDTA.TElib.fa\n";
print "		\t\tFamily names of intact TEs have been updated by $HQlib: $genome.EDTA.intact.gff3\n" if $HQlib ne '';
print "\t\t\t\tComparing to the provided library, EDTA found these novel TEs: $genome.EDTA.TElib.novel.fa
	\t\t\tThe provided library has been incorporated into the final library: $genome.EDTA.TElib.fa\n\n" if $HQlib ne '';
chdir "..";


#####################################
###### Post-library annotation ######
#####################################

ANNO:
if ($anno == 1){
	# report status
	chomp ($date = `date`);
	print "$date\tPerform post-EDTA analysis for whole-genome annotation:\n\n";

	# Make the post-library annotation working directory
	`mkdir $genome.EDTA.anno` unless -e "$genome.EDTA.anno" && -d "$genome.EDTA.anno";
	chdir "$genome.EDTA.anno";
	`rm ./* 2>/dev/null` if $overwrite == 1;
	`cp ../$genome.EDTA.final/$genome.EDTA.TElib.fa ./`;
	`cp ../$genome.EDTA.final/$genome.EDTA.intact.gff3 ./`;
	`cp ../$exclude ./` if $exclude ne '';
	`ln -s ../$genome $genome` unless -e $genome;

	# annotate TEs using RepeatMasker
	if ($rmout ne ''){
		print STDERR "$date\tA RepeatMasker result file $rmout is provided! Will use this file without running RepeatMasker.\n\n";
		if (-e "$genome.out"){
			my $old_rmout = `ls -l $genome.out|perl -nle 'my (\$month, \$day, \$time) = (split)[6,7,8]; \$time =~ s/://; print "\${month}_\${day}_\$time"'`;
			chomp $old_rmout;
			print "\t\t\t\t$genome.out exists in the $genome.EDTA.anno folder, renamed file to ${genome}_$old_rmout.out\n\n";
			`mv $genome.out ${genome}_$old_rmout.out`;
			}
		`ln -s $rmout $genome.out`;
		} else {
		print STDERR "$date\tHomology-based annotation of TEs using $genome.EDTA.TElib.fa from scratch.\n\n";
		`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.EDTA.TElib.fa $genome 2>/dev/null`;
		}
	die "ERROR: RepeatMasker results not found in $genome.out!\n\n" unless -s "$genome.out" or -s "$genome.mod.out";

	# exclude regions from TE annotation and make whole-genome TE annotation
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv 30 -minscore 300 -minlen 80 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
	`mv $genome.out.new $genome.EDTA.RM.out`;
	`perl $RMout2bed $genome.EDTA.RM.out > $genome.EDTA.RM.bed`; # a regular enriched bed
	`perl $bed2gff $genome.EDTA.RM.bed TE_homo > $genome.EDTA.RM.gff3`;
	`perl $gff2bed $genome.EDTA.RM.gff3 homology > $genome.EDTA.RM.bed`; # add the last column to this bed

	# combine homology-based and strutrual-based annotation (partly overlapping)
	`perl $gff2bed $genome.EDTA.intact.gff3 structural > $genome.EDTA.intact.bed`;
	`perl $combine_overlap $genome.EDTA.intact.bed $genome.EDTA.intact.bed.cmb 5`;
	`perl $get_frag $genome.EDTA.RM.bed $genome.EDTA.intact.bed.cmb $threads`;
	`perl $keep_nest $genome.EDTA.intact.bed $genome.EDTA.RM.bed $threads`;
	`grep homology $genome.EDTA.intact.bed-$genome.EDTA.RM.bed > $genome.EDTA.intact.bed-$genome.EDTA.RM.bed.homo`;
	`sort -suV $genome.EDTA.intact.bed-$genome.EDTA.RM.bed.homo $genome.EDTA.RM.bed-$genome.EDTA.intact.bed.cmb > $genome.EDTA.homo.bed`;
	`perl $bed2gff $genome.EDTA.homo.bed TE_homo > $genome.EDTA.homo.gff3`;
	`cat $genome.EDTA.intact.gff3 $genome.EDTA.homo.gff3 > $genome.EDTA.TEanno.gff3.raw`;
	`grep -v '^#' $genome.EDTA.TEanno.gff3.raw | sort -sV -k1,1 -k4,4 | perl -0777 -ne '\$date=\`date\`; \$date=~s/\\s+\$//; print "##gff-version 3\\n##date \$date\\n##Identity: Sequence identity (0-1) between the library sequence and the target region.\\n##ltr_identity: Sequence identity (0-1) between the left and right LTR regions.\\n##tsd: target site duplication.\\n##seqid source sequence_ontology start end score strand phase attributes\\n\$_"' - > $genome.EDTA.TEanno.gff3`;
	`rm $genome.EDTA.TEanno.gff3.raw 2>/dev/null`;

	# make non-overlapping annotation
	`perl $gff2bed $genome.EDTA.TEanno.gff3 structural > $genome.EDTA.TEanno.bed`;
	`perl $split_overlap $genome.EDTA.TEanno.bed $genome.EDTA.TEanno.split.bed`;
	`perl $bed2gff $genome.EDTA.TEanno.split.bed > $genome.EDTA.TEanno.split.gff3`;

	# make summary table for the non-overlapping annotation
	`perl $count_base $genome > $genome.stats`;
	`perl -nle 'my (\$chr, \$s, \$e, \$anno, \$dir, \$supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 \$chr \$s \$e NA \$dir \$anno \$supfam"' $genome.EDTA.TEanno.split.bed > $genome.EDTA.TEanno.out`;
	`perl $buildSummary -maxDiv 40 -stats $genome.stats $genome.EDTA.TEanno.out > $genome.EDTA.TEanno.sum 2>/dev/null`;
	my $tot_TE = `grep Total $genome.EDTA.TEanno.sum|grep %|awk '{print \$4}'`;
	chomp $tot_TE;

	# make low-threshold masked genome for MAKER
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 1 -misschar N -threads $threads -exclude $exclude` unless -s "$genome.MAKER.masked" and $overwrite == 0;
	`mv $genome.new.masked $genome.MAKER.masked`;
	my $maker_TE = `perl $count_base $genome.MAKER.masked`;
	$maker_TE = (split /\s+/, $maker_TE)[3];
	$maker_TE = sprintf("%.2f%%", $maker_TE*100);

	# check results and report status
	die "ERROR: TE annotation results not found in $genome.EDTA.TEanno.gff3!\n\n" unless -s "$genome.EDTA.TEanno.gff3";
	print "ERROR: The masked genome for MAKER annotation is not found in $genome.MAKER.masked!\n\n" unless -s "$genome.MAKER.masked";
	chomp ($date = `date`);
	print "$date\tTE annotation using the EDTA library has finished! Check out:\n";
	print "\t\t\t\tWhole-genome TE annotation (total TE: $tot_TE): $genome.EDTA.TEanno.gff3\n";
	print "\t\t\t\tWhole-genome TE annotation summary: $genome.EDTA.TEanno.sum\n";
	print "\t\t\t\tLow-threshold TE masking for MAKER gene annotation (masked: $maker_TE): $genome.MAKER.masked\n\n";
	`cp $genome.MAKER.masked $genome.EDTA.TEanno.gff3 $genome.EDTA.TEanno.sum ../`;

	# evaluate the annotation consistency
	if ($evaluate == 1){
		# report status
		chomp ($date = `date`);
		print "$date\tEvaluate the level of inconsistency for whole-genome TE annotation (slow step):\n\n";

		# extract whole-genome TE and perform all-v-all blast, then summarize the results
		`awk '{if (\$5~/[0-9]+/ && \$1>300 && \$7-\$6>80) print \$11"\t"\$5":"\$6".."\$7}' $genome.EDTA.TEanno.out | perl $call_seq - -C $genome > $genome.EDTA.TE.fa`;
		`perl $cleanup_nested -in $genome.EDTA.TE.fa -threads $threads -minlen 80 -miniden 80 -cov 0.95 -blastplus $blastplus 2>/dev/null`;
		`for i in nested all redun; do perl $count_nested -in $genome.EDTA.TE.fa.stat -cat \$i > $genome.EDTA.TE.fa.stat.\$i.sum; done`;

		# check results and report status
		die "ERROR: TE annotation stats results not found in $genome.EDTA.TE.fa.stat!\n\n" unless -s "$genome.EDTA.TE.fa.stat";
		chomp ($date = `date`);
		print "$date\tEvaluation of TE annotation finished! Check out these files:\n
				Overall: $genome.EDTA.TE.fa.stat.all.sum
				Nested: $genome.EDTA.TE.fa.stat.nested.sum
				Non-nested: $genome.EDTA.TE.fa.stat.redun.sum\n\n";
		}

	print "\t\t\t\tIf you want to learn more about the formatting and information of these files, please visit:
	\t\t\t\thttps://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A\n\n";

	}
