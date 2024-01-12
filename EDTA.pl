#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(strftime);
use Cwd qw(abs_path);

my $version = "v2.2.0";
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
#v2.2 01/05/2024

print "
#########################################################
##### Extensive de-novo TE Annotator (EDTA) $version  #####
##### Shujun Ou (shujun.ou.1\@gmail.com)             #####
#########################################################
\n\nParameters: @ARGV\n\n\n";


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
	--rmlib	[File]	Provide the RepeatModeler library containing classified TEs to enhance
			the sensitivity especially for LINEs. If no file is provided (default),
			EDTA will generate such file for you.
	--sensitive [0|1]	Use RepeatModeler to identify remaining TEs (1) or not (0,
				default). This step may help to recover some TEs.
	--anno [0|1]	Perform (1) or not perform (0, default) whole-genome TE annotation
			after TE library construction.
	--rmout	[File]	Provide your own homology-based TE annotation instead of using the
			EDTA library for masking. File is in RepeatMasker .out format. This
			file will be merged with the structural-based TE annotation. (--anno 1
			required). Default: use the EDTA library for annotation.
	--maxdiv [0-100]	Maximum divergence (0-100%, default: 40) of repeat fragments comparing to 
				library sequences.
	--evaluate [0|1]	Evaluate (1) classification consistency of the TE annotation.
				(--anno 1 required). Default: 1.
	--exclude [File]	Exclude regions (bed format) from TE masking in the MAKER.masked
				output. Default: undef. (--anno 1 required).
	--force	[0|1]	When no confident TE candidates are found: 0, interrupt and exit
			(default); 1, use rice TEs to continue.
	--u [float]	Neutral mutation rate to calculate the age of intact LTR elements.
			Intact LTR age is found in this file: *EDTA_raw/LTR/*.pass.list.
			Default: 1.3e-8 (per bp per year, from rice).
	--repeatmodeler [path]	The directory containing RepeatModeler (default: read from ENV)
	--repeatmasker	[path]	The directory containing RepeatMasker (default: read from ENV)
	--annosine	[path]	The directory containing AnnoSINE_v2 (default: read from ENV)
	--ltrretriever	[path]	The directory containing LTR_retriever (default: read from ENV)
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
my $HQlib = ''; #curated library
my $RMlib = 'null'; #RepeatModeler library, classified
my $cds = ''; #a fasta file containing cds of this genome.
my $sensitive = 0; #0, will not run RepeatModeler to get remaining TEs (default). 1, run RepeatModeler
my $anno = 0; #0, will not annotate whole-genome TE (default). 1, annotate with RepeatMasker
my $rmout = ''; #a RM .out file for custom homology-based annotation.
my $evaluate = 1; #1 will evaluate the consistancy of the TE annotation
my $exclude = ''; #a bed file exclude from TE annotation
my $force = 0; #if there is no confident TE found in EDTA_raw, 1 will use rice TEs as raw lib, 0 will error and interrupt.
my $miu = 1.3e-8; #mutation rate, per bp per year, from rice
my $threads = 4;
my $maxdiv = 40; # maximum divergence from lib sequences for fragmented repeats
my $script_path = $FindBin::Bin;
my $EDTA_raw = "$script_path/EDTA_raw.pl";
my $EDTA_process = "$script_path/EDTA_processK.pl";
my $cleanup_proteins = "$script_path/util/cleanup_proteins.pl";
my $cleanup_TE = "$script_path/util/cleanup_TE.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/util/cleanup_nested.pl";
my $count_nested = "$script_path/util/count_nested.pl";
my $count_base = "$script_path/util/count_base.pl";
my $make_masked = "$script_path/util/make_masked.pl";
my $make_gff3 = "$script_path/util/make_gff3_with_RMout.pl";
my $protlib = "$script_path/database/alluniRefprexp082813";
my $rice_LTR = "$script_path/database/rice7.0.0.liban.LTR";
my $rice_SINE = "$script_path/database/rice7.0.0.liban.SINE";
my $rice_LINE = "$script_path/database/rice7.0.0.liban.LINE";
my $rice_TIR = "$script_path/database/rice7.0.0.liban.TIR";
my $rice_helitron = "$script_path/database/rice7.0.0.liban.Helitron";
my $rename_TE = "$script_path/util/rename_TE.pl";
#my $rename_RM = "$script_path/util/rename_RM_TE.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
my $buildSummary = "$script_path/util/buildSummary.pl"; #modified from RepeatMasker. Robert M. Hubley (rhubley@systemsbiology.org)
my $filter_gff = "$script_path/util/filter_gff3.pl";
my $combine_RMrows = "$script_path/util/combine_RMrows.pl";
my $RMout2bed = "$script_path/util/RMout2bed.pl";
my $gff2RMout = "$script_path/util/gff2RMout.pl";
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
my $annosine = "";

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
                  'rmlib=s'       	=> \$RMlib,
                  'cds=s'                => \$cds,
                  'protlib=s'            => \$protlib,
                  'sensitive=s'          => \$sensitive,
		  'anno=i'               => \$anno,
		  'rmout=s'              => \$rmout,
		  'maxdiv=i'		 => \$maxdiv,
		  'evaluate=i'           => \$evaluate,
		  'exclude=s'            => \$exclude,
		  'force=i'              => \$force,
		  'u=s'                  => \$miu,
		  'repeatmodeler=s'      => \$repeatmodeler,
		  'repeatmasker=s'       => \$repeatmasker,
		  'tesorter=s'           => \$TEsorter,
		  'blast=s'              => \$blastplus,
		  'annosine=s'		 => \$annosine,
		  'ltrretriever=s'	 => \$LTR_retriever,
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

# get $maxdiv
$maxdiv = $maxdiv*100 if $maxdiv < 1;
$maxdiv =~ s/%//g;

# check bolean
if ($maxdiv < 0 or $maxdiv > 100){die "The expected value for the div parameter is 0 - 100!\n"}
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"}
if ($sensitive != 0 and $sensitive != 1){ die "The expected value for the sensitive parameter is 0 or 1!\n"}
if ($anno != 0 and $anno != 1){ die "The expected value for the anno parameter is 0 or 1!\n"}
if ($evaluate != 0 and $evaluate != 1){ die "The expected value for the evaluate parameter is 0 or 1!\n"}
if ($force != 0 and $force != 1){ die "The expected value for the force parameter is 0 or 1!\n"}
if ($miu !~ /[0-9\.e\-]+/){ die "The expected value for the u parameter is float value without units!\n"}
if ($debug != 0 and $debug != 1){ die "The expected value for the debug parameter is 0 or 1!\n"}
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"}


# define RepeatMasker -pa parameter
#my $rm_threads = int($threads/4);
my $rm_threads = $threads;

chomp (my $date = `date`);
print "$date\tDependency checking:\n";

# check files and dependencies
die "The script EDTA_raw.pl is not found in $EDTA_raw!\n" unless -s $EDTA_raw;
die "The script EDTA_processK13.pl is not found in $EDTA_process!\n" unless -s $EDTA_process;
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
die "The rice SINE sequence library is not found in $rice_SINE!\n" unless -s $rice_SINE;
die "The rice LINE sequence library is not found in $rice_LINE!\n" unless -s $rice_LINE;
die "The rice TIR sequence library is not found in $rice_TIR!\n" unless -s $rice_TIR;
die "The rice Helitron sequence library is not found in $rice_helitron!\n" unless -s $rice_helitron;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script buildSummary.pl is not found in $buildSummary!\n" unless -s $buildSummary;
die "The script filter_gff3.pl is not found in $filter_gff!\n" unless -s $filter_gff;
die "The script RMout2bed.pl is not found in $RMout2bed!\n" unless -s $RMout2bed;
die "The script gff2RMout.pl is not found in $gff2RMout!\n" unless -s $gff2RMout;
die "The script combine_RMrows.pl is not found in $combine_RMrows!\n" unless -s $combine_RMrows;
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
# AnnoSINE
chomp ($annosine=`which AnnoSINE_v2 2>/dev/null`) if $annosine eq '';
$annosine =~ s/\s+$//;
$annosine = dirname($annosine) unless -d $annosine;
$annosine="$annosine/" if $annosine ne '' and $annosine !~ /\/$/;
die "Error: AnnoSINE is not found in the AnnoSINE path $annosine!\n" unless (-X "${annosine}AnnoSINE_v2");
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

# remove sequence annotations (content after the first space in sequence names) and replace special characters with _, convert non-ATGC bases into Ns
`perl -nle 'my \$info=(split)[0]; \$info=~s/[\\~!@#\\\$%\\^&\\*\\(\\)\\+\\\-\\=\\?\\[\\]\\{\\}\\:;",\\<\\/\\\\\|]+/_/g; \$info=~s/_+/_/g; \$info=~s/[^ATGCN]/N/gi unless /^>/; print \$info' $genome > $genome.$rand.mod`;

# try to shortern sequences
my $id_len_max = 13; # allowed longest length of a sequence ID in the input file
if ($id_len > $id_len_max){
	chomp ($date = `date`);
	print "$date\tThe longest sequence ID in the genome contains $id_len characters, which is longer than the limit ($id_len_max)\n";
	print "\t\t\t\tTrying to reformat seq IDs...\n\t\t\t\tAttempt 1...\n";
	`perl -lne 'chomp; if (s/^>+//) {s/^\\s+//; \$_=(split)[0]; s/(.{1,$id_len_max}).*/>\$1/g;} print "\$_"' $genome.$rand.mod > $genome.$rand.temp`;
	my $new_id = `grep \\> $genome.$rand.temp|sort -u|wc -l`;
	chomp ($date = `date`);
	if ($old_id == $new_id){
		$id_mode = 1;
		`mv $genome.$rand.temp $genome.mod`;
		`rm $genome.$rand.mod 2>/dev/null`;
		print "$date\tSeq ID conversion successful!\n\n";
		} else {
		print "\t\t\t\tAttempt 2...\n";
		`perl -ne 'chomp; if (/^>/) {\$_=">\$1" if /([0-9]+)/;} print "\$_\n"' $genome.$rand.mod > $genome.$rand.temp`;
		$new_id = `grep \\> $genome.$rand.temp|sort -u|wc -l`;
		if ($old_id == $new_id){
			$id_mode = 2;
			`mv $genome.$rand.temp $genome.mod`;
			`rm $genome.$rand.mod 2>/dev/null`;
			print "$date\tSeq ID conversion successful!\n\n";
			} else {
			`rm $genome.$rand.temp $genome.$rand.mod 2>/dev/null`;
			die "$date\tERROR: Fail to convert seq IDs to <= $id_len_max characters! Please provide a genome with shorter seq IDs.\n\n";
			}
		}
	} else {
	`mv $genome.$rand.mod $genome.mod`;
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

# check $RMlib
if ($RMlib ne 'null'){
	if (-s $RMlib){
		print "\tA RepeatModeler library $RMlib is provided via --rmlib. Please make sure this is a RepeatModeler2 generated and classified library (some levels of unknown classification is OK).\n\n";
		chomp ($RMlib = `realpath $RMlib`);
		`ln -s $RMlib $genome.RM2.raw.fa` unless -e "$genome.RM2.raw.fa";
		#`cp $RMlib $genome.RM2.raw.fa` unless -s "$genome.RM2.raw.fa";
		$RMlib = "$genome.RM2.raw.fa";
		} else {
		die "\tERROR: The RepeatModeler library $RMlib you specified is not found!\n\n";
		}
	}# else {
	#	`touch $genome.RM2.raw.fa 2>/dev/null`;
	#}

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


############################################################
####### Get raw LTR/SINE/LINE/TIR/Helitron candidates ######
############################################################

ALL:

# report status
chomp ($date = `date`);
print "$date\tObtain raw TE libraries using various structure-based programs: \n";

# Get raw TE candidates
`perl $EDTA_raw --genome $genome --overwrite $overwrite --species $species --u $miu --threads $threads --genometools $genometools --ltrretriever $LTR_retriever --blastplus $blastplus --tesorter $TEsorter --GRF $GRF --trf_path $trf --repeatmasker $repeatmasker --repeatmodeler $repeatmodeler --annosine $annosine --convert_seq_name 0 --rmlib $RMlib`;

chdir "$genome.EDTA.raw";

# Force to use rice TEs when raw.fa is empty
if ($force eq 1){
	`cp $rice_LTR $genome.LTR.raw.fa` unless -s "$genome.LTR.raw.fa";
	`cp $rice_LTR $genome.LTR.intact.raw.fa` unless -s "$genome.LTR.intact.raw.fa";
	`cp $rice_LINE $genome.LINE.raw.fa` unless -s "$genome.LINE.raw.fa";
	`cp $rice_SINE $genome.SINE.raw.fa` unless -s "$genome.SINE.raw.fa";
	`cp $rice_TIR $genome.TIR.intact.raw.fa` unless -s "$genome.TIR.intact.raw.fa";
	`cp $rice_helitron $genome.Helitron.intact.raw.fa` unless -s "$genome.Helitron.intact.raw.fa";
	}

# check results and report status
die "ERROR: Raw LTR results not found in $genome.EDTA.raw/$genome.LTR.raw.fa and $genome.EDTA.raw/$genome.LTR.intact.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of intact LTRs in your genome. Consider to use the --force 1 parameter to overwrite this check\n" unless (-s "$genome.LTR.raw.fa" and -s "$genome.LTR.intact.raw.fa");
die "ERROR: Raw SINE results not found in $genome.EDTA.raw/$genome.SINE.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of SINEs in your genome.\n" unless -e "$genome.SINE.raw.fa"; # allow empty file
die "ERROR: Raw LINE results not found in $genome.EDTA.raw/$genome.LINE.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of LINEs in your genome.\n" unless -e "$genome.LINE.raw.fa"; # allow empty file
die "ERROR: Raw TIR results not found in $genome.EDTA.raw/$genome.TIR.intact.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of intact TIRs in your genome. Consider to use the --force 1 parameter to overwrite this check\n" unless -s "$genome.TIR.intact.raw.fa";
die "ERROR: Raw Helitron results not found in $genome.EDTA.raw/$genome.Helitron.intact.raw.fa\n\tIf you believe the program is working properly, this may be caused by the lack of intact Helitrons in your genome. Consider to use the --force 1 parameter to overwrite this check\n" unless -s "$genome.Helitron.intact.raw.fa";

# combine intact TEs
`cat $genome.LTR.intact.raw.fa $genome.TIR.intact.raw.fa $genome.Helitron.intact.raw.fa > $genome.EDTA.intact.raw.fa`;
`cat $genome.TIR.intact.raw.bed $genome.Helitron.intact.raw.bed | perl $bed2gff - TE_struc > $genome.EDTA.intact.gff3.temp`;
`cat $genome.LTR.intact.raw.gff3 >> $genome.EDTA.intact.gff3.temp`;
`sort -sV -k1,1 -k4,4 $genome.EDTA.intact.gff3.temp | grep -v '^#' > $genome.EDTA.intact.raw.gff3; rm $genome.EDTA.intact.gff3.temp`;

chomp ($date = `date`);
print "$date\tObtain raw TE libraries finished.
\t\t\t\tAll intact TEs found by EDTA: \n\t\t\t\t\t$genome.EDTA.intact.raw.fa \n\t\t\t\t\t$genome.EDTA.intact.raw.gff3\n\n";
chdir "..";


############################################################
####### Filter LTR/SINE/LINE/TIR/Helitron candidates #######
############################################################

FILTER:

# report status
chomp ($date = `date`);
print "$date\tPerform EDTA advance filtering for raw TE candidates and generate the stage 1 library: \n\n";

# remove existing results
`rm ./$genome.EDTA.combine/* 2>/dev/null` if $overwrite == 1;

# Filter raw TE candidates and the make stage 1 library
`perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -ltrint $genome.EDTA.raw/$genome.LTR.intact.raw.fa -line $genome.EDTA.raw/$genome.LINE.raw.fa -sine $genome.EDTA.raw/$genome.SINE.raw.fa -tir $genome.EDTA.raw/$genome.TIR.intact.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.intact.raw.fa -repeatmasker $repeatmasker -blast $blastplus -threads $threads`;

# check results, remove intermediate files, and report status
die "ERROR: Stage 1 library not found in $genome.EDTA.combine/$genome.EDTA.fa.stg1" unless -s "$genome.EDTA.combine/$genome.EDTA.fa.stg1";
chdir "$genome.EDTA.combine";
`rm ./$genome.LTR.raw.fa*Q* ./$genome.LTR.intact.raw.fa*Q* ./$genome.TIR.intact.raw.fa*Q* ./$genome.Helitron.intact.raw.fa*Q* ./$genome.TIR.Helitron.fa*Q* $genome*tbl $genome*out $genome*cleanup $genome*RMoutput $genome*stg1.raw* $genome.LTR.raw.fa-* $genome.LTR.intact.raw.fa-* $genome.TIR.intact.raw.fa-* $genome.Helitron.intact.raw.fa-* $genome.LINE_LTR.raw.fa $genome.LTR.SINE.LINE.fa *.ndb *.not *.ntf *.nto *.cat.gz *.cat *.masked *.ori.out *.nhr *.nin *.nsq 2>/dev/null` unless $debug eq 1;

chdir "..";
chomp ($date = `date`);
print "$date\tEDTA advance filtering finished.\n\n";


####################################
###### Final purge CDS in TEs ######
####################################

FINAL:

# report status
chomp ($date = `date`);
print "$date\tPerform EDTA final steps to generate a non-redundant comprehensive TE library.\n\n";

# Make the final working directory
`mkdir $genome.EDTA.final` unless -e "$genome.EDTA.final" && -d "$genome.EDTA.final";
chdir "$genome.EDTA.final";
`rm ./* 2>/dev/null` if $overwrite == 1;
`cp ../$genome.EDTA.raw/$genome.RM2.fa ./`;
`cp ../$genome.EDTA.combine/$genome.EDTA.fa.stg1 ./`;
`cp ../$cds ./` if $cds ne '';
`cp ../$HQlib ./` if $HQlib ne '';
`cp ../$genome.EDTA.combine/$genome.EDTA.intact.fa.cln ./$genome.EDTA.intact.fa.cln`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.raw.fa ./`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.raw.gff3 ./`;
`cp ../$exclude ./` if $exclude ne '';

# identify remaining TEs in the filtered RM2 library
if ($sensitive == 1 and -s "$genome.RM2.fa"){
	print "\t\t\t\tIdentify TE families from RepeatModeler results that are missed by structure-based methods.\n\n";
	chomp ($date = `date`);
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.EDTA.fa.stg1 $genome.RM2.fa 2>/dev/null`;
	`cp $genome.RM2.fa $genome.RM2.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	# clean up tandem and coding sequences in the RM2 library
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.RM2.fa.masked > $genome.RM2.fa.stg1`;
	`perl $cleanup_proteins -seq $genome.RM2.fa.stg1 -rmdnate 0 -rmline 0 -rmprot 1 -protlib $protlib -blast $blastplus -threads $threads`;
	if (-s "$genome.RM2.fa.stg1.clean"){
		`cat $genome.EDTA.fa.stg1 $genome.RM2.fa.stg1.clean > $genome.EDTA.raw.fa`;
		} else {
		print "\t\t\t\tNo extra repeat sequences found in the RepeatModeler output.\n\n";
		`cp $genome.EDTA.fa.stg1 $genome.EDTA.raw.fa`;
		}
	} else {
	print "\t\t\t\tSkipping the RepeatModeler results (--sensitive 0).\n\t\t\t\tRun EDTA.pl --step final --sensitive 1 if you want to add RepeatModeler results.\n\n";
	`cp $genome.EDTA.fa.stg1 $genome.EDTA.raw.fa`;
	}

# remove CDS in the non-redundant library and intact TEs
if (-s "$cds"){
	# report status
	chomp ($date = `date`);

	# cleanup TE-related sequences in the CDS file with TEsorter
	print "$date\tClean up TE-related sequences in the CDS file with TEsorter.\n\n";
	`perl $cleanup_TE -cds $cds -minlen 300 -tesorter $TEsorter -repeatmasker $repeatmasker -t $threads -rawlib $genome.EDTA.raw.fa 2>/dev/null`;
	`rm ./$cds ./$cds.code.r* 2>/dev/null` unless $debug eq 1;
	die "\tERROR: The $cds file is empty after TE clean up. Please check the file and $cds.code.noTE.\n\n" unless -s "$cds.code.noTE";
	$cds = "$cds.code.noTE";

	# remove cds-related sequences in the EDTA library
	print "\t\t\t\tRemove CDS-related sequences in the EDTA library.\n\n";
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.raw.fa 2>/dev/null`;
	`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.3 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.EDTA.raw.fa.masked > $genome.EDTA.raw.fa.cln`;

	# remove cds-related sequences in intact TEs
	print "\t\t\t\tRemove CDS-related sequences in intact TEs.\n\n";
	$rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.intact.fa.cln 2>/dev/null`;
	`cp $genome.EDTA.intact.fa.cln $genome.EDTA.intact.fa.cln.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.8 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 0 -f $genome.EDTA.intact.fa.cln.masked > $genome.EDTA.intact.fa.cln.rmCDS`;
	`perl $output_by_list 1 $genome.EDTA.intact.fa.cln 1 $genome.EDTA.intact.fa.cln.masked.cleanup -ex -FA > $genome.EDTA.intact.fa.cln2`;
	} else {
	print "\t\t\t\tSkipping the CDS cleaning step (--cds [File]) since no CDS file is provided or it's empty.\n\n";
	copy_file("$genome.EDTA.raw.fa", "./$genome.EDTA.raw.fa.cln");
	copy_file("$genome.EDTA.intact.fa.cln", "./$genome.EDTA.intact.fa.cln2");
	}

# Final rounds of redundancy removal and make final EDTA library
`perl $cleanup_nested -in $genome.EDTA.raw.fa.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blastplus 2>/dev/null`;

# rename all TEs in the EDTA library
`perl $rename_TE $genome.EDTA.raw.fa.cln.cln > $genome.EDTA.TElib.fa`;
#`perl $rename_TE $genome.EDTA.raw.fa.cln.cln | perl $format_TElib - > $genome.EDTA.TElib.fa`;

# identify novel TEs using the user provided $HQlib
if ($HQlib ne ''){
	# report status
	chomp ($date = `date`);
	print "$date\tCombine the high-quality TE library $HQlib with the EDTA library:\n\n";

	# remove known TEs in the EDTA library
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $HQlib $genome.EDTA.TElib.fa 2>/dev/null`;
	`cp $genome.EDTA.TElib.fa $genome.EDTA.TElib.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 0 -f $genome.EDTA.TElib.fa.masked > $genome.EDTA.TElib.novel.fa`;
	rename "$genome.EDTA.TElib.fa", "$genome.EDTA.TElib.ori.fa";
	`cat $HQlib $genome.EDTA.TElib.novel.fa > $genome.EDTA.TElib.fa`;
	copy_file("$genome.EDTA.TElib.novel.fa", "..");
	}

# reclassify intact TEs with the TE lib #113
`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $genome.EDTA.TElib.fa $genome.EDTA.intact.fa.cln2 2>/dev/null` unless (-s "$genome.EDTA.intact.fa.cln2.out" and $overwrite == 0);
die "ERROR: The masked file for $genome.EDTA.intact.fa.cln2 is not found! The RepeatMasker annotation on this file may be failed. Please check the $genome.EDTA.TElib.fa file for sequence naming formats especially when you provide a library via --curatedlib.\n" unless -s "$genome.EDTA.intact.fa.cln2.out";
`perl $reclassify -seq $genome.EDTA.intact.fa.cln2 -RM $genome.EDTA.intact.fa.cln2.out -cov 80 -len 80 -iden 60`; # 80-80-60

# remove inconsistently classified intact TEs and generate the final intact TEs
`perl $output_by_list 1 $genome.EDTA.intact.fa.cln2.rename 1 $genome.EDTA.intact.fa.cln2.false.list -ex -FA > $genome.EDTA.intact.fa`;

## generate clean intact gff3
# update the family names in the intact.raw.gff3 file
`perl $rename_by_list $genome.EDTA.intact.raw.gff3 $genome.EDTA.intact.fa.cln2.rename.list 1 > $genome.EDTA.intact.raw.gff3.rename`;
`sed 's/.*Name=//; s/;Classifica.*//' $genome.EDTA.intact.raw.gff3.rename | sort -u > $genome.EDTA.intact.raw.gff3.rename.famlist`;
# get a dirty list of intact.gff
`grep \\> $genome.EDTA.intact.fa | sed 's/>//; s/#.*//' | perl $output_by_list 1 $genome.EDTA.intact.raw.gff3.rename.famlist 1 - -ex | awk '{print "Name\\t"\$1"\\nParent\\t"\$1"\\nID\\t"\$1}' > $genome.EDTA.intact.raw.gff3.rename.dirtlist`;
# first attempt purging the gff3
`perl $filter_gff $genome.EDTA.intact.raw.gff3.rename $genome.EDTA.intact.raw.gff3.rename.dirtlist > $genome.EDTA.intact.gff3`;
# remake the remove list an purge again
`perl -nle 'my \$id = \$1 if /=(repeat_region[0-9]+);/; print "Parent\\t\$id\nName\\t\$id" if defined \$id' $genome.EDTA.intact.raw.gff3.rename.removed >> $genome.EDTA.intact.raw.gff3.rename.dirtlist`;
`perl $filter_gff $genome.EDTA.intact.raw.gff3.rename $genome.EDTA.intact.raw.gff3.rename.dirtlist > $genome.EDTA.intact.gff3`;

# check results
die "ERROR: Final TE library not found in $genome.EDTA.TElib.fa" unless -s "$genome.EDTA.TElib.fa";
die "ERROR: Intact TE annotation not found in $genome.EDTA.intact.gff3" unless -s "$genome.EDTA.intact.gff3";
copy_file("$genome.EDTA.TElib.fa", "..");
copy_file("$genome.EDTA.intact.fa", "..");
copy_file("$genome.EDTA.intact.gff3", "..");

# remove intermediate files
`rm $genome.EDTA.intact.fa.cln.* $genome.EDTA.raw.fa.* $genome.EDTA.TElib.fa.* $genome.LTR.TIR.Helitron.fa.stg1.* $genome.masked *.cat.gz 2>/dev/null` if $debug eq 0;

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
	`rm $genome.EDTA.TElib.fa* 2>/dev/null`; # clean up libraries
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
		`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div $maxdiv -lib $genome.EDTA.TElib.fa $genome 2>/dev/null`;
		}
	die "ERROR: RepeatMasker results not found in $genome.out!\n\n" unless -s "$genome.out" or -s "$genome.mod.out";

	# exclude regions from TE annotation and make whole-genome TE annotation
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv $maxdiv -minscore 300 -minlen 80 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
	# combine RepeatMasker lines that appears to be the same element
	`perl $combine_RMrows -rmout $genome.out.new -maxdiv 3.5 -maxgap 35`;
	`mv $genome.out.new.cmb $genome.EDTA.RM.out`;
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
	`perl $gff2RMout $genome.EDTA.TEanno.split.gff3 $genome.EDTA.TEanno.split.out`;

	# make summary table for the non-overlapping annotation
	`perl $count_base $genome > $genome.stats`;
	`perl -nle 'my (\$chr, \$s, \$e, \$anno, \$dir, \$supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 \$chr \$s \$e NA \$dir \$anno \$supfam"' $genome.EDTA.TEanno.split.bed > $genome.EDTA.TEanno.out`;
	`perl $buildSummary -maxDiv 40 -stats $genome.stats $genome.EDTA.TEanno.out > $genome.EDTA.TEanno.sum 2>/dev/null`;
	my $tot_TE = `grep Total $genome.EDTA.TEanno.sum|grep %|awk '{print \$4}'`;
	chomp $tot_TE;

	# make low-threshold masked genome for MAKER
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 1 -misschar N -threads $threads -exclude $exclude` unless (-s "$genome.MAKER.masked" and $overwrite == 0);
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
	`cp $genome.MAKER.masked ../`; # make no backup for this file
	copy_file("$genome.EDTA.TEanno.gff3", "..");
	copy_file("$genome.EDTA.TEanno.sum", "..");

	# evaluate the annotation consistency
	if ($evaluate == 1){
		# report status
		chomp ($date = `date`);
		print "$date\tEvaluate the level of inconsistency for whole-genome TE annotation:\n\n";

		# extract whole-genome TE and perform all-v-all blast, then summarize the results
		`awk '{if (\$5~/[0-9]+/ && \$1>300 && \$7-\$6>80) print \$11"\t"\$5":"\$6".."\$7}' $genome.EDTA.TEanno.out | perl $call_seq - -C $genome > $genome.EDTA.TE.fa`;
		`perl $cleanup_nested -in $genome.EDTA.TE.fa -threads $threads -minlen 80 -miniden 80 -cov 0.95 -blastplus $blastplus -iter 1 -maxcount 100000 2>/dev/null`;
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


##########################
###### Subroutines #######
##########################

sub copy_file {
	my ($file, $path) = ($_[0], $_[1]);
	# Generate new name with current date and time
	my $new_name = $file . "_" . strftime("%Y%m%d_%H%M%S", localtime);
	
	# resolve symlinks and existing files
	if (-l "$path/$file") {
		# File is a symbolic link
		unlink "$path/$file" or die "ERROR: Failed to remove symbolic link for $path/$file\n\n";
		} elsif (-f "$path/$file") {
        	# File is a regular file
        	rename "$path/$file", "$path/$new_name" or die "ERROR: Failed to rename file: $path/$file\n\n";
        	}

        # copy file to the path if it's not the current path
	`cp $file $path` if abs_path($path) ne abs_path('.');
	}
