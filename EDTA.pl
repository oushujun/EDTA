#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(strftime);
use Cwd qw(abs_path);
use File::Path qw(rmtree);

my $version = "v2.3.3";
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
#v2.2.2 08/08/2024
#v2.3 03/09/2026
#v2.3.1 05/10/2026
#v2.3.2 07/10/2026
#v2.3.3 08/31/2026

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
	--modules [all|plant|list]	Which raw TE discovery modules to run. Default:
					all (ltr, sine, line, tir, helitron). plant is a
					shortcut for ltr,tir,helitron — it skips SINE and
					LINE, which annotate <2% of most plant genomes but
					cost the most time (AnnoSINE + RepeatModeler). You
					may also give an explicit comma list, e.g.
					ltr,tir,helitron,line. Excluded modules leave empty
					library files; downstream stages run unchanged.
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
	--maker	[0|1]	Produce (1) or not (0, default) the low-threshold MAKER.masked
				genome for MAKER gene annotation. (--anno 1 required).
	--exclude [File]	Exclude regions (bed format) from TE masking in the MAKER.masked
				output. Default: undef. (--anno 1 and --maker 1 required).
	--force	[0|1]	(default: 0) 0: When no confident TE candidates are found, interrupt and exit.
			             1: Use rice TEs to continue.
	--wholeelement [0|1]	Keep LTR retrotransposons as whole elements in the library instead
				of splitting into LTR and INT regions (default: 0).
	--u [float]	Neutral mutation rate to calculate the age of intact LTR elements.
			Intact LTR age is found in this file: *EDTA_raw/LTR/*.pass.list.
			Default: 1.3e-8 (per bp per year, from rice).
	--repeatmodeler [path]	The directory containing RepeatModeler (default: read from ENV)
	--repeatmasker	[path]	The directory containing RepeatMasker (default: read from ENV)
	--annosine	[path]	The directory containing AnnoSINE_v2 (default: read from ENV)
	--tirlearner	[path]	The directory containing TIR-Learner (default: read from ENV)
	--ltrretriever	[path]	The directory containing LTR_retriever (default: read from ENV)
	--check_dependencies Check if dependencies are fullfiled and quit
	--threads|-t [int]	Number of theads to run this script (default: 4)
	--tmpdir [Dir]		Directory for the temporary files of this run and all its
				child tools (default: .EDTA.tmp.<pid> in the working
				directory). Set EDTA_TMPDIR_KEEP=1 to keep the
				inherited TMPDIR instead.
	--debug	 [0|1]	Retain intermediate files (default: 0)
	--help|-h 	Display this help info
\n";

# pre-defined
my $genome = '';
my $check_dependencies = undef;
my $species = "others";
my $modules = "all"; #which raw TE modules to run: all, plant (= ltr,tir,helitron), or a comma list
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
my $maker = 0; #0, will not produce the low-threshold MAKER.masked genome (default). 1, produce it.
my $force = 0; #if there is no confident TE found in EDTA_raw, 1 will use rice TEs as raw lib, 0 will error and interrupt.
my $miu = 1.3e-8; #mutation rate, per bp per year, from rice
my $threads = 4;
my $tmpdir = ''; #private scratch dir for TMPDIR isolation; default: .EDTA.tmp.$$
my $maxdiv = 40; # maximum divergence from lib sequences for fragmented repeats
my $script_path = $FindBin::Bin;
my $EDTA_raw = "$script_path/EDTA_raw.pl";
my $EDTA_process = "$script_path/EDTA_processK.pl";
my $cleanup_proteins = "$script_path/bin/cleanup_proteins.pl";
my $cleanup_TE = "$script_path/bin/cleanup_TE.pl";
my $cleanup_tandem = "$script_path/bin/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/bin/cleanup_nested.pl";
my $count_nested = "$script_path/bin/count_nested.pl";
my $evaluation = "$script_path/bin/evaluation.pl";
my $count_base = "$script_path/bin/count_base.pl";
my $make_masked = "$script_path/bin/make_masked.pl";
my $make_gff3 = "$script_path/bin/make_gff3_with_RMout.pl";
my $protlib = "$script_path/database/alluniRefprexp082813";
my $rice_LTR = "$script_path/database/rice7.0.0.liban.LTR";
my $rice_SINE = "$script_path/database/rice7.0.0.liban.SINE";
my $rice_LINE = "$script_path/database/rice7.0.0.liban.LINE";
my $rice_TIR = "$script_path/database/rice7.0.0.liban.TIR";
my $rice_helitron = "$script_path/database/rice7.0.0.liban.Helitron";
my $rename_TE = "$script_path/bin/rename_TE.pl";
my $update_LTRbound = "$script_path/bin/update_LTRbound.pl";
my $ltrbound_denovo = "$script_path/bin/ltrbound_denovo.py";
my $seqid_codec = "$script_path/bin/seqid_codec.pl";
#my $rename_RM = "$script_path/bin/rename_RM_TE.pl";
my $call_seq = "$script_path/bin/call_seq_by_list.pl";
my $buildSummary = "$script_path/bin/buildSummary.pl"; #modified from RepeatMasker. Robert M. Hubley (rhubley@systemsbiology.org)
my $filter_gff = "$script_path/bin/filter_gff3.pl";
my $combine_RMrows = "$script_path/bin/combine_RMrows.pl";
my $RMout2bed = "$script_path/bin/RMout2bed.pl";
my $gff2RMout = "$script_path/bin/gff2RMout.pl";
my $bed2gff = "$script_path/bin/bed2gff.pl";
my $gff2bed = "$script_path/bin/gff2bed.pl";
my $gff2gtf = "$script_path/bin/gff2gtf.pl";
my $get_frag = "$script_path/bin/get_frag.pl";
my $keep_nest = "$script_path/bin/keep_nest.pl";
my $combine_overlap = "$script_path/bin/combine_overlap.pl";
my $split_overlap = "$script_path/bin/split_overlap.pl";
my $reclassify = "$script_path/bin/classify_by_lib_RM.pl";
my $rename_by_list = "$script_path/bin/rename_by_list.pl";
my $output_by_list = "$script_path/bin/output_by_list.pl";
my $format_TElib = "$script_path/bin/format_TElib.pl";
my $format_gff3 = "$script_path/bin/format_gff3.pl";
my $add_id = "$script_path/bin/add_id.pl";
my $div_table = "$script_path/bin/div_table2.pl";
my $div_plot = "$script_path/bin/div_plot2.R";
my $density_table = "$script_path/bin/density_table.py";
my $density_plot = "$script_path/bin/density_plot.R";
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
my $TIR_Learner = "";

my $wholeelement = 0; #0, split LTR library into LTR/INT (default); 1, keep whole elements
my $beta2 = 0; #0, beta2 is not ready. 1, developer mode.
#my $reanno = 0; #0, use existing whole-genome RM results (beta); 1, de novo Repeatmasker using the EDTA library (default)
my $debug = 0;
my $help = undef;

# read parameters
if ( !GetOptions( 'genome=s'            => \$genome,
                  'species=s'           => \$species,
                  'modules=s'           => \$modules,
                  'step=s'              => \$step,
                  'overwrite=i'         => \$overwrite,
                  'curatedlib=s'        => \$HQlib,
                  'rmlib=s'       	=> \$RMlib,
                  'cds=s'                => \$cds,
                  'protlib=s'            => \$protlib,
                  'sensitive=i'          => \$sensitive,
		  'anno=i'               => \$anno,
		  'rmout=s'              => \$rmout,
		  'maxdiv=f'		 => \$maxdiv,
		  'evaluate=i'           => \$evaluate,
		  'exclude=s'            => \$exclude,
		  'maker=i'              => \$maker,
		  'force=i'              => \$force,
		  'u=s'                  => \$miu,
		  'repeatmodeler=s'      => \$repeatmodeler,
		  'repeatmasker=s'       => \$repeatmasker,
		  'tesorter=s'           => \$TEsorter,
		  'blast=s'              => \$blastplus,
		  'annosine=s'		 => \$annosine,
		  'tirlearner=s'	 => \$TIR_Learner,
		  'ltrretriever=s'	 => \$LTR_retriever,
		  'threads|t=i'          => \$threads,
		  'tmpdir=s'             => \$tmpdir,
		  'wholeelement=i'       => \$wholeelement,
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
$maxdiv = int($maxdiv*10+0.5)/10;
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"}
if ($sensitive != 0 and $sensitive != 1){ die "The expected value for the sensitive parameter is 0 or 1!\n"}
if ($anno != 0 and $anno != 1){ die "The expected value for the anno parameter is 0 or 1!\n"}
if ($evaluate != 0 and $evaluate != 1){ die "The expected value for the evaluate parameter is 0 or 1!\n"}
if ($force != 0 and $force != 1){ die "The expected value for the force parameter is 0 or 1!\n"}
if ($maker != 0 and $maker != 1){ die "The expected value for the maker parameter is 0 or 1!\n"}
if ($miu !~ /^[0-9.eE+-]+$/ or $miu !~ /^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$/){ die "The expected value for the u parameter is float value without units!\n"}
if ($debug != 0 and $debug != 1){ die "The expected value for the debug parameter is 0 or 1!\n"}
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"}
if ($threads < 1){ die "The expected value for the threads parameter is an integer >= 1!\n"}

# --- TMPDIR isolation: keep descendants off the system /tmp ----------------
# Unless explicitly kept, point TMPDIR at a private scratch dir in the
# working directory so no tool (python tempfile, sort spill, blast temp)
# can fill the machine's /tmp. --tmpdir selects a custom location
# (e.g. a node-local SSD); EDTA_TMPDIR_KEEP=1 keeps the environment as-is.
my $edta_own_tmp = 0;
unless (defined $ENV{EDTA_TMPDIR_KEEP} and $ENV{EDTA_TMPDIR_KEEP} eq '1'){
	if ($tmpdir ne '' and -d $tmpdir){
		$ENV{TMPDIR} = $tmpdir; # user-selected dir, never removed by EDTA
		} else {
		$ENV{TMPDIR} = abs_path(".")."/.EDTA.tmp.$$"; # run-private default scratch
		$edta_own_tmp = 1;
		}
	mkdir($ENV{TMPDIR}) unless -d $ENV{TMPDIR};
	}


# define RepeatMasker -pa parameter
#my $rm_threads = int($threads/4);
my $rm_threads = $threads;

chomp (my $date = `date`);
print "$date\tDependency checking:\n";

# check files and dependencies
die "The script EDTA_raw.pl is not found in $EDTA_raw!\n" unless -s $EDTA_raw;
die "The script EDTA_processK.pl is not found in $EDTA_process!\n" unless -s $EDTA_process;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;
die "The script cleanup_TE.pl is not found in $cleanup_TE!\n" unless -s $cleanup_TE;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script count_nested.pl is not found in $count_nested!\n" unless -s $count_nested;
die "The script evaluation.pl is not found in $evaluation!\n" unless -s $evaluation;
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
die "The script gff2gtf.pl is not found in $gff2gtf!\n" unless -s $gff2gtf;
die "The script get_frag.pl is not found in $get_frag!\n" unless -s $get_frag;
die "The script keep_nest.pl is not found in $keep_nest!\n" unless -s $keep_nest;
die "The script combine_overlap.pl is not found in $combine_overlap!\n" unless -s $combine_overlap;
die "The script split_overlap.pl is not found in $split_overlap!\n" unless -s $split_overlap;
die "The script classify_by_lib_RM.pl is not found in $reclassify!\n" unless -s $reclassify;
die "The script rename_by_list.pl is not found in $rename_by_list!\n" unless -s $rename_by_list;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script format_gff3.pl is not found in $format_gff3!\n" unless -s $format_gff3;
die "The script add_id.pl is not found in $add_id!\n" unless -s $add_id;
die "The script div_table2.pl is not found in $div_table!\n" unless -s $div_table;
die "The script div_plot2.R is not found in $div_plot!\n" unless -s $div_plot;
die "The script density_table.py is not found in $density_table!\n" unless -s $density_table;
die "The script density_plot.R is not found in $density_plot!\n" unless -s $density_plot;

# GenomeTools
chomp ($genometools=`command -v gt 2>/dev/null`) if $genometools eq '';
$genometools =~ s/\s+$//;
$genometools = dirname($genometools) unless -d $genometools;
$genometools="$genometools/" if $genometools ne '' and $genometools !~ /\/$/;
die "Error: gt is not found in the genometools path $genometools!\n" unless -X "${genometools}gt";
# AnnoSINE
chomp ($annosine=`command -v AnnoSINE_v2 2>/dev/null`) if $annosine eq '';
$annosine =~ s/\s+$//;
$annosine = dirname($annosine) unless -d $annosine;
$annosine="$annosine/" if $annosine ne '' and $annosine !~ /\/$/;
die "Error: AnnoSINE is not found in the AnnoSINE path $annosine!\n" unless (-X "${annosine}AnnoSINE_v2");
# LTR_retriever
chomp ($LTR_retriever=`command -v LTR_retriever 2>/dev/null`) if $LTR_retriever eq '';
$LTR_retriever =~ s/\s+$//;
$LTR_retriever = dirname($LTR_retriever) unless -d $LTR_retriever;
$LTR_retriever="$LTR_retriever/" if $LTR_retriever ne '' and $LTR_retriever !~ /\/$/;
die "Error: LTR_retriever is not found in the LTR_retriever path $LTR_retriever!\n" unless -X "${LTR_retriever}LTR_retriever";
# RepeatMasker
my $rand=int(rand(1000000));
chomp ($repeatmasker=`command -v RepeatMasker 2>/dev/null`) if $repeatmasker eq '';
$repeatmasker =~ s/\s+$//;
$repeatmasker = dirname($repeatmasker) unless -d $repeatmasker;
$repeatmasker="$repeatmasker/" if $repeatmasker ne '' and $repeatmasker !~ /\/$/;
die "Error: RepeatMasker is not found in the RepeatMasker path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";
`cp $script_path/database/dummy060817.fa ./dummy060817.fa.$rand`;
my $RM_test=`${repeatmasker}RepeatMasker -e ncbi -q -pa 1 -no_is -nolow dummy060817.fa.$rand -lib dummy060817.fa.$rand 2>/dev/null`;
`rm dummy060817.fa.$rand* 2>/dev/null`;
die "Error: The RMblast engine is not installed in RepeatMasker!\n" unless $RM_test=~s/done//gi;
# RepeatModeler
chomp ($repeatmodeler=`command -v RepeatModeler 2>/dev/null`) if $repeatmodeler eq '';
$repeatmodeler =~ s/\s+$//;
$repeatmodeler = dirname($repeatmodeler) unless -d $repeatmodeler;
$repeatmodeler="$repeatmodeler/" if $repeatmodeler ne '' and $repeatmodeler !~ /\/$/;
die "Error: RepeatModeler is not found in the RepeatModeler path $repeatmodeler!\n" unless -X "${repeatmodeler}RepeatModeler";
# makeblastdb, blastn, blastx
chomp ($blastplus=`command -v makeblastdb 2>/dev/null`) if $blastplus eq '';
$blastplus =~ s/\s+$//;
$blastplus = dirname($blastplus) unless -d $blastplus;
$blastplus="$blastplus/" if $blastplus ne '' and $blastplus !~ /\/$/;
die "Error: makeblastdb is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}makeblastdb";
die "Error: blastn is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";
die "Error: blastx is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastx";
# TEsorter
chomp ($TEsorter=`command -v TEsorter 2>/dev/null`) if $TEsorter eq '';
$TEsorter =~ s/\s+$//;
$TEsorter = dirname($TEsorter) unless -d $TEsorter;
$TEsorter="$TEsorter/" if $TEsorter ne '' and $TEsorter !~ /\/$/;
die "Error: TEsorter is not found in the TEsorter path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
# mdust
chomp ($mdust=`command -v mdust 2>/dev/null`) if $mdust eq '';
$mdust =~ s/\s+$//;
$mdust = dirname($mdust) unless -d $mdust;
$mdust = "$mdust/" if $mdust ne '' and $mdust !~ /\/$/;
die "Error: mdust is not found in the mdust path $mdust!\n" unless -X "${mdust}mdust";
# trf
chomp ($trf=`command -v trf 2>/dev/null`) if $trf eq '';
$trf=~s/\n$//;
die "Error: Tandem Repeat Finder is not found in the TRF path $trf!\n" unless $trf ne '' and -X $trf;
# GRF (only needed by the TIR module; EDTA_raw skips TIR with a warning when absent)
chomp ($GRF = `command -v grf-main 2>/dev/null`) if $GRF eq '';
$GRF =~ s/\n$//;
if ($GRF eq '' or !-X $GRF){
	print STDERR "Warning: The Generic Repeat Finder (GRF) is not found in the GRF path: $GRF!\n\t\tThe TIR module will be skipped.\n\n";
	}

print "\tAll passed!\n\n";
exit if $check_dependencies;

# make a softlink to the user-provided files
my $genome_file = basename($genome);
softlink_file($genome, $genome_file);
$genome = $genome_file;

# check if duplicated sequences found (single pass for total and unique ID counts)
my $id_counts = `grep -a \\> $genome|sort|uniq -c|awk '{t+=\$1} END{print t+0" "NR}'`;
my ($raw_id, $old_id) = $id_counts =~ /(\d+)\s+(\d+)/;
if ($raw_id > $old_id){
	chomp ($date = `date`);
	die "$date\tERROR: Identical sequence IDs found in the provided genome! Please resolve this issue and try again.\n";
	}

# Normalize genome: clean IDs, replace special chars, convert non-ATGC to N,
# and encode sequence IDs with short base-62 codes if they are too long for
# the rmblastn 50-char ID limit (accounting for LTR_retriever coordinate appending).
my $seqid_mapfile = "";
if (-s "$genome.mod" and $overwrite == 0){
	# Resume: use existing normalized genome
	$genome = "$genome.mod";
	if (-s "$genome.seqid.map"){
		# Verify integrity: seq count in .mod must match lines in mapfile
		my $mod_count = `grep -ac \\> $genome`;
		chomp $mod_count;
		my $map_count = `wc -l < $genome.seqid.map`;
		chomp $map_count;
		if ($mod_count != $map_count){
			chomp ($date = `date`);
			die "$date\tERROR: Integrity check failed for $genome.seqid.map " .
				"($mod_count sequences vs $map_count map entries). " .
				"Please rerun with --overwrite 1.\n";
			}
		$seqid_mapfile = abs_path("$genome.seqid.map");
		chomp ($date = `date`);
		print "$date\tExisting encoded genome $genome and mapping file found, integrity verified.\n\n";
		} else {
		chomp ($date = `date`);
		print "$date\tExisting normalized genome $genome found (no encoding needed).\n\n";
		}
	} else {
	# Fresh run: clean and encode genome
	chomp ($date = `date`);
	print "$date\tCleaning and normalizing sequence IDs...\n";
	`perl $seqid_codec encode_fasta $genome $genome.mod.seqid.map $genome.mod`;
	if ($? != 0){
		die "$date\tERROR: Genome normalization failed. Check error messages above.\n";
		}
	$genome = "$genome.mod";
	if (-s "$genome.seqid.map"){
		$seqid_mapfile = abs_path("$genome.seqid.map");
		print "\tSequence IDs encoded. Mapping file: $genome.seqid.map\n\n";
		} else {
		print "\tSequence IDs are short enough, no encoding needed.\n\n";
		}
	# Verify unique ID count
	my $new_id = `grep -a \\> $genome|sort -u|wc -l`;
	chomp $new_id;
	if ($old_id != $new_id){
		chomp ($date = `date`);
		die "$date\tERROR: Seq ID normalization produced non-unique IDs. Please check your genome file.\n";
		}
	}

# check $HQlib
if ($HQlib ne ''){
	if (-s $HQlib){
		print "\tA custom library $HQlib is provided via --curatedlib. Please make sure this is a manually curated library but not machine generated.\n\n";
		chomp ($HQlib = `realpath $HQlib`);
		my $HQlib_file = basename($HQlib);
		softlink_file($HQlib, $HQlib_file);
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
		softlink_file($RMlib, "$genome.RM2.raw.fa");
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
		softlink_file($cds, $cds_file);
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
		softlink_file($exclude, $exclude_file);
		$exclude = $exclude_file;
		} else {
		die "\tERROR: The exclusion BED file $exclude you specified is not found!\n\n";
		}
	}

$step = uc $step;
my %valid_steps = map {$_ => 1} qw(ALL FILTER FINAL ANNO);
die "ERROR: Invalid --step value \"$step\". Valid choices are: all, filter, final, anno.\n" unless $valid_steps{$step};

# --modules: which raw TE discovery modules to run. "plant" is a shortcut that
# skips the SINE and LINE modules (in most plant genomes they annotate <2% of
# the sequence, while their discovery — AnnoSINE and RepeatModeler — costs the
# most wall time). Modules excluded here leave empty library files behind so
# the filter/final/anno stages proceed unchanged.
my %valid_modules = map {$_ => 1} qw/ltr sine line tir helitron/;
$modules = "ltr,tir,helitron" if $modules =~ /^plant$/i;
$modules = "all" if $modules =~ /^all$/i;
if ($modules ne "all"){
	my @mods = split /\s*,\s*/, $modules;
	die "ERROR: Invalid --modules value \"$modules\". Valid choices are: all, plant, or a comma-separated list of ltr, sine, line, tir, helitron.\n"
		unless @mods and not grep { not $valid_modules{$_} } @mods;
	print "\tNote: running raw modules only for: @mods (--modules).\n\n";
	}
goto $step;


############################################################
####### Get raw LTR/SINE/LINE/TIR/Helitron candidates ######
############################################################

ALL:

# report status
chomp ($date = `date`);
print "$date\tObtain raw TE libraries using various structure-based programs: \n";

# Get raw TE candidates
system("EDTA_DUP_CHECK_DONE=1 perl $EDTA_raw --genome $genome --overwrite $overwrite --species $species --type $modules --u $miu --threads $threads --genometools $genometools --ltrretriever $LTR_retriever --blastplus $blastplus --tesorter $TEsorter --GRF $GRF --trf_path $trf --repeatmasker $repeatmasker --repeatmodeler $repeatmodeler --annosine $annosine --tirlearner $TIR_Learner --convert_seq_name 0 --rmlib $RMlib --wholeelement $wholeelement")==0 or die "EDTA_raw.pl failed with exit code ".($? >> 8)."\n";

chdir "$genome.EDTA.raw" or die "Cannot enter $genome.EDTA.raw: $!\n";

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
`sort -T . -sV -k1,1 -k4,4 $genome.EDTA.intact.gff3.temp | grep -v '^#' > $genome.EDTA.intact.raw.gff3; rm $genome.EDTA.intact.gff3.temp`;

chomp ($date = `date`);
print "$date\tObtain raw TE libraries finished.
\tAll intact TEs found by EDTA: \n\t\t$genome.EDTA.intact.raw.fa \n\t\t$genome.EDTA.intact.raw.gff3\n\n";
chdir "..";


############################################################
####### Filter LTR/SINE/LINE/TIR/Helitron candidates #######
############################################################

FILTER:

# report status
chomp ($date = `date`);
print "$date\tPerform EDTA advance filtering for raw TE candidates and generate the stage 1 library: \n\n";

# remove existing results
`rm ./$genome.EDTA.combine/* ./$genome.EDTA.combine/.step*.done 2>/dev/null` if $overwrite == 1;

# Filter raw TE candidates and the make stage 1 library.
# Guarded for restartability (added 2026-08-24): EDTA.pl used to call
# EDTA_processK.pl unconditionally and only check for the stage 1 library
# afterwards, so a chained resubmit re-ran the entire filtering stage even when
# it had already completed. Individual steps inside that stage now resume via
# the .step*.done stamps in $genome.EDTA.combine/.
if (-s "$genome.EDTA.combine/$genome.EDTA.fa.stg1" and $overwrite == 0){
	print "$date\tExisting stage 1 library $genome.EDTA.combine/$genome.EDTA.fa.stg1 found!\n\t\tWill keep this file without rerunning this module.\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
system("perl $EDTA_process -genome $genome -ltr $genome.EDTA.raw/$genome.LTR.raw.fa -ltrint $genome.EDTA.raw/$genome.LTR.intact.raw.fa -line $genome.EDTA.raw/$genome.LINE.raw.fa -sine $genome.EDTA.raw/$genome.SINE.raw.fa -tir $genome.EDTA.raw/$genome.TIR.intact.raw.fa -helitron $genome.EDTA.raw/$genome.Helitron.intact.raw.fa -repeatmasker $repeatmasker -blast $blastplus -threads $threads")==0 or die "EDTA_processK.pl failed with exit code ".($? >> 8)."\n";
	}

# check results, remove intermediate files, and report status
die "ERROR: Stage 1 library not found in $genome.EDTA.combine/$genome.EDTA.fa.stg1" unless -s "$genome.EDTA.combine/$genome.EDTA.fa.stg1";
chdir "$genome.EDTA.combine" or die "Cannot enter $genome.EDTA.combine: $!\n";
`rm ./$genome.LTR.raw.fa*Q* ./$genome.LTR.intact.raw.fa*Q* ./$genome.TIR.intact.raw.fa*Q* ./$genome.Helitron.intact.raw.fa*Q* ./$genome.TIR.Helitron.fa*Q* $genome*tbl $genome*out $genome*cleanup $genome*RMoutput $genome*stg1.raw* $genome.LTR.raw.fa-* $genome.LTR.intact.raw.fa-* $genome.TIR.intact.raw.fa-* $genome.Helitron.intact.raw.fa-* $genome.LINE_LTR.raw.fa $genome.LTR.SINE.LINE.fa *.ndb *.not *.ntf *.nto *.cat.gz *.cat *.masked *.ori.out *.nhr *.nin *.nsq 2>/dev/null` unless $debug eq 1;

chdir "..";
chomp ($date = `date`);
print "$date\tEDTA advance filtering finished.\n\n";


####################################
###### Final purge CDS in TEs ######
####################################

FINAL:

# whole-stage resume guard: the final stage must be rebuilt all-or-nothing. A partial
# rebuild (e.g. re-running the combine steps while skipping cleanup_nested/rename_TE
# via their own guards) would mix renamed and unrenamed states and change the library.
if ($overwrite == 0 and -s "$genome.EDTA.TElib.fa" and -s "$genome.EDTA.intact.fa" and -s "$genome.EDTA.intact.gff3"){
	chomp ($date = `date`);
	print "$date\tExisting final TE library and intact TEs found, skipping the final stage (--overwrite 0).\n\n";
	goto ANNO;
	}

# report status
chomp ($date = `date`);
print "$date\tPerform EDTA final steps to generate a non-redundant comprehensive TE library.\n\n";

# Make the final working directory
`mkdir $genome.EDTA.final` unless -e "$genome.EDTA.final" && -d "$genome.EDTA.final";
die "Cannot create directory $genome.EDTA.final: $!\n" unless -d "$genome.EDTA.final";
chdir "$genome.EDTA.final" or die "Cannot enter $genome.EDTA.final: $!\n";
`rm ./* 2>/dev/null` if $overwrite == 1;
`cp ../$genome.EDTA.raw/$genome.RM2.fa ./`;
`cp ../$genome.EDTA.combine/$genome.EDTA.fa.stg1 ./`;
`cp ../$cds ./` if $cds ne '';
`cp ../$HQlib ./` if $HQlib ne '';
`cp ../$genome.EDTA.combine/$genome.EDTA.intact.fa.cln ./$genome.EDTA.intact.fa.cln`;
#`cp ../$genome.EDTA.raw/$genome.EDTA.intact.raw.fa ./`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.raw.gff3 ./`;
`cp ../$exclude ./` if $exclude ne '';

# identify remaining TEs in the filtered RM2 library
if ($sensitive == 1 and -s "$genome.RM2.fa"){
	print "\tFilter RepeatModeler results that are ignored in the raw step.\n\n";
	chomp ($date = `date`);
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -nolow -div 40 -lib $genome.EDTA.fa.stg1 $genome.RM2.fa 2>&1`;
	my $rm_exit = $? >> 8;
	die "ERROR: RepeatMasker failed on $genome.RM2.fa (exit code $rm_exit):\n$rm_status\n\n" if $? != 0 and not -e "$genome.RM2.fa.masked";
	`cp $genome.RM2.fa $genome.RM2.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`cp $genome.RM2.fa $genome.RM2.fa.masked` unless -e "$genome.RM2.fa.masked";
	# clean up tandem and coding sequences in the RM2 library
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.RM2.fa.masked > $genome.RM2.fa.stg1`;
	`perl $cleanup_proteins -seq $genome.RM2.fa.stg1 -rmdnate 0 -rmline 0 -rmprot 1 -protlib $protlib -blast $blastplus -threads $threads`;
	if (-s "$genome.RM2.fa.stg1.clean"){
		`cat $genome.EDTA.fa.stg1 $genome.RM2.fa.stg1.clean > $genome.EDTA.raw.fa`;
		} else {
		print "\t\tNo extra repeat sequences found in the RepeatModeler output.\n\n";
		`cp $genome.EDTA.fa.stg1 $genome.EDTA.raw.fa`;
		}
	} elsif ($sensitive == 1){
	print "\tSkipping the RepeatModeler results: $genome.RM2.fa not found.\n\t\tThis file is generated by the raw step when EDTA is run with --sensitive 1 (--step all).\n\n";
	`cp $genome.EDTA.fa.stg1 $genome.EDTA.raw.fa`;
	} else {
	print "\tSkipping the RepeatModeler results (--sensitive 0).\n\t\tRun EDTA.pl --step final --sensitive 1 if you want to add RepeatModeler results.\n\n";
	`cp $genome.EDTA.fa.stg1 $genome.EDTA.raw.fa`;
	}

# remove CDS in the non-redundant library and intact TEs
if (-s "$cds"){
	# report status
	chomp ($date = `date`);

	# cleanup TE-related sequences in the CDS file with TEsorter
	print "$date\tClean up TE-related sequences in the CDS file with TEsorter.\n\n";
	my $cds_clean_err = `perl $cleanup_TE -cds $cds -minlen 300 -tesorter $TEsorter -repeatmasker $repeatmasker -t $threads -rawlib $genome.EDTA.raw.fa 2>&1`;
	die "ERROR: cleanup_TE failed (exit code ".($? >> 8)."):\n$cds_clean_err\n" if $? != 0;
	`rm ./$cds ./$cds.code.r* 2>/dev/null` unless $debug eq 1;
	die "\tERROR: The $cds file is empty after TE clean up. Please check the file and $cds.code.noTE.\n\n" unless -s "$cds.code.noTE";
	$cds = "$cds.code.noTE";

	# remove cds-related sequences in the EDTA library
	print "\tRemove CDS-related sequences in the EDTA library.\n\n";
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.raw.fa 2>&1`;
	my $rm_exit = $? >> 8;
	die "ERROR: RepeatMasker failed on $genome.EDTA.raw.fa (exit code $rm_exit):\n$rm_status\n\n" if $? != 0 and not -e "$genome.EDTA.raw.fa.masked";
	`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.masked` unless -e "$genome.EDTA.raw.fa.masked";
	`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.3 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.EDTA.raw.fa.masked > $genome.EDTA.raw.fa.cln.tmp.$$ && mv $genome.EDTA.raw.fa.cln.tmp.$$ $genome.EDTA.raw.fa.cln`;
	die "cleanup_tandem failed for $genome.EDTA.raw.fa.cln\n" if $? != 0;

	# remove cds-related sequences in intact TEs
	print "\tRemove CDS-related sequences in intact TEs.\n\n";
	$rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -nolow -div 40 -cutoff 225 -lib $cds $genome.EDTA.intact.fa.cln 2>&1`;
	$rm_exit = $? >> 8;
	die "ERROR: RepeatMasker failed on $genome.EDTA.intact.fa.cln (exit code $rm_exit):\n$rm_status\n\n" if $? != 0 and not -e "$genome.EDTA.intact.fa.cln.masked";
	`cp $genome.EDTA.intact.fa.cln $genome.EDTA.intact.fa.cln.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`cp $genome.EDTA.intact.fa.cln $genome.EDTA.intact.fa.cln.masked` unless -e "$genome.EDTA.intact.fa.cln.masked";
	`perl $cleanup_tandem -misschar N -Nscreen 1 -nc 1000 -nr 0.8 -minlen 80 -maxlen 5000000 -trf 0 -cleanN 0 -f $genome.EDTA.intact.fa.cln.masked > $genome.EDTA.intact.fa.cln.rmCDS`;
	`perl $output_by_list 1 $genome.EDTA.intact.fa.cln 1 $genome.EDTA.intact.fa.cln.masked.cleanup -ex -FA > $genome.EDTA.intact.fa.cln2`;
	} else {
	print "\tSkipping the CDS cleaning step (--cds [File]) since no CDS file is provided or it's empty.\n\n";
	`cp $genome.EDTA.raw.fa $genome.EDTA.raw.fa.cln`;
	`cp $genome.EDTA.intact.fa.cln $genome.EDTA.intact.fa.cln2`;
	}

# Final rounds of redundancy removal and make final EDTA library
# resume guard added 2026-08-28: cleanup_nested took 41.0 h on the 11.1 Gb oat genome
# (job 20121980) and was unguarded, so any wall-clock kill later in FINAL redid all of it.
# A killed cleanup_nested leaves $genome.EDTA.raw.fa.cln.iter* snapshots behind (they are
# unlinked only after the final .cln.cln is fully written), so leftover iter files mark an
# incomplete run: do not let the -s guard accept a truncated .cln.cln.
if (-s "$genome.EDTA.raw.fa.cln.cln" and $overwrite == 0 and !glob "$genome.EDTA.raw.fa.cln.iter*"){
	print "\tExisting $genome.EDTA.raw.fa.cln.cln found, skipping cleanup_nested (--overwrite 0).\n\n";
} else {
	my $nested_err = `perl $cleanup_nested -in $genome.EDTA.raw.fa.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blastplus 2>&1`;
	die "ERROR: cleanup_nested failed (exit code ".($? >> 8)."):\n$nested_err\n" if $? != 0;
}

# rename all TEs in the EDTA library
# resume guard added 2026-08-28: keeps the TE library byte-identical across a resumed
# FINAL, so a library already handed to an external RepeatMasker run stays valid.
if (-s "$genome.EDTA.TElib.fa" and $overwrite == 0
    and (!$wholeelement or -s "$genome.EDTA.TElib.fa.rename_map")){
	print "\tExisting $genome.EDTA.TElib.fa found, skipping rename_TE (--overwrite 0).\n\n";
} elsif ($wholeelement){
	`perl $rename_TE $genome.EDTA.raw.fa.cln.cln --map $genome.EDTA.TElib.fa.rename_map > $genome.EDTA.TElib.fa.tmp.$$ && mv $genome.EDTA.TElib.fa.tmp.$$ $genome.EDTA.TElib.fa`;
	die "rename_TE failed for $genome.EDTA.TElib.fa\n" if $? != 0;
} else {
	`perl $rename_TE $genome.EDTA.raw.fa.cln.cln > $genome.EDTA.TElib.fa.tmp.$$ && mv $genome.EDTA.TElib.fa.tmp.$$ $genome.EDTA.TElib.fa`;
	die "rename_TE failed for $genome.EDTA.TElib.fa\n" if $? != 0;
}
#`perl $rename_TE $genome.EDTA.raw.fa.cln.cln | perl $format_TElib - > $genome.EDTA.TElib.fa`;

# Build the LTR boundary file (id, total_len, lLTR_len, rLTR_len) for the final library.
# Primary path measures boundaries directly from the shipped sequences by self-alignment
# (bin/ltrbound_denovo.py): correct by construction, so it stays right after advance
# filtering trims a library sequence. ~1/3 of LTR entries resolve a terminal-repeat pair;
# the rest genuinely lost one and are correctly omitted. Fallback without python3 reuses
# LTR_retriever's numbers, but only for entries filtering left unchanged.
# NOTE: cwd is $genome.EDTA.final here (chdir above), and $genome.EDTA.raw is its sibling.
if ($wholeelement and -s "$genome.EDTA.TElib.fa"){
	my $py = `command -v python3 2>/dev/null`;
	if (-s $ltrbound_denovo and $py ne ''){
		local $ENV{EDTA_BLASTN} = "${blastplus}blastn";
		`python3 $ltrbound_denovo $genome.EDTA.TElib.fa $genome.EDTA.TElib.LTRbound $threads`;
	} elsif (-s "../$genome.EDTA.raw/$genome.LTRlib.fa.LTRbound" and -s "$genome.EDTA.TElib.fa.rename_map"){
		`perl $update_LTRbound $genome.EDTA.TElib.fa.rename_map ../$genome.EDTA.raw/$genome.LTRlib.fa.LTRbound $genome.EDTA.TElib.fa > $genome.EDTA.TElib.LTRbound`;
	}
}

# identify novel TEs using the user provided $HQlib
if ($HQlib ne ''){
	# report status
	chomp ($date = `date`);
	print "$date\tCombine the high-quality TE library $HQlib with the EDTA library:\n\n";

	# remove known TEs in the EDTA library
	my $rm_status = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -nolow -div 40 -lib $HQlib $genome.EDTA.TElib.fa 2>&1`;
	my $rm_exit = $? >> 8;
	die "ERROR: RepeatMasker failed on $genome.EDTA.TElib.fa (exit code $rm_exit):\n$rm_status\n\n" if $? != 0 and not -e "$genome.EDTA.TElib.fa.masked";
	`cp $genome.EDTA.TElib.fa $genome.EDTA.TElib.fa.masked` if $rm_status =~ /No repetitive sequences were detected/i;
	`cp $genome.EDTA.TElib.fa $genome.EDTA.TElib.fa.masked` unless -e "$genome.EDTA.TElib.fa.masked";
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 0 -f $genome.EDTA.TElib.fa.masked > $genome.EDTA.TElib.novel.fa`;
	rename "$genome.EDTA.TElib.fa", "$genome.EDTA.TElib.ori.fa";
	`cat $HQlib $genome.EDTA.TElib.novel.fa > $genome.EDTA.TElib.fa`;
	copy_file("$genome.EDTA.TElib.novel.fa", "..");
	}

# reclassify intact TEs with the TE lib #113
`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -nolow -div 40 -lib $genome.EDTA.TElib.fa $genome.EDTA.intact.fa.cln2 2>/dev/null` unless (-s "$genome.EDTA.intact.fa.cln2.out" and $overwrite == 0);
die "ERROR: The masked file for $genome.EDTA.intact.fa.cln2 is not found! The RepeatMasker annotation on this file may be failed. Please check the $genome.EDTA.TElib.fa file for sequence naming formats especially when you provide a library via --curatedlib.\n" unless -s "$genome.EDTA.intact.fa.cln2.out";
`perl $reclassify -seq $genome.EDTA.intact.fa.cln2 -RM $genome.EDTA.intact.fa.cln2.out -cov 80 -len 80 -iden 60`; # 80-80-60

# remove inconsistently classified intact TEs and generate the final intact TEs
`perl $output_by_list 1 $genome.EDTA.intact.fa.cln2.rename 1 $genome.EDTA.intact.fa.cln2.false.list -ex -FA > $genome.EDTA.intact.fa`;


## generate clean intact gff3
my $intact_gff_head = "##This file follows the ENSEMBL standard: https://useast.ensembl.org/info/website/upload/gff3.html
##Column 3: Sequence Ontology of repeat features. Please refer to the SO database for more details: http://www.sequenceontology.org/. In cases where the SO database does not have the repeat feature, tentative SO names are used, with a full list included in EDTA/bin/TE_Sequence_Ontology.txt (Enhancement notes), and the sequence_ontology in Column 9 uses the closest parent SO.
##Column 9: 
##      ID: unique ID for this feature in the genome.
##      classification: Same as Column 3 but formatted following the RepeatMasker naming convention.
##      sequence_ontology: Sequence Ontology ID of the feature.
##      identity: Sequence identity (0-1) between terminal sequences for structurally annotated TIR elements.
##      ltr_identity: Sequence identity (0-1) between the left and right LTR regions for structurally annotated LTR elements.
##      Name: Repeat family name. Some may be shown as coordinates, which are single-copy and structrually identified elements that are not included in the repeat library.
##      method=structural: Indicate this entry is produced by structural annotation.
##      motif/TSD/TIR: structural features of structurally annotated LTR and TIR elements.
##For more details about this file, please refer to the EDTA wiki: https://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A
##seqid source sequence_ontology start end score strand phase attributes";

# update the family names in the intact.raw.gff3 file
`perl $rename_by_list $genome.EDTA.intact.raw.gff3 $genome.EDTA.intact.fa.cln2.rename.list 1 > $genome.EDTA.intact.raw.gff3.rename`;
`sed 's/.*Name=//i; s/;Classifica.*//i' $genome.EDTA.intact.raw.gff3.rename | sort -u > $genome.EDTA.intact.raw.gff3.rename.famlist`;

# get a dirty list of intact.gff
`grep -a \\> $genome.EDTA.intact.fa | sed 's/>//; s/#.*//' | perl $output_by_list 1 $genome.EDTA.intact.raw.gff3.rename.famlist 1 - -ex | awk '{print "Name\\t"\$1"\\nParent\\t"\$1"\\nID\\t"\$1}' > $genome.EDTA.intact.raw.gff3.rename.dirtlist`;

# first attempt purging the gff3 (only its .removed side-effect is used, stdout discarded)
`perl $filter_gff $genome.EDTA.intact.raw.gff3.rename $genome.EDTA.intact.raw.gff3.rename.dirtlist > /dev/null`;

# remake the remove list and purge again
`perl -nle 'my \$id = \$1 if /=(repeat_region[0-9]+);/; print "Parent\\t\$id\nName\\t\$id" if defined \$id' $genome.EDTA.intact.raw.gff3.rename.removed >> $genome.EDTA.intact.raw.gff3.rename.dirtlist`;
`echo "##gff-version 3\n##date $date\n##This file contains repeats annotated by EDTA $version based on structural features.\n$intact_gff_head" > $genome.EDTA.intact.gff3`;
`perl $filter_gff $genome.EDTA.intact.raw.gff3.rename $genome.EDTA.intact.raw.gff3.rename.dirtlist >> $genome.EDTA.intact.gff3`;

# format gff3
`perl $format_gff3 $genome.EDTA.intact.gff3 > gff3.temp.$$.gff3 && mv gff3.temp.$$.gff3 $genome.EDTA.intact.gff3`;
die "format_gff3 failed to produce $genome.EDTA.intact.gff3\n" if $? != 0;

# add TE_IDs to the intact.fa sequence IDs
`perl $add_id -fa $genome.EDTA.intact.fa -gff $genome.EDTA.intact.gff3 > $genome.EDTA.intact.fa.renamed; mv $genome.EDTA.intact.fa.renamed $genome.EDTA.intact.fa`;

# check results
die "ERROR: Final TE library not found in $genome.EDTA.TElib.fa" unless -s "$genome.EDTA.TElib.fa";
die "ERROR: Intact TE annotation not found in $genome.EDTA.intact.gff3" unless -s "$genome.EDTA.intact.gff3";
copy_file("$genome.EDTA.TElib.fa", "..");
copy_file("$genome.EDTA.TElib.LTRbound", "..") if $wholeelement and -s "$genome.EDTA.TElib.LTRbound";
copy_file("$genome.EDTA.intact.fa", "..");
copy_file("$genome.EDTA.intact.gff3", "..");

# Decode sequence IDs in user-facing output files (parent directory copies)
if ($seqid_mapfile ne '' and -s $seqid_mapfile){
	unlink "../$genome.EDTA.intact.fa.decoded";
	system("perl $seqid_codec decode_fasta ../$genome.EDTA.intact.fa $seqid_mapfile ../$genome.EDTA.intact.fa.decoded && mv ../$genome.EDTA.intact.fa.decoded ../$genome.EDTA.intact.fa")==0 or die "Failed to decode ../$genome.EDTA.intact.fa: $?\n";
	unlink "../$genome.EDTA.intact.gff3.decoded";
	system("perl $seqid_codec decode_text ../$genome.EDTA.intact.gff3 $seqid_mapfile ../$genome.EDTA.intact.gff3.decoded && mv ../$genome.EDTA.intact.gff3.decoded ../$genome.EDTA.intact.gff3")==0 or die "Failed to decode ../$genome.EDTA.intact.gff3: $?\n";
	unlink "../$genome.EDTA.TElib.fa.decoded";
	system("perl $seqid_codec decode_fasta ../$genome.EDTA.TElib.fa $seqid_mapfile ../$genome.EDTA.TElib.fa.decoded && mv ../$genome.EDTA.TElib.fa.decoded ../$genome.EDTA.TElib.fa")==0 or die "Failed to decode ../$genome.EDTA.TElib.fa: $?\n";
	if (-s "../$genome.EDTA.TElib.novel.fa"){
		unlink "../$genome.EDTA.TElib.novel.fa.decoded";
		system("perl $seqid_codec decode_fasta ../$genome.EDTA.TElib.novel.fa $seqid_mapfile ../$genome.EDTA.TElib.novel.fa.decoded && mv ../$genome.EDTA.TElib.novel.fa.decoded ../$genome.EDTA.TElib.novel.fa")==0 or die "Failed to decode ../$genome.EDTA.TElib.novel.fa: $?\n";
		}
	}

# remove intermediate files, but keep the resume-guard files $genome.EDTA.raw.fa.cln.cln
# and $genome.EDTA.TElib.fa.rename_map that the rm globs would otherwise delete; the
# keep-prefixed temp names do not match any of the globs below
if ($debug eq 0){
	rename "$genome.EDTA.raw.fa.cln.cln", "keep.$genome.EDTA.raw.fa.cln.cln" if -e "$genome.EDTA.raw.fa.cln.cln";
	rename "$genome.EDTA.TElib.fa.rename_map", "keep.$genome.EDTA.TElib.fa.rename_map" if -e "$genome.EDTA.TElib.fa.rename_map";
	`rm $genome.EDTA.intact.fa.cln.* $genome.EDTA.raw.fa.* $genome.EDTA.TElib.fa.* $genome.LTR.TIR.Helitron.fa.stg1.* $genome.masked *.cat.gz 2>/dev/null`;
	rename "keep.$genome.EDTA.raw.fa.cln.cln", "$genome.EDTA.raw.fa.cln.cln" if -e "keep.$genome.EDTA.raw.fa.cln.cln";
	rename "keep.$genome.EDTA.TElib.fa.rename_map", "$genome.EDTA.TElib.fa.rename_map" if -e "keep.$genome.EDTA.TElib.fa.rename_map";
	}

# report status
chomp ($date = `date`);
print "$date\tEDTA final stage finished! You may check out:
		The final EDTA TE library: $genome.EDTA.TElib.fa\n";
print "		Family names of intact TEs have been updated by $HQlib: $genome.EDTA.intact.gff3\n" if $HQlib ne '';
print "		Comparing to the provided library, EDTA found these novel TEs: $genome.EDTA.TElib.novel.fa
		The provided library has been incorporated into the final library: $genome.EDTA.TElib.fa\n\n" if $HQlib ne '';
chdir "..";

# Decode sequence IDs in user-facing deliverables only; intermediate files in the
# working directories stay encoded so partial-rerun resume state is never mixed.
# decode_text is safe on FASTA files because DNA sequence lines cannot match the _J code pattern.
if ($seqid_mapfile ne '' and -s $seqid_mapfile){
	for my $f ("$genome.EDTA.final/$genome.EDTA.TElib.fa", "$genome.EDTA.final/$genome.EDTA.intact.fa", "$genome.EDTA.final/$genome.EDTA.intact.gff3", "$genome.EDTA.final/$genome.EDTA.TElib.novel.fa"){
		next unless -s $f;
		unlink "$f.decoded";
		system("perl $seqid_codec decode_text $f $seqid_mapfile $f.decoded && mv $f.decoded $f")==0 or die "Failed to decode $f: $?\n";
		}
	# Also decode RM2.raw.fa in the parent directory
	if (-s "$genome.RM2.raw.fa"){
		unlink "$genome.RM2.raw.fa.decoded";
		system("perl $seqid_codec decode_text $genome.RM2.raw.fa $seqid_mapfile $genome.RM2.raw.fa.decoded && mv $genome.RM2.raw.fa.decoded $genome.RM2.raw.fa")==0 or die "Failed to decode $genome.RM2.raw.fa: $?\n";
		}
	}


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
	die "Cannot create directory $genome.EDTA.anno: $!\n" unless -d "$genome.EDTA.anno";
	chdir "$genome.EDTA.anno" or die "Cannot enter $genome.EDTA.anno: $!\n";
	`rm ./* 2>/dev/null` if $overwrite == 1;
	`rm $genome.EDTA.TElib.fa* 2>/dev/null`; # clean up libraries
	`cp ../$genome.EDTA.TElib.fa ./`;
	`cp ../$genome.EDTA.intact.gff3 ./`;
	`cp ../$exclude ./` if $exclude ne '';
	`ln -s ../$genome $genome` unless -e $genome;

	my $gff_head = "##This file follows the ENSEMBL standard: https://useast.ensembl.org/info/website/upload/gff3.html
##Column 3: Sequence Ontology of repeat features. Please refer to the SO database for more details: http://www.sequenceontology.org/. In cases where the SO database does not have the repeat feature, tentative SO names are used, with a full list included in EDTA/bin/TE_Sequence_Ontology.txt (Enhancement notes), and the sequence_ontology in Column 9 uses the closest parent SO.
##Column 7: The Smith-Waterman score generated by RepeatMasker, only available for homology entries.
##Column 9: 
##	ID: unique ID for this feature in the genome.
##	classification: Same as Column 3 but formatted following the RepeatMasker naming convention.
##	sequence_ontology: Sequence Ontology ID of the feature.
##	identity: Sequence identity (0-1) between the library sequence and the target region.
##	ltr_identity: Sequence identity (0-1) between the left and right LTR regions for structurally annotated LTR elements.
##	Name: Repeat family name. Some may be shown as coordinates, which are single-copy and structrually identified elements that are not included in the repeat library.
##	method: Indicate if this entry is produced by structural annotation or homology annotation.
##	motif/TSD/TIR: structural features of structurally annotated LTR and TIR elements.
##For more details about this file, please refer to the EDTA wiki: https://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A
##seqid source sequence_ontology start end score strand phase attributes";

	# annotate TEs using RepeatMasker
	if ($rmout ne ''){
		print STDERR "$date\tA RepeatMasker result file $rmout is provided! Will use this file without running RepeatMasker.\n\n";
		if (-e "$genome.out"){
			my $old_rmout = `ls -l $genome.out|perl -nle 'my (\$month, \$day, \$time) = (split)[5,6,7]; \$time =~ s/://; print "\${month}_\${day}_\$time"'`;
			chomp $old_rmout;
			print "\t$genome.out exists in the $genome.EDTA.anno folder, renamed file to ${genome}_$old_rmout.out\n\n";
			`mv $genome.out ${genome}_$old_rmout.out`;
			}
		`ln -s $rmout $genome.out`;
		} else {
		print STDERR "$date\tHomology-based annotation of TEs using $genome.EDTA.TElib.fa from scratch.\n\n";
		my $rm_anno_err;
		$rm_anno_err = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -nolow -div $maxdiv -lib $genome.EDTA.TElib.fa $genome 2>&1` unless (-s "$genome.out" and $overwrite == 0);
		die "ERROR: RepeatMasker failed on $genome (exit code ".($? >> 8)."):\n$rm_anno_err\n\n" if defined $rm_anno_err and $? != 0 and not -s "$genome.out";
		}
	die "ERROR: RepeatMasker results not found in $genome.out!\n\n" unless -s "$genome.out";

	# Decode RepeatMasker output and genome so all downstream processing produces decoded IDs
	if ($seqid_mapfile ne '' and -s $seqid_mapfile){
		unlink "$genome.out.decoded";
		system("perl $seqid_codec decode_text $genome.out $seqid_mapfile $genome.out.decoded && mv $genome.out.decoded $genome.out")==0 or die "Failed to decode $genome.out: $?\n";
		# Replace genome symlink with decoded FASTA (downstream scripts read seq IDs from the genome);
		# decode to a temp name first so a failed decode cannot leave the anno dir without a genome.
		# The stamp skips the genome-scale rewrite on --overwrite 0 resumes; it is written only
		# after the decoded genome has been moved into place.
		unless (-e ".genome_decoded.stamp" and $overwrite == 0){
			unlink $genome, "./$genome.decoded";
			if (system("perl $seqid_codec decode_fasta ../$genome $seqid_mapfile ./$genome.decoded") == 0 and -s "./$genome.decoded"){
				rename "./$genome.decoded", "./$genome" or die "ERROR: Failed to move the decoded genome ./$genome.decoded into place: $!\n";
				open my $genome_stamp, ">", ".genome_decoded.stamp" or die "ERROR: Cannot create .genome_decoded.stamp: $!\n";
				close $genome_stamp;
				} else {
				die "ERROR: Failed to decode sequence IDs in $genome (exit code ".($? >> 8)."). The encoded genome is still available at ../$genome\n";
				}
			}
		}

	# exclude regions from TE annotation and make whole-genome TE annotation
	`perl $make_masked -genome $genome -rmout $genome.out -maxdiv $maxdiv -minscore 300 -minlen 80 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
	# combine RepeatMasker lines that appears to be the same element
	`perl $combine_RMrows -rmout $genome.out.new -maxdiv 3.5 -maxgap 35`;
	`mv $genome.out.new.cmb $genome.EDTA.RM.out`;
	die "combine_RMrows failed to produce $genome.EDTA.RM.out\n" unless -s "$genome.EDTA.RM.out";
	`perl $RMout2bed $genome.EDTA.RM.out > $genome.EDTA.RM.bed`; # a regular enriched bed
	`perl $bed2gff $genome.EDTA.RM.bed TE_homo > $genome.EDTA.RM.gff3`;
	`perl $gff2bed $genome.EDTA.RM.gff3 homology > $genome.EDTA.RM.bed`; # add the last column to this bed

	# combine homology-based and strutrual-based annotation (partly overlapping)
	`perl $gff2bed $genome.EDTA.intact.gff3 structural > $genome.EDTA.intact.bed`;
	`perl $combine_overlap $genome.EDTA.intact.bed $genome.EDTA.intact.bed.cmb 5`;
	`perl $get_frag $genome.EDTA.RM.bed $genome.EDTA.intact.bed.cmb $threads`;
	`perl $keep_nest $genome.EDTA.intact.bed $genome.EDTA.RM.bed $threads`;
	`grep homology $genome.EDTA.intact.bed-$genome.EDTA.RM.bed > $genome.EDTA.intact.bed-$genome.EDTA.RM.bed.homo`;
	`sort -T . -suV $genome.EDTA.intact.bed-$genome.EDTA.RM.bed.homo $genome.EDTA.RM.bed-$genome.EDTA.intact.bed.cmb > $genome.EDTA.homo.bed`;
	`perl $bed2gff $genome.EDTA.homo.bed TE_homo > $genome.EDTA.homo.gff3`;
	`cat $genome.EDTA.intact.gff3 $genome.EDTA.homo.gff3 > $genome.EDTA.TEanno.gff3.raw`;
	# write the header first and append the sorted body, instead of slurping the whole GFF into memory
	chomp (my $anno_date = `date`);
	$anno_date =~ s/\s+$//;
	`printf "##gff-version 3\n##date $anno_date\n##This file contains repeats annotated by EDTA $version with both structural and homology methods. Repeats can be overlapping due to nested insertions.\n$gff_head\n" > $genome.EDTA.TEanno.gff3`;
	`grep -v '^#' $genome.EDTA.TEanno.gff3.raw | sort -T . -sV -k1,1 -k4,4 >> $genome.EDTA.TEanno.gff3`;
	`perl $format_gff3 $genome.EDTA.TEanno.gff3 > gff3.temp.$$.gff3 && mv gff3.temp.$$.gff3 $genome.EDTA.TEanno.gff3`;
	die "format_gff3 failed to produce $genome.EDTA.TEanno.gff3\n" if $? != 0;
	`perl $gff2gtf --gff $genome.EDTA.TEanno.gff3 --remove repeat_region,long_terminal_repeat,target_site_duplication --out $genome.EDTA.TEanno.gtf`;
	`rm $genome.EDTA.TEanno.gff3.raw 2>/dev/null`;

	# make non-overlapping annotation
	`perl $gff2bed $genome.EDTA.TEanno.gff3 structural > $genome.EDTA.TEanno.bed`;
	`perl $split_overlap $genome.EDTA.TEanno.bed $genome.EDTA.TEanno.split.bed`;
	`echo "##gff-version 3\n##date $date\n##This file contains all repeats annotated by EDTA $version in the split format (non-overlapping). Repeats can be broken into pieces by nested insertions.\n$gff_head" > $genome.EDTA.TEanno.split.gff3`;
	`perl $bed2gff $genome.EDTA.TEanno.split.bed | grep -v '^#' >> $genome.EDTA.TEanno.split.gff3`;
	`perl $format_gff3 $genome.EDTA.TEanno.split.gff3 > gff3.temp.$$.gff3 && mv gff3.temp.$$.gff3 $genome.EDTA.TEanno.split.gff3`;
	die "format_gff3 failed to produce $genome.EDTA.TEanno.split.gff3\n" if $? != 0;
	`perl $gff2gtf --gff $genome.EDTA.TEanno.split.gff3 --remove repeat_region,long_terminal_repeat,target_site_duplication --out $genome.EDTA.TEanno.split.gtf`;
	`perl $gff2RMout $genome.EDTA.TEanno.split.gff3 $genome.EDTA.TEanno.split.out`;

	# make plots
	`perl $div_table $genome.EDTA.TEanno.bed $genome $genome`;
	my $div_plot_err = `Rscript $div_plot $genome.div_long $genome 2>&1`;
	print STDERR "Warning: divergence plot failed (Rscript exit code ".($? >> 8)."), continuing without plots:\n$div_plot_err\n" if $? != 0;
	`python3 $density_table -genome $genome -gff $genome.EDTA.TEanno.split.gff3 > $genome.EDTA.TEanno.split.density`;
	my $density_plot_err = `Rscript $density_plot $genome.EDTA.TEanno.split.density 2>&1`;
	print STDERR "Warning: density plot failed (Rscript exit code ".($? >> 8)."), continuing without plots:\n$density_plot_err\n" if $? != 0;
	`mv chromosome_density_plots.pdf $genome.EDTA.TEanno.density_plots.pdf`;

	# make summary table for the non-overlapping annotation
	`perl $count_base $genome > $genome.stats`;
	`perl -nle 'my (\$chr, \$s, \$e, \$anno, \$dir, \$supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 \$chr \$s \$e NA \$dir \$anno \$supfam"' $genome.EDTA.TEanno.split.bed > $genome.EDTA.TEanno.out`;
	my $summary_err = `perl $buildSummary -maxDiv $maxdiv -stats $genome.stats $genome.EDTA.TEanno.out 2>&1 > $genome.EDTA.TEanno.sum`;
	die "ERROR: buildSummary failed (exit code ".($? >> 8)."):\n$summary_err\n\n" if $? != 0;
	my $tot_TE = `grep Total $genome.EDTA.TEanno.sum|grep %|awk '{print \$4}'`;
	chomp $tot_TE;

	# make low-threshold masked genome for MAKER (optional; off by default, enable with --maker 1)
	my $maker_TE = '';
	if ($maker == 1){
		# only promote $genome.new.masked when this make_masked call actually ran;
		# the unguarded make_masked above also writes $genome.new.masked with different thresholds
		unless (-s "$genome.MAKER.masked" and $overwrite == 0){
			`perl $make_masked -genome $genome -rmout $genome.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 1 -misschar N -threads $threads -exclude $exclude`;
			`mv $genome.new.masked $genome.MAKER.masked` if -s "$genome.new.masked";
			}
		$maker_TE = `perl $count_base $genome.MAKER.masked`;
		my ($maker_ratio) = (split /\s+/, $maker_TE)[3];
		if (defined $maker_ratio and $maker_ratio =~ /^[0-9.eE+-]+$/){
			$maker_TE = sprintf("%.2f%%", $maker_ratio*100);
			} else {
			warn "\tWARNING: could not parse the masking ratio from the $count_base output for $genome.MAKER.masked\n";
			$maker_TE = 0;
			}
	}

	# check results and report status
	chomp (my $anno_body = `grep -vc '^#' $genome.EDTA.TEanno.gff3`);
	$anno_body = 0 unless $anno_body =~ /^\d+$/;
	die "ERROR: TE annotation results not found in $genome.EDTA.TEanno.gff3!\n\n" if $anno_body == 0;
	print "ERROR: The masked genome for MAKER annotation is not found in $genome.MAKER.masked!\n\n" if ($maker == 1 and !-s "$genome.MAKER.masked");
	chomp ($date = `date`);
	print "$date\tTE annotation using the EDTA library has finished! Check out:\n";
	print "\t\tWhole-genome TE annotation (total TE: $tot_TE): $genome.EDTA.TEanno.gff3 $genome.EDTA.TEanno.gtf\n";
	print "\t\tWhole-genome TE annotation summary: $genome.EDTA.TEanno.sum\n";
	print "\t\tWhole-genome TE divergence plot: ${genome}_divergence_plot.pdf\n";
	print "\t\tWhole-genome TE density plot: $genome.EDTA.TEanno.density_plots.pdf\n";
	print "\t\tLow-threshold TE masking for MAKER gene annotation (masked: $maker_TE): $genome.MAKER.masked\n\n" if $maker == 1;

	# copy results out
	`cp $genome.MAKER.masked ../` if $maker == 1; # make no backup for this file
	copy_file("$genome.EDTA.TEanno.gff3", "..");
	copy_file("$genome.EDTA.TEanno.gtf", "..");
	copy_file("$genome.EDTA.TEanno.sum", "..");
	copy_file("${genome}_divergence_plot.pdf", "..");
	copy_file("$genome.EDTA.TEanno.density_plots.pdf", "..");

	# Decode user-facing deliverables in the EDTA.anno directory and parent directory copies
	if ($seqid_mapfile ne '' and -s $seqid_mapfile){
		for my $f (glob "$genome.EDTA.TEanno.*"){
			next unless -f $f and $f !~ /\.pdf$/;
			unlink "$f.decoded";
			system("perl $seqid_codec decode_text $f $seqid_mapfile $f.decoded && mv $f.decoded $f")==0 or die "Failed to decode $f: $?\n";
			}
		# Decode parent directory copies
		if (-s "../$genome.MAKER.masked"){
			unlink "../$genome.MAKER.masked.decoded";
			system("perl $seqid_codec decode_fasta ../$genome.MAKER.masked $seqid_mapfile ../$genome.MAKER.masked.decoded && mv ../$genome.MAKER.masked.decoded ../$genome.MAKER.masked")==0 or die "Failed to decode ../$genome.MAKER.masked: $?\n";
			}
		for my $decode_f ("$genome.EDTA.TEanno.gff3", "$genome.EDTA.TEanno.gtf", "$genome.EDTA.TEanno.sum"){
			next unless -s "../$decode_f";
			unlink "../$decode_f.decoded";
			system("perl $seqid_codec decode_text ../$decode_f $seqid_mapfile ../$decode_f.decoded && mv ../$decode_f.decoded ../$decode_f")==0 or die "Failed to decode ../$decode_f: $?\n";
			}
		}

	# evaluate the annotation consistency
	if ($evaluate == 1){
		# report status
		chomp ($date = `date`);
		print "$date\tEvaluate the level of inconsistency for whole-genome TE annotation:\n\n";

		# extract whole-genome TE, all-v-all blast, and summarize consistency.
		# Single source of truth: bin/evaluation.pl (the all-v-all blast is checkpointed and resumes
		# partial runs). -out preserves the $genome.EDTA.TE.fa* output naming reported below.
		my $eval_err = `perl $evaluation -genome $genome -anno $genome.EDTA.TEanno.out -out $genome.EDTA.TE.fa -maxcount 100000 -mincov 0.95 -blast $blastplus -threads $threads -overwrite $overwrite 2>&1`;
		die "ERROR: evaluation.pl failed (exit code ".($? >> 8)."):\n$eval_err\n\n" if $? != 0;

		# check results and report status
		die "ERROR: TE annotation stats results not found in $genome.EDTA.TE.fa.stat!\n\n" unless -s "$genome.EDTA.TE.fa.stat";
		chomp ($date = `date`);
		print "$date\tEvaluation of TE annotation finished! Check out these files:
		(each contains unfiltered + divergence-titrated reports at <=40%, <=30%, <=20%, <=10%, <=5%)\n
		Overall: $genome.EDTA.TE.fa.stat.all.sum
		Nested: $genome.EDTA.TE.fa.stat.nested.sum
		Non-nested: $genome.EDTA.TE.fa.stat.redun.sum\n\n";
		}

	print "\t\tIf you want to learn more about the formatting and information of these files, please visit:
	\t\thttps://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A\n\n";

	}


# clean up the run-private scratch dir: its contents are temporary by
# definition (tool caches, sort spill), so remove it wholesale on the
# natural exit path
if ($edta_own_tmp and -d $ENV{TMPDIR}){
	rmtree($ENV{TMPDIR}, { error => \my $err } );
	}


##########################
###### Subroutines #######
##########################

sub copy_file {
	my ($file, $path) = ($_[0], $_[1]);
	# Generate new name with the last modified date and time
	my $mod_time = (stat($file))[9];
	my $new_name = $file . "_" . strftime("%Y%m%d_%H%M%S", localtime($mod_time));
	
	# resolve symlinks and existing files
	if (-l "$path/$file") {
		# File is a symbolic link
		unlink "$path/$file" or die "ERROR: Failed to remove symbolic link for $path/$file\n\n";
		} elsif (-f "$path/$file") {
        	# File is a regular file
        	rename "$path/$file", "$path/$new_name" or die "ERROR: Failed to rename file: $path/$file\n\n";
        	}

        # copy file to the path if it's not the current path; die if an existing file
	# fails to land (e.g. full disk), but skip silently when the source is absent
	# (e.g. optional plot PDFs when R is missing)
	if (abs_path($path) ne abs_path('.') and -e $file){
		`cp $file $path`;
		die "ERROR: Failed to copy $file to $path\n" if $? != 0 or -z "$path/$file";
		}
	}

sub softlink_file {
	my ($src, $dst) = ($_[0], $_[1]);
	# check "same file" BEFORE touching anything: a valid symlink (or the real file
	# itself) named $dst in the working directory must be kept, not unlinked and
	# re-linked to itself (ln -s name name => ELOOP). Only a dangling symlink or a
	# link pointing elsewhere gets replaced; a different real file is fatal.
	if (-e $dst or -l $dst){
		my ($a, $b) = (eval { abs_path($dst) }, eval { abs_path($src) });
		if (defined $a and defined $b and $a eq $b){
			return; # dst already resolves to src: nothing to do
			}
		if (-l $dst){
			unlink $dst; # dangling or pointing elsewhere: replace below
			} else {
			die "ERROR: $dst already exists in the working directory and is not $src.\n\tPlease remove it or run EDTA in a clean directory.\n\n";
			}
		}
	my $src_abs = abs_path($src);
	$src_abs = $src unless defined $src_abs;
	`ln -s $src_abs $dst`;
	die "ERROR: failed to create softlink $dst -> $src: $!\n" unless -e $dst;
	}
