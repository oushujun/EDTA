#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use File::Spec; # for obtaining the real path of a file
use Cwd qw(abs_path); # for resolving the genome softlink
use File::Path qw(rmtree); # for run-private scratch cleanup
use Pod::Usage;

########################################################
##### Perform initial searches for TE candidates    ####
##### Shujun Ou (shujun.ou.1@gmail.com, 07/16/2020) ####
########################################################

## Input:
#	$genome

## Output:
#	$genome.LTR.raw.fa, $genome.LTR.intact.raw.fa, $genome.LTR.intact.raw.gff3
#	$genome.TIR.intact.fa, $genome.TIR.intact.raw.gff3
#	$genome.Helitron.intact.raw.fa, $genome.Helitron.intact.raw.gff3
#	$genome.LINE.raw.fa, $genome.SINE.raw.fa

my $usage = "\nObtain raw TE libraries using various structure-based programs

perl EDTA_raw.pl [options]
	--genome	[File]	The genome FASTA
	--species [rice|maize|others]	Specify the species for identification
					of TIR candidates. Default: others
	--type	[ltr|sine|line|tir|helitron|all|comma list]
					Specify which type of raw TE candidates
					you want to get. Default: all. Accepts a
					comma-separated subset, e.g. ltr,tir,helitron
					(passed by EDTA.pl --modules plant). Excluded
					modules leave empty library files behind.
	--rmlib	[FASTA]	The RepeatModeler library, classified output.
	--overwrite	[0|1]	If previous results are found, decide to
				overwrite (1, rerun) or not (0, default).
	--convert_seq_name	[0|1]	Convert long sequence name to <= 15
					characters and remove annotations (1,
					default) or use the original (0)
	--u [float]	Neutral mutation rate to calculate the age of intact LTR elements.
			Intact LTR age is found in this file: *EDTA_raw/LTR/*.pass.list.
			Default: 1.3e-8 (per bp per year, from rice).
	--genometools	[path]	Path to the GenomeTools program. (default: find from ENV)
	--annosine	[path]	Path to the AnnoSINE program. (default: find from EDTA/bin)
	--ltrretriever	[path]	Path to the LTR_retriever program. (default: find from ENV)
	--blastplus	[path]	Path to the BLAST+ program. (default: find from ENV)
	--tesorter	[path]	Path to the TEsorter program. (default: find from ENV)
	--GRF		[path]	Path to the GRF program. (default: find from ENV)
	--trf_path	[path]	Path to the TRF program. (default: find from ENV)
	--mdust		[path]	Path to the mdust program. (default: find from ENV)
	--repeatmasker	[path]	Path to the RepeatMasker program. (default: find from ENV)
	--repeatmodeler	[path]	Path to the RepeatModeler2 program. (default: find from ENV)
	--wholeelement	[0|1]	Keep LTR-RTs as whole elements instead of
				splitting into LTR/INT regions. Default: 0
	--threads|-t	[int]	Number of theads to run this script. Default: 4
	--parallel_modules	[0|1|2]	Run the raw TE modules (LTR, SINE, LINE, TIR,
				Helitron) sequentially (0), all concurrently with
				weighted thread budgets (1), or staged (2, default):
				light modules first with all threads, then the
				heaviest module (LINE) alone with the full budget.
				Parallel modes are used only when more than one
				module is active and threads >= number of modules.
	--module_weights	[str]	Relative expected durations of the raw modules,
				used to size per-module thread budgets so modules
				finish together (bigger weight = longer module).
				Default ltr=2.5,sine=1.2,line=4,tir=1.5,helitron=1.
	--help|-h	Display this help info
\n";

# pre-defined
my $genome = '';
my $species = 'others';
my $type = 'all';
my $RMlib = 'null';
my $overwrite = 0; #0, no rerun. 1, rerun even old results exist.
my $convert_name = 1; #0, use original seq names; 1 shorten names.
my $wholeelement = 0; #0, split LTR library into LTR/INT (default); 1, keep whole elements
my $maxint = 5000; #maximum interval length (bp) between TIRs (for GRF in TIR-Learner)
my $miu = 1.3e-8; #mutation rate, per bp per year, from rice
my $threads = 4;
my $parallel_modules = 2; #0, run the raw TE modules sequentially; 1, all concurrently (weighted budgets); 2, staged: light modules first, then the heaviest alone with all threads
my $module_threads = $threads; #thread budget for one module, set before modules run
my $module_weights = ''; #user override, e.g. "ltr=2.5,tir=1.5,helitron=1,sine=1.2,line=4"
my $script_path = $FindBin::Bin;
my $cleanup_misclas = "$script_path/bin/cleanup_misclas.pl";
my $get_range = "$script_path/bin/get_range.pl";
my $rename_LTR = "$script_path/bin/rename_LTR_skim.pl";
my $rename_RM = "$script_path/bin/rename_RM_TE.pl";
my $filter_gff = "$script_path/bin/filter_gff3.pl";
my $rename_tirlearner = "$script_path/bin/rename_tirlearner.pl";
my $call_seq = "$script_path/bin/call_seq_by_list.pl";
my $seqid_codec = "$script_path/bin/seqid_codec.pl";
my $output_by_list = "$script_path/bin/output_by_list.pl";
my $cleanup_tandem = "$script_path/bin/cleanup_tandem.pl";
my $get_ext_seq = "$script_path/bin/get_ext_seq.pl";
my $format_helitronscanner = "$script_path/bin/format_helitronscanner_out.pl";
my $flank_filter = "$script_path/bin/flanking_filter.pl";
my $make_bed = "$script_path/bin/make_bed_with_intact.pl";
my $bed2gff = "$script_path/bin/bed2gff.pl";
my $genometools = ''; #path to the genometools program
my $repeatmasker = ''; #path to the RepeatMasker program
my $repeatmodeler = ''; #path to the RepeatModeler program
my $LTR_retriever = ''; #path to the LTR_retriever program
my $TEsorter = ''; #path to the TEsorter program
my $blastplus = ''; #path to the blastn program
my $mdust = ''; #path to mdust
my $trf = ''; #path to trf
my $GRF = ''; #path to GRF
my $annosine = ""; #path to the AnnoSINE program
my $TIR_Learner = ""; #path to TIR-Learner program
my $LTR_FINDER = ""; #path to LTR_FINDER_parallel program  #tianyulu
my $LTR_HARVEST = ""; #path to LTR_HARVEST_parallel program  #tianyulu
my $HelitronScanner = ""; #path to HelitronScanner program  #tianyulu
my $HelitronScanner_Runner = "$script_path/bin/run_helitron_scanner.py";

my $help = undef;

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^--genome$/i and $ARGV[$k+1] !~ /^-/;
	$species = $ARGV[$k+1] if /^--species$/i and $ARGV[$k+1] !~ /^-/;
	$type = lc $ARGV[$k+1] if /^--type$/i and $ARGV[$k+1] !~ /^-/;
	$RMlib = $ARGV[$k+1] if /^--rmlib$/i and $ARGV[$k+1] !~ /^-/;
	$overwrite = $ARGV[$k+1] if /^--overwrite$/i and $ARGV[$k+1] !~ /^-/;
	$convert_name = $ARGV[$k+1] if /^--convert_seq_name$/i and $ARGV[$k+1] !~ /^-/;
	$miu = $ARGV[$k+1] if /^--u$/i and $ARGV[$k+1] !~ /^-/;
	$genometools = $ARGV[$k+1] if /^--genometools/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^--repeatmasker$/i and $ARGV[$k+1] !~ /^-/;
	$repeatmodeler = $ARGV[$k+1] if /^--repeatmodeler$/i and $ARGV[$k+1] !~ /^-/;
	$annosine = $ARGV[$k+1] if /^--annosine$/i and $ARGV[$k+1] !~ /^-/;
	$TIR_Learner = $ARGV[$k+1] if /^--tirlearner$/i and $ARGV[$k+1] !~ /^-/;
	$LTR_retriever = $ARGV[$k+1] if /^--ltrretriever/i and $ARGV[$k+1] !~ /^-/;
	$TEsorter = $ARGV[$k+1] if /^--tesorter$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^--blastplus$/i and $ARGV[$k+1] !~ /^-/;
	$mdust = $ARGV[$k+1] if /^--mdust$/i and $ARGV[$k+1] !~ /^-/;
	$trf = $ARGV[$k+1] if /^--trf_path$/i and $ARGV[$k+1] !~ /^-/;
	$GRF = $ARGV[$k+1] if /^--GRF$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^--threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$module_weights = $ARGV[$k+1] if /^--module_weights$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
	$parallel_modules = $ARGV[$k+1] if /^--parallel_modules$/i and $ARGV[$k+1] !~ /^-/;
	$wholeelement = $ARGV[$k+1] if /^--wholeelement$/i and $ARGV[$k+1] !~ /^-/;
	$help = 1 if /^--help$|^-h$/i;
	$k++;
	}

# check files and parameters
if ($help){
	pod2usage( {
		-verbose => 0,
		-exitval => 0,
		-message => "$usage\n" } );
	}

if (!-s $genome){
	pod2usage( {
		-message => "At least 1 parameter is required:\n1) Input fasta file: --genome\n".
		"\n$usage\n\n",
		-verbose => 0,
		-exitval => 2 } );
	}

if ($species){
	$species =~ s/rice/Rice/i;
	$species =~ s/maize/Maize/i;
	$species =~ s/others/others/i;
	die "The expected value for the species parameter is Rice or Maize or others!\n" unless $species eq "Rice" or $species eq "Maize" or $species eq "others";
	}

die "The expected value for the type parameter is all, or a comma-separated list of: ltr, sine, line, tir, helitron!\n"
	unless $type eq "all" or $type =~ /^(ltr|sine|line|tir|helitron)(\s*,\s*(ltr|sine|line|tir|helitron))*$/;

# which raw TE modules will run
my %module_active;
if ($type eq "all"){
	%module_active = map {$_ => 1} qw/ltr sine line tir helitron/;
	} else {
	%module_active = map {$_ => 1} split /\s*,\s*/, $type;
	}
my %skip_module; #modules skipped because a module-specific dependency is missing

# check bolean
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"};
if ($convert_name != 0 and $convert_name != 1){ die "The expected value for the convert_seq_name parameter is 0 or 1!\n"};
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"};
if ($parallel_modules != 0 and $parallel_modules != 1 and $parallel_modules != 2){ die "The expected value for the parallel_modules parameter is 0, 1, or 2!\n"};
if ($miu !~ /^[0-9.eE+-]+$/){ die "The expected value for the u parameter is float value without units!\n"}

# --- TMPDIR isolation: keep descendants off the system /tmp ----------------
# Only set a private scratch when no TMPDIR is in effect at all: an inherited
# TMPDIR (EDTA.pl's isolated scratch, or the user's environment) is kept
# as-is; EDTA_TMPDIR_KEEP=1 keeps whatever the environment provided.
my $raw_own_tmp = 0;
unless ((defined $ENV{EDTA_TMPDIR_KEEP} and $ENV{EDTA_TMPDIR_KEEP} eq '1')
	or (defined $ENV{TMPDIR} and length $ENV{TMPDIR})){
	$ENV{TMPDIR} = abs_path(".")."/.EDTA.raw.tmp.$$";
	$raw_own_tmp = 1;
	mkdir($ENV{TMPDIR}) unless -d $ENV{TMPDIR};
	}

chomp (my $date = `date`);
print STDERR "$date\tEDTA_raw: Check dependencies, prepare working directories.\n\n";

# check files and dependencies
die "The script get_range.pl is not found in $get_range!\n" unless -s $get_range;
die "The script rename_LTR_skim.pl is not found in $rename_LTR!\n" unless -s $rename_LTR;
die "The script filter_gff3.pl is not found in $filter_gff!\n" unless -s $filter_gff;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script get_ext_seq.pl is not found in $get_ext_seq!\n" unless -s $get_ext_seq;
die "The script format_helitronscanner_out.pl is not found in $format_helitronscanner!\n" unless -s $format_helitronscanner;
die "The script flanking_filter.pl is not found in $flank_filter!\n" unless -s $flank_filter;
die "The script bed2gff.pl is not found in $bed2gff!\n" unless -s $bed2gff;
die "The script make_bed_with_intact.pl is not found in $make_bed!\n" unless -s $make_bed;

# GenomeTools
chomp ($genometools=`command -v gt 2>/dev/null`) if $genometools eq '';
$genometools =~ s/\s+$//;
$genometools = dirname($genometools) unless -d $genometools;
$genometools="$genometools/" if $genometools ne '' and $genometools !~ /\/$/;
die "Error: gt is not found in the genometools path $genometools!\n" unless -X "${genometools}gt";
# RepeatMasker
my $rand=int(rand(1000000));
chomp ($repeatmasker=`command -v RepeatMasker 2>/dev/null`) if $repeatmasker eq '';
$repeatmasker =~ s/\s+$//;
$repeatmasker = dirname($repeatmasker) unless -d $repeatmasker;
$repeatmasker="$repeatmasker/" if $repeatmasker ne '' and $repeatmasker !~ /\/$/;
die "Error: RepeatMasker is not found in the RepeatMasker path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";
`cp \"$script_path/database/dummy060817.fa\" ./dummy060817.fa.$rand`;  #tianyulu
my $RM_test=`${repeatmasker}RepeatMasker -e ncbi -q -pa 1 -no_is -nolow dummy060817.fa.$rand -lib dummy060817.fa.$rand 2>/dev/null`;
`rm dummy060817.fa.$rand*`;
die "Error: The RMblast engine is not installed in RepeatMasker!\n" unless $RM_test=~s/done//gi;
# RepeatModeler
chomp ($repeatmodeler=`command -v RepeatModeler 2>/dev/null`) if $repeatmodeler eq '';
$repeatmodeler =~ s/\s+$//;
$repeatmodeler = dirname($repeatmodeler) unless -d $repeatmodeler;
$repeatmodeler="$repeatmodeler/" if $repeatmodeler ne '' and $repeatmodeler !~ /\/$/;
die "Error: RepeatModeler is not found in the RepeatModeler path $repeatmodeler!\n" unless -X "${repeatmodeler}RepeatModeler";
# AnnoSINE
chomp ($annosine=`command -v AnnoSINE_v2 2>/dev/null`) if $annosine eq '';
$annosine =~ s/\s+$//;
$annosine = dirname($annosine) unless -d $annosine;
$annosine="$annosine/" if $annosine ne '' and $annosine !~ /\/$/;
if (!-X "${annosine}AnnoSINE_v2"){
	if ($module_active{sine}){
		print STDERR "Warning: AnnoSINE is not found in the AnnoSINE path $annosine!\n\t\tThe SINE module will be skipped.\n\n";
		$skip_module{sine} = 1;
		}
	}
# LTR_retriever
chomp ($LTR_retriever=`command -v LTR_retriever 2>/dev/null`) if $LTR_retriever eq '';
$LTR_retriever =~ s/\s+$//;
$LTR_retriever = dirname($LTR_retriever) unless -d $LTR_retriever;
$LTR_retriever="$LTR_retriever/" if $LTR_retriever ne '' and $LTR_retriever !~ /\/$/;
die "Error: LTR_retriever is not found in the LTR_retriever path $LTR_retriever!\n" unless -X "${LTR_retriever}LTR_retriever";
# TEsorter
chomp ($TEsorter=`command -v TEsorter 2>/dev/null`) if $TEsorter eq '';
$TEsorter =~ s/\s+$//;
$TEsorter = dirname($TEsorter) unless -d $TEsorter;
$TEsorter="$TEsorter/" if $TEsorter ne '' and $TEsorter !~ /\/$/;
die "Error: TEsorter is not found in the TEsorter path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
# makeblastdb, blastn
chomp ($blastplus=`command -v makeblastdb 2>/dev/null`) if $blastplus eq '';
$blastplus =~ s/\s+$//;
$blastplus = dirname($blastplus) unless -d $blastplus;
$blastplus="$blastplus/" if $blastplus ne '' and $blastplus !~ /\/$/;
die "Error: makeblastdb is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}makeblastdb";
die "Error: blastn is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";
# mdust
chomp ($mdust=`command -v mdust 2>/dev/null`) if $mdust eq '';
$mdust =~ s/\s+$//;
$mdust = dirname($mdust) unless -d $mdust;
$mdust = "$mdust/" if $mdust ne '' and $mdust !~ /\/$/;
die "Error: mdust is not found in the mdust path $mdust!\n" unless -X "${mdust}mdust";
# trf
chomp ($trf=`command -v trf 2>/dev/null`) if $trf eq '';
$trf=~s/\n$//;
die "Error: Tandem Repeat Finder is not found in the TRF path $trf!\n" if $trf eq '' or !-X $trf;
# GRF
chomp ($GRF = `command -v grf-main 2>/dev/null`) if $GRF eq '';
$GRF =~ s/\n$//;
my $grfp= dirname ($GRF);
$grfp =~ s/\n$//;
$grfp="$grfp/" if $grfp ne '' and $grfp !~ /\/$/;
if ($GRF eq '' or !-X "${grfp}grf-main"){
	if ($module_active{tir}){
		print STDERR "Warning: The Generic Repeat Finder (GRF) is not found in the GRF path: $grfp!\n\t\tThe TIR module will be skipped.\n\n";
		$skip_module{tir} = 1;
		}
	}

# TIR-Learner  #tianyuLu
# Remove any trailing whitespace
$TIR_Learner =~ s/\s+$//;
if ($TIR_Learner eq "") {
	# Find TIR-Learner path and remove any trailing newline
	chomp ($TIR_Learner=`command -v TIR-Learner 2>/dev/null`);
} else {
	# # Extract directory name from path if path is not a directory
	# If path is directory
	if (-d $TIR_Learner) {
		# Add trailing slash if path not already end with slash
		$TIR_Learner .= "/" if $TIR_Learner !~ /\/$/;
		$TIR_Learner = "python3 $TIR_Learner/TIR-Learner.py";
	}
}
my $TIR_Learner_exe = (split /\s+/, $TIR_Learner)[-1];
chomp ($TIR_Learner_exe = `command -v $TIR_Learner_exe 2>/dev/null`) if defined $TIR_Learner_exe and $TIR_Learner_exe ne '' and $TIR_Learner_exe !~ /\//;
if ($TIR_Learner eq "" or !-X $TIR_Learner_exe){
	if ($module_active{tir}){
		print STDERR "Warning: TIR-Learner is not found in the path $TIR_Learner!\n\t\tThe TIR module will be skipped.\n\n";
		$skip_module{tir} = 1;
		}
	}


# LTR_FINDER_parallel  #tianyuLu
$LTR_FINDER =~ s/\s+$//;
if ($LTR_FINDER eq "") {
    chomp ($LTR_FINDER=`command -v LTR_FINDER_parallel 2>/dev/null`);
} else {
    if (-d $LTR_FINDER) {
        $LTR_FINDER .= "/" if $LTR_FINDER !~ /\/$/;
        $LTR_FINDER = "perl $LTR_FINDER/LTR_FINDER_parallel";
    }
}
my $LTR_FINDER_exe = (split /\s+/, $LTR_FINDER)[-1];
chomp ($LTR_FINDER_exe = `command -v $LTR_FINDER_exe 2>/dev/null`) if defined $LTR_FINDER_exe and $LTR_FINDER_exe ne '' and $LTR_FINDER_exe !~ /\//;
die "Error: LTR_FINDER_parallel is not found in the path $LTR_FINDER!\n" if $module_active{ltr} and ($LTR_FINDER eq "" or !-X $LTR_FINDER_exe);

# LTR_HARVEST_parallel  #tianyuLu
$LTR_HARVEST =~ s/\s+$//;
if ($LTR_HARVEST eq "") {
    chomp ($LTR_HARVEST=`command -v LTR_HARVEST_parallel 2>/dev/null`);
} else {
    if (-d $LTR_HARVEST) {
        $LTR_HARVEST .= "/" if $LTR_HARVEST !~ /\/$/;
        $LTR_HARVEST = "perl $LTR_HARVEST/LTR_HARVEST_parallel";
    }
}
my $LTR_HARVEST_exe = (split /\s+/, $LTR_HARVEST)[-1];
chomp ($LTR_HARVEST_exe = `command -v $LTR_HARVEST_exe 2>/dev/null`) if defined $LTR_HARVEST_exe and $LTR_HARVEST_exe ne '' and $LTR_HARVEST_exe !~ /\//;
die "Error: LTR_HARVEST_parallel is not found in the path $LTR_HARVEST!\n" if $module_active{ltr} and ($LTR_HARVEST eq "" or !-X $LTR_HARVEST_exe);

# HelitronScanner  #tianyuLu
$HelitronScanner =~ s/\s+$//;
if ($HelitronScanner eq "") {
    chomp ($HelitronScanner=`command -v HelitronScanner 2>/dev/null`);
} else {
    if (-d $HelitronScanner) {
        $HelitronScanner .= "/" if $HelitronScanner !~ /\/$/;
    }
}
if ($HelitronScanner eq "" or !(-d $HelitronScanner or -X $HelitronScanner)){
	if ($module_active{helitron}){
		print STDERR "Warning: HelitronScanner is not found in the path $HelitronScanner!\n\t\tThe Helitron module will be skipped.\n\n";
		$skip_module{helitron} = 1;
		}
	}

# make a softlink to the genome
my $genome_file = basename($genome);
unlink $genome_file if -l $genome_file;
if (-e $genome_file){
	die "Error: $genome_file already exists in the working directory and is not the input genome! Please rename or remove it and try again.\n" unless abs_path($genome_file) eq abs_path($genome);
	} else {
	`ln -s $genome $genome_file`;
	die "Error: failed to create the genome softlink $genome_file!\n" unless -l $genome_file;
	}
$genome = $genome_file;

# check $RMlib
if ($RMlib ne 'null'){
	if (-e $RMlib){
		print "\tA RepeatModeler library $RMlib is provided via --rmlib. Please make sure this is a RepeatModeler2 generated and classified library (some levels of unknown classification is OK).\n\n";
		chomp ($RMlib = `realpath $RMlib`);
		`ln -s $RMlib $genome.RM2.raw.fa` unless -s "$genome.RM2.raw.fa";
		$RMlib = "$genome.RM2.raw.fa";
		} else {
		die "\tERROR: The RepeatModeler library $RMlib you specified is not found!\n\n";
		}
	}

# check if duplicated sequences found
# EDTA.pl sets EDTA_DUP_CHECK_DONE=1 once it has checked the original genome for
# duplicate IDs; standalone runs (env unset) always perform the check here
my ($raw_id, $old_id);
unless ($ENV{EDTA_DUP_CHECK_DONE}){
	chomp (my $id_counts = `grep -a \\> $genome|sort|uniq -c|awk '{n++; s+=\$1} END{print s"\\t"n}'`);
	($raw_id, $old_id) = $id_counts =~ /^([0-9]+)\t([0-9]+)$/;
	die "Error: failed to read sequence IDs from $genome!\n" unless defined $raw_id;
	if ($raw_id > $old_id){
		chomp ($date = `date`);
		die "$date\tERROR: Identical sequence IDs found in the provided genome! Please resolve this issue and try again.\n";
		}
	}

if ($convert_name == 1){
if (-s "$genome.mod"){
	$genome = "$genome.mod";
	chomp ($date = `date`);
	print "$date\tExisting normalized genome $genome found, will use it.\n\n";
	} else {
	# Normalize genome: clean IDs, encode if needed
	chomp ($date = `date`);
	print "$date\tCleaning and normalizing sequence IDs...\n";
	`perl $seqid_codec encode_fasta $genome $genome.mod.seqid.map $genome.mod`;
	if ($? != 0){
		die "$date\tERROR: Genome normalization failed. Check error messages above.\n";
		}
	$genome = "$genome.mod";
	if (-s "$genome.seqid.map"){
		print "\tSequence IDs encoded. Mapping file: $genome.seqid.map\n\n";
		} else {
		print "\tSequence IDs are short enough, no encoding needed.\n\n";
		}
	# Verify unique ID count
	my $new_id = `grep -a \\> $genome|sort -u|wc -l`;
	chomp $new_id;
	# $old_id is not computed above when the original-genome check was skipped; recover it so this check always runs
	chomp ($old_id = `grep -a \\> $genome_file|sort -u|wc -l`) unless defined $old_id;
	chomp $old_id;
	if ($old_id != $new_id){
		chomp ($date = `date`);
		die "$date\tERROR: Seq ID normalization produced non-unique IDs. Please check your genome file.\n";
		}
	}
}
# pre-set parameters
my $genome_file_real_path=File::Spec->rel2abs($genome); # the genome file with real path

# Make working directories
`mkdir $genome.EDTA.raw` unless -e "$genome.EDTA.raw" && -d "$genome.EDTA.raw";
`mkdir $genome.EDTA.raw/LTR` unless -e "$genome.EDTA.raw/LTR" && -d "$genome.EDTA.raw/LTR";
`mkdir $genome.EDTA.raw/SINE` unless -e "$genome.EDTA.raw/SINE" && -d "$genome.EDTA.raw/SINE";
`mkdir $genome.EDTA.raw/LINE` unless -e "$genome.EDTA.raw/LINE" && -d "$genome.EDTA.raw/LINE";
`mkdir $genome.EDTA.raw/TIR` unless -e "$genome.EDTA.raw/TIR" && -d "$genome.EDTA.raw/TIR";
`mkdir $genome.EDTA.raw/Helitron` unless -e "$genome.EDTA.raw/Helitron" && -d "$genome.EDTA.raw/Helitron";
foreach my $dir ("$genome.EDTA.raw", "$genome.EDTA.raw/LTR", "$genome.EDTA.raw/SINE", "$genome.EDTA.raw/LINE", "$genome.EDTA.raw/TIR", "$genome.EDTA.raw/Helitron"){
	die "Cannot create directory $dir: $!\n" unless -d $dir;
	}

# touch the expected (empty) result files of modules skipped due to missing dependencies
if ($skip_module{sine}){
	`touch $genome.EDTA.raw/$genome.SINE.raw.fa` unless -e "$genome.EDTA.raw/$genome.SINE.raw.fa";
	}
if ($skip_module{tir}){
	foreach my $ext ("fa", "gff3", "bed"){
		`touch $genome.EDTA.raw/$genome.TIR.intact.raw.$ext` unless -e "$genome.EDTA.raw/$genome.TIR.intact.raw.$ext";
		}
	}
if ($skip_module{helitron}){
	foreach my $ext ("fa", "gff3", "bed"){
		`touch $genome.EDTA.raw/$genome.Helitron.intact.raw.$ext` unless -e "$genome.EDTA.raw/$genome.Helitron.intact.raw.$ext";
		}
	}
# modules excluded via --type still need their expected (empty) outputs so the
# downstream EDTA.pl / EDTA_processK.pl stages can proceed with empty libraries
if (not $module_active{sine}){
	`touch $genome.EDTA.raw/$genome.SINE.raw.fa` unless -e "$genome.EDTA.raw/$genome.SINE.raw.fa";
	}
if (not $module_active{line}){
	`touch $genome.EDTA.raw/$genome.LINE.raw.fa` unless -e "$genome.EDTA.raw/$genome.LINE.raw.fa";
	`touch $genome.EDTA.raw/$genome.RM2.fa` unless -e "$genome.EDTA.raw/$genome.RM2.fa";
	}


##################################################
###### Dispatch the raw TE discovery modules ######
##################################################

my %module_sub = (
	ltr => \&run_LTR_module,
	sine => \&run_SINE_module,
	line => \&run_LINE_module,
	tir => \&run_TIR_module,
	helitron => \&run_HEL_module,
	);
my @active_modules = grep { $module_active{$_} and not $skip_module{$_} } qw/ltr sine line tir helitron/;

# relative expected durations (bigger weight = longer module): thread budgets are
# sized from these so parallel modules finish at roughly the same time.
my %module_weight = (ltr=>2.5, sine=>1.2, line=>4, tir=>1.5, helitron=>1);
if ($module_weights ne ''){
	foreach my $pair (split /,/, $module_weights){
		my ($mod, $w) = $pair =~ /^\s*(\w+)\s*=\s*([0-9.]+)\s*$/;
		die "Error: malformed --module_weights entry \"$pair\" (expected module=number).\n"
			unless defined $mod and defined $w and exists $module_weight{$mod} and $w > 0;
		$module_weight{$mod} = $w;
		}
	}

# fork one child per module in $mods with weighted thread budgets over $budget_threads,
# wait for all, and die if any failed (never lets the parent continue on module failure)
sub run_module_batch {
	my ($mods, $budget_threads) = @_;
	return unless scalar @$mods;
	my $sum_w = 0;
	$sum_w += $module_weight{$_} for @$mods;
	my %module_budget;
	my $used = 0;
	foreach my $mod (@$mods){
		my $raw = $budget_threads * $module_weight{$mod} / $sum_w;
		$module_budget{$mod} = int($raw);
		$module_budget{$mod} = 1 if $module_budget{$mod} < 1;
		$used += $module_budget{$mod};
		}
	if ($used < $budget_threads){
		# hand out remaining threads by largest fractional remainder, then by weight
		my @byfrac = sort {
			my $fa = $budget_threads*$module_weight{$a}/$sum_w; $fa -= int($fa);
			my $fb = $budget_threads*$module_weight{$b}/$sum_w; $fb -= int($fb);
			$fb <=> $fa or $module_weight{$b} <=> $module_weight{$a}
			} @$mods;
		my $i = 0;
		while ($used < $budget_threads){
			$module_budget{$byfrac[$i % scalar(@byfrac)]}++;
			$used++;
			$i++;
			}
		}
	chomp ($date = `date`);
	print STDERR "$date\tRunning ".scalar(@$mods)." module(s) (".join(", ", @$mods).") in parallel. Thread budgets: "
		.join(", ", map {"$_=$module_budget{$_}"} @$mods).".\n";
	print STDERR "\tTip: cores sit idle after short modules finish; a higher --threads improves utilization on a multi-core node.\n"
		if $budget_threads < 4 * scalar(@$mods);
	my %module_pid;
	foreach my $mod (@$mods){
		my $pid = fork();
		die "Error: cannot fork a child process for the $mod module: $!\n" unless defined $pid;
		if ($pid == 0){
			# child: run only this module then exit; only the parent continues below
			$module_threads = $module_budget{$mod};
			print STDERR "\tRunning the $mod module with $module_threads threads (PID $$).\n";
			eval { $module_sub{$mod}->() };
			if ($@){
				print STDERR "\tError: the $mod module failed: $@\n";
				exit 1;
				}
			exit 0;
			}
		$module_pid{$pid} = $mod;
		}
	print STDERR "\n";
	my $failed = 0;
	foreach my $pid (sort keys %module_pid){
		waitpid($pid, 0);
		my $module_status = $?; # capture before any other command overwrites $?
		if ($module_status != 0){
			chomp ($date = `date`);
			print STDERR "$date\tError: the $module_pid{$pid} module (PID $pid) exited with a nonzero status ($module_status)!\n";
			$failed = 1;
			}
		}
	die "Error: One or more raw TE modules failed. Please check the error messages above.\n\n" if $failed;
	}

if ($parallel_modules == 1 and scalar(@active_modules) > 1 and $threads >= scalar(@active_modules)){
	# mode 1: run all active modules concurrently, weighted thread budgets
	run_module_batch(\@active_modules, $threads);
	} elsif ($parallel_modules == 2 and scalar(@active_modules) > 1 and $threads >= scalar(@active_modules)){
	# mode 2 (staged): the light modules run first (weighted split of ALL threads, so they
	# finish quickly), then the single heaviest module (LINE/RepeatModeler by default
	# weights) runs alone with the FULL thread budget — its thread scaling is sublinear,
	# so maximizing its budget beats running it concurrently at a reduced share
	my ($heaviest) = sort { $module_weight{$b} <=> $module_weight{$a} or $a cmp $b } @active_modules;
	my @light = grep { $_ ne $heaviest } @active_modules;
	chomp ($date = `date`);
	print STDERR "$date\tStaged schedule: stage 1 = light modules (".join(", ", @light).
		"), stage 2 = the $heaviest module alone with all $threads threads.\n";
	run_module_batch(\@light, $threads) if scalar @light;
	chomp ($date = `date`);
	print STDERR "$date\tStage 1 finished. Starting the $heaviest module with the full $threads threads.\n";
	run_module_batch([$heaviest], $threads);
	} else {
	# run active modules one after another, each with the full thread budget
	$module_threads = $threads;
	foreach my $mod (@active_modules){
		$module_sub{$mod}->();
		}
	}


###########################
###### LTR_retriever ######
###########################

sub run_LTR_module {

chomp ($date = `date`);
print STDERR "$date\tStart to find LTR candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/LTR" or die "Cannot enter $genome.EDTA.raw/LTR: $!\n";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
chomp ($date = `date`);
if ($overwrite eq 0 and -s "$genome.LTR.raw.fa"){
	print STDERR "$date\tExisting result file $genome.LTR.raw.fa found!\n\t\tWill keep this file without rerunning this module.\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify LTR retrotransposon candidates from scratch.\n\n";

# run LTRharvest and LTR_FINDER_parallel concurrently (independent searches)
my $half_threads = int($module_threads/2);
$half_threads = 1 if $half_threads < 1;
my %ltr_searcher_pid;
if ($overwrite eq 0 and -s "$genome.harvest.combine.scn"){
	print STDERR "$date\tExisting raw result $genome.harvest.scn found!\n\t\tWill use this for further analyses.\n\n";
	} else {
	my $harvest_pid = fork();
	die "Error: cannot fork a child process for LTR_HARVEST_parallel: $!\n" unless defined $harvest_pid;
	if ($harvest_pid == 0){
		# `perl $LTR_HARVEST -seq $genome -threads $half_threads -gt $genometools -size 1000000 -time 300`;
		my $rc = system("$LTR_HARVEST -seq $genome -threads $half_threads -gt $genometools -size 1000000 -time 300 2>> $genome.harvest.log");  #tianyulu #575
		exit($rc == 0 ? 0 : 1);
		}
	$ltr_searcher_pid{$harvest_pid} = "LTR_HARVEST_parallel";
	}

# run LTR_FINDER_parallel
if ($overwrite eq 0 and -s "$genome.finder.combine.scn"){
	print STDERR "$date\tExisting raw result $genome.finder.combine.scn found!\n\t\tWill use this for further analyses.\n\n";
	} else {
	my $finder_pid = fork();
	die "Error: cannot fork a child process for LTR_FINDER_parallel: $!\n" unless defined $finder_pid;
	if ($finder_pid == 0){
		# `perl $LTR_FINDER -seq $genome -threads $half_threads -harvest_out -size 1000000 -time 300`;
		my $rc = system("$LTR_FINDER -seq $genome -threads $half_threads -harvest_out -size 1000000 -time 300 2>> $genome.finder.log");  #tianyulu #575
		exit($rc == 0 ? 0 : 1);
		}
	$ltr_searcher_pid{$finder_pid} = "LTR_FINDER_parallel";
	}
foreach my $pid (sort keys %ltr_searcher_pid){
	waitpid($pid, 0);
	die "Error: $ltr_searcher_pid{$pid} exited with a nonzero status ($?)! Check the log files in $genome.EDTA.raw/LTR.\n" if $? != 0;
	}

# run LTR_retriever
my $status = 0;
if ($overwrite eq 0 and (-s "$genome.mod.LTRlib.fa" or -s "$genome.LTRlib.fa")){
	print STDERR "$date\tExisting LTR_retriever result found!\n\t\tWill use this for further analyses.\n\n";
	} else {
	`cat $genome.harvest.combine.scn $genome.finder.combine.scn > $genome.rawLTR.scn`;
	my $we_flag = $wholeelement ? "-wholeelement" : "";
	$status = system("${LTR_retriever}LTR_retriever -genome $genome -inharvest $genome.rawLTR.scn -u $miu -threads $module_threads -noanno $we_flag -trf_path $trf -blastplus $blastplus -repeatmasker $repeatmasker -salvage 1 -convert_seq_name 0");
	}

# get full-length LTR from pass.list (or use LTR_retriever's whole-element output)
if ($wholeelement and -s "$genome.LTR.intact.fa"){
	# whole-element mode: LTR_retriever already produced annotated intact elements
	`cp $genome.LTR.intact.fa $genome.LTR.intact.fa.ori`;
} else {
	`awk '{if (\$1 !~ /#/) print \$1"\\t"\$1}' $genome.pass.list | perl $call_seq - -C $genome > $genome.LTR.intact.fa.ori`;
	`perl -i -nle 's/\\|.*//; print \$_' $genome.LTR.intact.fa.ori`;
	`perl $rename_LTR $genome.LTR.intact.fa.ori $genome.defalse > $genome.LTR.intact.fa.anno`;
	`mv $genome.LTR.intact.fa.anno $genome.LTR.intact.fa.ori`;
}

# remove simple repeats and candidates with simple repeats at terminals
`${mdust}mdust $genome.LTR.intact.fa.ori > $genome.LTR.intact.fa.ori.dusted`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.LTR.intact.fa.ori.dusted > $genome.LTR.intact.fa.ori.dusted.cln`;

# annotate and remove not LTR candidates
if (-s "$genome.LTR.intact.fa.ori.dusted.cln"){
	my $tesorter_err = `${TEsorter}TEsorter $genome.LTR.intact.fa.ori.dusted.cln --disable-pass2 -p $module_threads 2>&1`;
	die "TEsorter failed on $genome.LTR.intact.fa.ori.dusted.cln: $tesorter_err\n" unless -s "$genome.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv";
	`perl $cleanup_misclas $genome.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv`;
} elsif ($status == 0){
	print "\t\tLTR_retriever is finished without error, but no LTR is identified.\n\n";
	`touch $genome.LTR.intact.fa.ori.dusted.cln.cln`;
} else {
	print "\t\tLTR_retriever exited with error, please test run EDTA with EDTA/test/ to make sure the installation is correct.\n\n";
}
`mv $genome.LTR.intact.fa.ori.dusted.cln.cln $genome.LTR.intact.raw.fa`;
`mv $genome.LTR.intact.fa.ori.dusted.cln.cln.list $genome.LTR.intact.raw.fa.anno.list`;
`cp $genome.LTR.intact.raw.fa.anno.list ../`;

# generate annotated output and gff
`perl $output_by_list 1 $genome.LTR.intact.fa.ori 1 $genome.LTR.intact.raw.fa -FA -ex|grep -a \\>|perl -nle 's/>//; print "Name\\t\$_"' > $genome.LTR.intact.fa.ori.rmlist`;
`perl $filter_gff $genome.pass.list.gff3 $genome.LTR.intact.fa.ori.rmlist | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' > $genome.LTR.intact.raw.gff3`;
	}

# remove the short-ID genome symlink created above (line ~380). Done outside the
# else block so it also runs on resume (--overwrite 0), when that block is skipped
# and the symlink would otherwise be left behind. Nothing downstream reads it.
`rm -f $genome 2>/dev/null`;

# copy result files out
`touch $genome.LTRlib.fa` unless -e "$genome.LTRlib.fa";
`cp $genome.LTRlib.fa $genome.LTR.raw.fa`;
`cp $genome.LTRlib.fa ../$genome.LTR.raw.fa`;
`cp $genome.LTR.intact.fa ../$genome.LTR.intact.raw.fa` if -s "$genome.LTR.intact.fa";
`cp $genome.LTR.intact.gff3 ../$genome.LTR.intact.raw.gff3` if -s "$genome.LTR.intact.gff3";
`cp $genome.LTR.intact.raw.fa $genome.LTR.intact.raw.gff3 ../ 2>/dev/null`;
`cp $genome.LTRlib.fa.LTRbound ../$genome.LTRlib.fa.LTRbound 2>/dev/null` if $wholeelement;
chdir '../..' or die "Cannot return to the working directory from the LTR module: $!\n";

# check results
chomp ($date = `date`);
die "Error: LTR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.LTR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.LTR.raw.fa"){
	print STDERR "$date\tFinish finding LTR candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The LTR result file has 0 bp!\n\n";
	}

}


#############################
######    AnnoSINE     ######
#############################
sub run_SINE_module {

chomp ($date = `date`);
print STDERR "$date\tStart to find SINE candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/SINE" or die "Cannot enter $genome.EDTA.raw/SINE: $!\n";
`ln -s ../../$genome $genome` unless -s $genome;

# Remove existing results
`rm -rf Seed_SINE.fa Step* HMM_out 2>/dev/null` if $overwrite eq 1;

# run AnnoSINE_v2
my $status; # record status of AnnoSINE execution
if (-s "Seed_SINE.fa"){
	print STDERR "$date\tExisting result file Seed_SINE.fa found!\n\t\tWill keep this file without rerunning this module.\n\t\tPlease specify --overwrite 1 if you want to rerun AnnoSINE_v2.\n\n";
	} else { 
	$status = system("python3 ${annosine}AnnoSINE_v2 --temp_dir $genome_file_real_path.EDTA.raw/SINE/ -t $module_threads -a 2 --num_alignments 50000 -rpm 0 --copy_number 3 --shift 100 -auto 1 3 $genome ./ 2>> $genome.AnnoSINE.log");
	#$status = system("python3 ${annosine}AnnoSINE_v2 --temp_dir $genome_file_real_path.EDTA.raw/SINE/ -t $threads -a 2 --num_alignments 50000 -rpm 0 --copy_number 3 --shift 100 -auto 1 3 $genome ./ > /dev/null 2>&1");
	}

# filter and reclassify AnnoSINE candidates with TEsorter and make SINE library
if (-s "Seed_SINE.fa"){
	# annotate and remove non-SINE candidates
	`awk '{gsub(/Unknown/, "unknown"); print \$1}' Seed_SINE.fa > $genome.AnnoSINE.raw.fa`;
	if (-s "$genome.AnnoSINE.raw.fa"){
		`${TEsorter}TEsorter $genome.AnnoSINE.raw.fa --disable-pass2 -p $module_threads 2>/dev/null`;
		die "Error: TEsorter failed to generate a non-empty $genome.AnnoSINE.raw.fa.rexdb.cls.tsv for the SINE module!\n" unless -s "$genome.AnnoSINE.raw.fa.rexdb.cls.tsv";
		`perl $cleanup_misclas $genome.AnnoSINE.raw.fa.rexdb.cls.tsv`;

		# clean up tandem repeat
		`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.AnnoSINE.raw.fa.cln > $genome.SINE.raw.fa`;
		} else {
		print "\t\tAnnoSINE is finished without error, but no SINE candidate is obtained.\n\n";
		`touch $genome.SINE.raw.fa`;
		}
	}
elsif ($status == 0) {
	print "\t\tAnnoSINE is finished without error, but the Seed_SINE.fa file is not produced.\n\n";
       	`touch $genome.SINE.raw.fa`;
	}
else {
	print "\t\tAnnoSINE exited with error, please test run AnnoSINE to make sure it's working.\n\n";
	}

# copy result files out
`cp $genome.SINE.raw.fa ../`;
chdir '../..' or die "Cannot return to the working directory from the SINE module: $!\n";

# check results
chomp ($date = `date`);
die "Error: SINE results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.SINE.raw.fa";
if (-s "$genome.EDTA.raw/$genome.SINE.raw.fa"){
	print STDERR "$date\tFinish finding SINE candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The SINE result file has 0 bp!\n\n";
	}

}


#############################
######  RepeatModeler  ######
#############################

sub run_LINE_module {

chomp ($date = `date`);
print STDERR "$date\tStart to find LINE candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/LINE" or die "Cannot enter $genome.EDTA.raw/LINE: $!\n";
`ln -s ../../$genome $genome` unless -s $genome;
`cp ../../$RMlib $RMlib` if $RMlib ne 'null';

# Try to recover existing results or run RepeatModeler2
chomp ($date = `date`);
if ($overwrite eq 0 and -s $RMlib){
	if (-s "$genome-families.fa"){
		print STDERR "$date\tExisting result file $genome-families.fa found!\n\t\tWill not use the provided RepeatModeler2 library since --overwrite 0.\n\t\tPlease specify --overwrite 1 if you want to use the provided --rmlib file.\n\n";
		} else {
		`cp $RMlib "$genome-families.fa" 2>/dev/null`;
		}
	}

if ($overwrite eq 0 and -s "$genome-families.fa"){
	print STDERR "$date\tExisting result file $genome-families.fa found!\n\t\tWill keep this file without rerunning this module.\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	# run RepeatModeler2
	print STDERR "$date\tIdentify LINE retrotransposon candidates from scratch.\n\n";
	my $status; # record status of RepeatModeler execution
	`${repeatmodeler}BuildDatabase -name $genome $genome`;
	$status = system("${repeatmodeler}RepeatModeler -engine ncbi -threads $module_threads -database $genome  > repeatmodeler.log 2>&1");
	if ($status != 0) {
		# Execute the old version of RepeatModeler
		warn "RepeatModeler failed with -threads, retrying with -pa...\n";
		$status = system("${repeatmodeler}RepeatModeler -engine ncbi -pa $module_threads -database $genome > repeatmodeler.log 2>&1");
		if ($status != 0) {
			print "ERROR: RepeatModeler did not run correctly. Please test run this command:
				${repeatmodeler}RepeatModeler -engine ncbi -pa $module_threads -database $genome
				ERROR\n";
			exit;
			}
		}
	`rm $genome*nal $genome*nhr $genome*nin $genome*nnd $genome*nni $genome*nog $genome*nsq $genome*njs $genome*translation 2>/dev/null`;
	}

# filter and reclassify RepeatModeler candidates with TEsorter and make LINE library
if (-s "$genome-families.fa"){
	# annotate and remove misclassified candidates
	`awk '{gsub(/Unknown/, "unknown"); print \$1}' $genome-families.fa > $genome.RM2.raw.fa` if -e "$genome-families.fa";
		my $tesorter_err = `${TEsorter}TEsorter $genome.RM2.raw.fa --disable-pass2 -p $module_threads 2>&1`;
		die "TEsorter failed on $genome.RM2.raw.fa: $tesorter_err\n" if -s "$genome.RM2.raw.fa" and not -s "$genome.RM2.raw.fa.rexdb.cls.tsv";
	`perl $cleanup_misclas $genome.RM2.raw.fa.rexdb.cls.tsv`;

	# reclassify clean candidates
	$tesorter_err = `${TEsorter}TEsorter $genome.RM2.raw.fa.cln --disable-pass2 -p $module_threads 2>&1`;
	die "TEsorter failed on $genome.RM2.raw.fa.cln: $tesorter_err\n" if -s "$genome.RM2.raw.fa.cln" and not -s "$genome.RM2.raw.fa.cln.rexdb.cls.tsv";
	`perl -nle 's/>\\S+\\s+/>/; print \$_' $genome.RM2.raw.fa.cln.rexdb.cls.lib > $genome.RM2.raw.fa.cln`;

        # clean up tandem repeat
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.RM2.raw.fa.cln > $genome.RM2.raw.fa.cln2`;
	`grep -P 'LINE|SINE' $genome.RM2.raw.fa.cln2 | perl $output_by_list 1 $genome.RM2.raw.fa.cln2 1 - -FA > $genome.LINE.raw.fa`;
	`grep -P 'LINE|SINE' $genome.RM2.raw.fa.cln2 | perl $output_by_list 1 $genome.RM2.raw.fa.cln2 1 - -FA -ex > $genome.RM2.fa`;
	} else {
	print "\t\tRepeatModeler is finished, but the $genome-families.fa file is not produced.\n\n";
	`touch $genome.RM2.raw.fa $genome.LINE.raw.fa $genome.RM2.fa`;
	}

# copy result files out
`cp $genome.LINE.raw.fa $genome.RM2.fa ../`; #update the filtered RM2 result in the EDTA/raw folder
`cp $genome.RM2.raw.fa ../../`; #update the raw RM2 result in the EDTA folder
chdir '../..' or die "Cannot return to the working directory from the LINE module: $!\n";

# check results
chomp ($date = `date`);
die "Error: LINE results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.LINE.raw.fa";
if (-s "$genome.EDTA.raw/$genome.LINE.raw.fa"){
	print STDERR "$date\tFinish finding LINE candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The LINE result file has 0 bp!\n\n";
	}
}


###########################
######  TIR-Learner  ######
###########################
sub run_TIR_module {

chomp ($date = `date`);
print STDERR "$date\tStart to find TIR candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/TIR" or die "Cannot enter $genome.EDTA.raw/TIR: $!\n";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
chomp ($date = `date`);
if ($overwrite eq 0 and (-s "$genome.TIR.intact.raw.fa" or -s "$genome.TIR.intact.fa")){
	print STDERR "$date\tExisting result file $genome.TIR.intact.raw.fa found!\n\t\tWill keep this file without rerunning this module.\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify TIR candidates from scratch.\n\n";
	print STDERR "Species: $species\n";

	# run TIR-Learner
	my $status = 1;
	if ($overwrite eq 0 and -s "./TIR-Learner-Result/TIR-Learner_FinalAnn.fa"){
		print STDERR "$date\tExisting raw result TIR-Learner_FinalAnn.fa found!\n\t\tWill use this for further analyses.\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
		} else {
		$status = system("$TIR_Learner -f $genome_file_real_path -s $species -p $module_threads -l $maxint -o $genome_file_real_path.EDTA.raw/TIR");
		}

	# clean raw predictions with flanking alignment
	`perl $rename_tirlearner ./TIR-Learner-Result/TIR-Learner_FinalAnn.fa | perl -nle 's/TIR-Learner_//gi; print \$_' > $genome.TIR`;
	`perl $get_ext_seq $genome $genome.TIR`;
	`perl $flank_filter -genome $genome -query $genome.TIR.ext30.fa -miniden 90 -mincov 0.9 -maxct 20 -blastplus $blastplus -t $module_threads -overwrite $overwrite`;

	# recover superfamily info
	`perl $output_by_list 1 $genome.TIR 1 $genome.TIR.ext30.fa.pass.fa -FA -MSU0 -MSU1 > $genome.TIR.ext30.fa.pass.fa.ori`;

	# remove simple repeats and candidates with simple repeats at terminals
	`${mdust}mdust $genome.TIR.ext30.fa.pass.fa.ori > $genome.TIR.ext30.fa.pass.fa.dusted`;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.TIR.ext30.fa.pass.fa.dusted > $genome.TIR.ext30.fa.pass.fa.dusted.cln`;

	# annotate and remove non-TIR candidates
	if (-s "$genome.TIR.ext30.fa.pass.fa.dusted.cln"){
		my $tesorter_err = `${TEsorter}TEsorter $genome.TIR.ext30.fa.pass.fa.dusted.cln --disable-pass2 -p $module_threads 2>&1`;
		die "TEsorter failed on $genome.TIR.ext30.fa.pass.fa.dusted.cln: $tesorter_err\n" unless -s "$genome.TIR.ext30.fa.pass.fa.dusted.cln.rexdb.cls.tsv";
		`perl $cleanup_misclas $genome.TIR.ext30.fa.pass.fa.dusted.cln.rexdb.cls.tsv`;
	} elsif ($status == 0) {
		print "\t\tTIR-Learner is finished without error, but no TIR is identified.\n\n";
		`touch $genome.TIR.ext30.fa.pass.fa.dusted.cln.cln`;
	} else {
		print "\t\tTIR-Learner exited with error, please test run EDTA with EDTA/test/ to make sure the installation is correct.\n\n";
	}
	`mv $genome.TIR.ext30.fa.pass.fa.dusted.cln.cln $genome.TIR.intact.raw.fa`;
	`cp $genome.TIR.ext30.fa.pass.fa.dusted.cln.cln.list $genome.TIR.intact.raw.fa.anno.list`;
	`cp $genome.TIR.intact.raw.fa.anno.list ../`;

	# get gff3 of intact TIR elements
	`perl -nle 's/\\-\\+\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class='DNA'; \$class='MITE' if \$e-\$s+1 <= 600; my (\$tir, \$iden, \$tsd)=(\$1, \$2/100, \$3) if \$anno=~/TIR:(.*)_([0-9.]+)_TSD:([a-z0-9._]+)_LEN/i; print "\$chr \$s \$e \$chr:\$s..\$e \$class/\$supfam structural \$iden . . . TSD=\$tsd;TIR=\$tir"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3 | perl $output_by_list 4 - 1 $genome.TIR.intact.raw.fa -MSU0 -MSU1 > $genome.TIR.intact.raw.bed`;
	`perl $bed2gff $genome.TIR.intact.raw.bed TIR > $genome.TIR.intact.raw.gff3`;
	}

# copy result files out
`cp $genome.TIR.intact.bed ../$genome.TIR.intact.raw.bed` if -s "$genome.TIR.intact.bed"; # recover <EDTA2.2 results
`cp $genome.TIR.intact.gff3 ../$genome.TIR.intact.raw.gff3` if -s "$genome.TIR.intact.gff3";
`cp $genome.TIR.intact.fa ../$genome.TIR.intact.raw.fa` if -s "$genome.TIR.intact.fa";
`cp $genome.TIR.intact.raw.fa $genome.TIR.intact.raw.gff3 $genome.TIR.intact.raw.bed ../ 2>/dev/null`;
chdir '../..' or die "Cannot return to the working directory from the TIR module: $!\n";

# check results
chomp ($date = `date`);
die "Error: TIR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.TIR.intact.raw.fa";
if (-s "$genome.EDTA.raw/$genome.TIR.intact.raw.fa"){
	print STDERR "$date\tFinish finding TIR candidates.\n\n";
	} else {
	print STDERR "Warning: The TIR result file has 0 bp!\n\n";
	}

}


#############################
###### HelitronScanner ######
#############################
sub run_HEL_module {

chomp ($date = `date`);
print STDERR "$date\tStart to find Helitron candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/Helitron" or die "Cannot enter $genome.EDTA.raw/Helitron: $!\n";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
chomp ($date = `date`);
if ($overwrite eq 0 and (-s "$genome.Helitron.intact.raw.fa" or -s "$genome.Helitron.intact.fa")){
	print STDERR "$date\tExisting result file $genome.Helitron.intact.raw.fa found!\n\t\tWill keep this file without rerunning this module.\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify Helitron candidates from scratch.\n\n";

# run HelitronScanner
my $status = 1;
if ($overwrite eq 0 and (-s "$genome.HelitronScanner.draw.hel.fa" and -s "$genome.HelitronScanner.draw.rc.hel.fa")){
#cat $genome.HelitronScanner.draw.hel.fa $genome.HelitronScanner.draw.rc.hel.fa
	print STDERR "$date\tExisting HelitronScanner result files $genome.HelitronScanner.draw.hel.fa $genome.HelitronScanner.draw.rc.hel.fa found!\n\t\tWill keep these files without rerunning HelitronScanner\n\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	$status = system("python3 $HelitronScanner_Runner --genome $genome --cpu $module_threads --hsdir \"$HelitronScanner\" 2>> $genome.HelitronScanner.log");
	}

# filter candidates based on repeatness of flanking regions
`perl $format_helitronscanner -genome $genome -sitefilter 1 -minscore 12 -keepshorter 1 -extlen 30 -extout 0`;
`perl $format_helitronscanner -genome $genome -sitefilter 1 -minscore 12 -keepshorter 1 -extlen 30 -extout 1`;
`perl $flank_filter -genome $genome -query $genome.HelitronScanner.filtered.ext.fa -miniden 90 -mincov 0.9 -maxct 5 -blastplus $blastplus -t $module_threads -overwrite $overwrite`; #more relaxed
#`perl $flank_filter -genome $genome -query $genome.HelitronScanner.filtered.ext.fa -miniden 80 -mincov 0.8 -maxct 5 -blastplus $blastplus -t $threads -overwrite $overwrite`; #more stringent

# remove simple repeats and candidates with simple repeats at terminals
`perl $output_by_list 1 $genome.HelitronScanner.filtered.fa 1 $genome.HelitronScanner.filtered.ext.fa.pass.fa -FA > $genome.HelitronScanner.filtered.fa.pass.fa`;
`${mdust}mdust $genome.HelitronScanner.filtered.fa.pass.fa > $genome.HelitronScanner.filtered.fa.pass.fa.dusted`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.HelitronScanner.filtered.fa.pass.fa.dusted | perl -nle 's/^(>.*)\\s+(.*)\$/\$1#DNA\\/Helitron\\t\$2/; print \$_' > $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln`;

# annotate and remove non-Helitron candidates
if (-s "$genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln"){
	my $tesorter_err = `${TEsorter}TEsorter $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln --disable-pass2 -p $module_threads 2>&1`;
	die "TEsorter failed on $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln: $tesorter_err\n" unless -s "$genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.rexdb.cls.tsv";
	`perl $cleanup_misclas $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.rexdb.cls.tsv`;
} elsif ($status == 0) {
	print "\t\tHelitronScanner is finished without error, but no Helitron is identified.\n\n";
	`touch $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln`;
} else {
	print "\t\tHelitronScanner exited with error, please test run EDTA with EDTA/test/ to make sure the installation is correct.\n\n";
}
`mv $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln $genome.Helitron.intact.raw.fa`;
`cp $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln.list $genome.Helitron.intact.raw.fa.anno.list`;
`cp $genome.Helitron.intact.raw.fa.anno.list ../`;

# get intact Helitrons and gff3
`perl $make_bed $genome.Helitron.intact.raw.fa > $genome.Helitron.intact.raw.bed`;
`perl $bed2gff $genome.Helitron.intact.raw.bed HEL > $genome.Helitron.intact.raw.gff3`;
	}

# copy result files out
`cp $genome.Helitron.intact.bed ../$genome.Helitron.intact.raw.bed` if -s "$genome.Helitron.intact.bed"; # recover <EDTA2.2 results
`cp $genome.Helitron.intact.gff3 ../$genome.Helitron.intact.raw.gff3` if -s "$genome.Helitron.intact.gff3";
`cp $genome.Helitron.intact.fa ../$genome.Helitron.intact.raw.fa` if -s "$genome.Helitron.intact.fa";
`cp $genome.Helitron.intact.raw.fa $genome.Helitron.intact.raw.gff3 $genome.Helitron.intact.raw.bed ../ 2>/dev/null`;
chdir '../..' or die "Cannot return to the working directory from the Helitron module: $!\n";

# check results
chomp ($date = `date`);
die "Error: Helitron results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.Helitron.intact.raw.fa";
if (-s "$genome.EDTA.raw/$genome.Helitron.intact.raw.fa"){
	print STDERR "$date\tFinish finding Helitron candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The Helitron result file has 0 bp!\n\n";
	}

}

chomp ($date = `date`);
print STDERR "$date\tExecution of EDTA_raw.pl is finished!\n\n";

# clean up the run-private scratch dir (only ever created for standalone runs;
# under EDTA.pl the inherited TMPDIR belongs to the parent)
if ($raw_own_tmp and -d $ENV{TMPDIR}){
	rmtree($ENV{TMPDIR}, { error => \my $err } );
	}
