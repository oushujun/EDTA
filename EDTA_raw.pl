#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

########################################################
##### Perform initial searches for TE candidates    ####
##### Shujun Ou (shujun.ou.1@gmail.com, 05/31/2019) ####
########################################################

## Input:
#	$genome

## Output:
#	$genome.LTR.raw.fa
#	$genome.TIR.raw.fa
#	$genome.MITE.raw.fa
#	$genome.Helitron.raw.fa

my $usage = "\nObtain raw TE libraries using various structure-based programs
	perl EDTA_raw.pl [options]
		-genome	[File]	The genome FASTA
		-species [Rice|Maize|others]	Specify the species for identification of TIR candidates. Default: others
		-type	[ltr|tir|mite|helitron|all]	Specify which type of raw TE candidates you want to get. Default: all
		-blastplus      [path]  Path to the blastn program. Defalut: read from \$ENV
		-mdust	[program]	The mdust program. Default: included in this package.
		-overwrite	[0|1]	If previous results are found, decide to overwrite (1, rerun) or not (0, default).
		-threads|-t	[int]	Number of theads to run this script. Default: 4
		-help|-h	Display this help info
\n";

# pre-defined
my $genome = '';
my $species = 'others';
my $type = 'all';
my $overwrite = 0; #0, no rerun. 1, rerun even old results exist.
my $maxint = 5000; #maximum interval length (bp) between TIRs (for GRF in TIR-Learner)
my $threads = 4;
my $script_path = $FindBin::Bin;
my $genometools = "$script_path/bin/genometools-1.5.10/bin/gt";
my $LTR_FINDER = "$script_path/bin/LTR_FINDER_parallel/LTR_FINDER_parallel";
my $LTR_retriever = "$script_path/bin/LTR_retriever/LTR_retriever";
my $TIR_Learner = "$script_path/bin/TIR-Learner1.23/TIR-Learner2.sh";
my $rename_tirlearner = "$script_path/util/rename_tirlearner.pl";
my $MITE_Hunter = "$script_path/bin/MITE-Hunter2/MITE_Hunter_manager.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $get_ext_seq = "$script_path/util/get_ext_seq.pl";
my $HelitronScanner = "$script_path/util/run_helitron_scanner.sh";
my $format_helitronscanner = "$script_path/util/format_helitronscanner_out.pl";
my $flank_filter = "$script_path/util/flanking_filter.pl";
my $mdust = "";
my $blastplus = ''; #path to the blastn program

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$species = $ARGV[$k+1] if /^-species$/i and $ARGV[$k+1] !~ /^-/;
	$type = lc $ARGV[$k+1] if /^-type$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^-blastplus$/i and $ARGV[$k+1] !~ /^-/;
	$mdust = $ARGV[$k+1] if /^-mdust$/i and $ARGV[$k+1] !~ /^-/;
	$overwrite = $ARGV[$k+1] if /^-overwrite$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
  }

my $date=`date`;
chomp ($date);
print STDERR "$date\tEDTA_raw: Check files and dependencies, prepare working directories.\n\n";

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "The GenomeTools is not found in $genometools!\n" unless -s $genometools;
die "The LTR_FINDER_parallel is not found in $LTR_FINDER!\n" unless -s $LTR_FINDER;
die "The LTR_retriever is not found in $LTR_retriever!\n" unless -s $LTR_retriever;
die "The TIR_Learner is not found in $TIR_Learner!\n" unless -s $TIR_Learner;
die "The MITE_Hunter is not found in $MITE_Hunter!\n" unless -s $MITE_Hunter;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script get_ext_seq.pl is not found in $get_ext_seq!\n" unless -s $get_ext_seq;
die "The HelitronScanner is not found in $HelitronScanner!\n" unless -s $HelitronScanner;
die "The script format_helitronscanner_out.pl is not found in $format_helitronscanner!\n" unless -s $format_helitronscanner;
die "The script flanking_filter.pl is not found in $flank_filter!\n" unless -s $flank_filter;
$blastplus=`which blastn 2>/dev/null` if $blastplus eq '';
$blastplus=~s/blastn\n//;
die "makeblastdb is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}makeblastdb";
die "blastn is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;

# Make working directories
`mkdir $genome.EDTA.raw` unless -e "$genome.EDTA.raw" && -d "$genome.EDTA.raw";
`mkdir $genome.EDTA.raw/LTR` unless -e "$genome.EDTA.raw/LTR" && -d "$genome.EDTA.raw/LTR";
`mkdir $genome.EDTA.raw/TIR` unless -e "$genome.EDTA.raw/TIR" && -d "$genome.EDTA.raw/TIR";
`mkdir $genome.EDTA.raw/MITE` unless -e "$genome.EDTA.raw/MITE" && -d "$genome.EDTA.raw/MITE";
`mkdir $genome.EDTA.raw/Helitron` unless -e "$genome.EDTA.raw/Helitron" && -d "$genome.EDTA.raw/Helitron";


###########################
###### LTR_retriever ######
###########################

if ($type eq "ltr" or $type eq "all"){

$date=`date`;
chomp ($date);
print STDERR "$date\tStart to find LTR candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/LTR";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
$date=`date`;
chomp ($date);
if ($overwrite eq 0 and -s "$genome.LTR.raw.fa"){
	print STDERR "$date\tExisting result file $genome.LTRlib.fa found! Will keep this file without rerunning this module.\n\tPlease specify -overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify LTR retrotransposon candidates from scratch.\n\n";

# run LTRharvest
`$genometools suffixerator -db $genome -indexname $genome -tis -suf -lcp -des -ssp -sds -dna -mirrored 2>/dev/null`;
`$genometools ltrharvest -index $genome -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $genome.harvest.scn 2>/dev/null`;
`rm $genome.des $genome.esq $genome.lcp $genome.llv $genome.md5 $genome.prj $genome.sds $genome.ssp $genome.suf 2>/dev/null`;

# run LTR_FINDER_parallel
`perl $LTR_FINDER -seq $genome -threads $threads -harvest_out -size 1000000 -time 300`;

# run LTR_retriever
`cat $genome.harvest.scn $genome.finder.combine.scn > $genome.rawLTR.scn`;
`perl $LTR_retriever -genome $genome -inharvest $genome.rawLTR.scn -threads $threads -noanno`;
`for i in *.mod.*; do mv \$i \$(echo \$i|sed 's/\\.mod\\././'); done > /dev/null 2>&1`;
	}

# copy result files out
`cp $genome.LTRlib.fa $genome.LTR.raw.fa`;
`cp $genome.LTRlib.fa ../$genome.LTR.raw.fa`;
chdir '../..';

# check results
$date=`date`;
chomp ($date);
die "Error: LTR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.LTR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.LTR.raw.fa"){
	print STDERR "$date\tFinish finding LTR candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The LTR result file has 0 bp!\n\n";
	}

}

###########################
######  TIR-Learner  ######
###########################

if ($type eq "tir" or $type eq "all"){

$date=`date`;
chomp ($date);
print STDERR "$date\tStart to find TIR candidates.\n\n";

# enter the working directory, filter out short sequences and create genome softlink
chdir "$genome.EDTA.raw/TIR";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
$date=`date`;
chomp ($date);
if ($overwrite eq 0 and -s "$genome.TIR.raw.fa"){
	print STDERR "$date\tExisting result file $genome.TIR.raw.fa found! Will keep this file without rerunning this module.\n\tPlease specify -overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify TIR candidates from scratch.\n\n";

	$species =~ s/rice/Rice/i;
	$species =~ s/maize/Maize/i;
	$species =~ s/others/others/i;
	print STDERR "Species: $species\n";

	# run TIR-Learner
	`sh $TIR_Learner -g $genome -s $species -t $threads -l $maxint`;
	`perl $rename_tirlearner ./TIR-Learner-Result/TIR-Learner_FinalAnn.fa | perl -nle 's/TIR-Learner_//g; print \$_' > $genome.TIR`;

	# clean raw predictions with flanking alignment
	`perl $get_ext_seq $genome $genome.TIR`;
	`perl $flank_filter -genome $genome -query $genome.TIR.ext30.fa -miniden 90 -mincov 0.9 -maxct 20 -blastplus $blastplus -t $threads`;

	# recover superfamily info
	`perl  $output_by_list 1 $genome.TIR 1 $genome.TIR.ext30.fa.pass.fa -FA -MSU0 -MSU1 > $genome.TIR.ext30.fa.pass.fa.ori`;

	# remove simple repeats and candidates with simple repeats at terminals
	`${mdust}mdust $genome.TIR.ext30.fa.pass.fa.ori > $genome.TIR.ext30.fa.pass.fa.dusted`;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.ext30.fa.pass.fa.dusted > $genome.TIR.ext30.fa.pass.fa.dusted.cln`;
	`cp $genome.TIR.ext30.fa.pass.fa.dusted.cln $genome.TIR.raw.fa`;

	}

# copy result files out
`cp $genome.TIR.raw.fa ../$genome.TIR.raw.fa`;
chdir '../..';

# check results
die "Error: TIR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.TIR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.TIR.raw.fa"){
	print STDERR "Finish finding TIR candidates.\n\n";
	} else {
	`touch "$genome.EDTA.raw/$genome.TIR.raw.fa"`;
	print STDERR "Warning: The TIR result file has 0 bp!\n\n";
	}

}


###########################
######  MITE-Hunter  ######
###########################

if ($type eq "mite" or $type eq "all"){

$date=`date`;
chomp ($date);
print STDERR "$date\tStart to find MITE candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/MITE";
`rm -rf genome*`;
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
$date=`date`;
chomp ($date);
if ($overwrite eq 0 and -s "$genome.MITE.raw.fa"){
	print STDERR "$date\tExisting result file $genome.MITE.raw.fa found! Will keep this file without rerunning this module.\n\tPlease specify -overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify MITE candidates from scratch.\n\n";

print STDERR "$date\tWarning: Because MITE-Hunter is too slow and only contribute limited new TIR candidates, it is taken down temporary until a better solution is found.\n\n";
if (0) { #debug line
# run MITE-Hunter
`perl $MITE_Hunter -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c $threads -i $genome > /dev/null 2>&1`;
`cat *_Step8_* > $genome.MITE.raw.fa`;
`rm $genome $genome.index $genome.nhr $genome.nin $genome.nsq`;
} #debug line
	}

# copy result files out
#`cp $genome.MITE.raw.fa ../$genome.MITE.raw.fa`;
if (-s "../$genome.TIR.raw.fa") { #debug line
	`cp ../$genome.TIR.raw.fa ../$genome.MITE.raw.fa`; #debug line
	} else { #debug line
	print STDERR "As a temporary fix, the TIR-Learner is used to mock the MITE-Hunter result. Please run -type tir first.\n\n"; #debug line
	} #debug line
chdir '../..';

# check results
$date=`date`;
chomp ($date);
die "Error: MITE results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.MITE.raw.fa";
if (-s "$genome.EDTA.raw/$genome.MITE.raw.fa"){
	print STDERR "$date\tFinish finding MITE candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The MITE result file has 0 bp!\n\n";
	}

}

#############################
###### HelitronScanner ######
#############################

if ($type eq "helitron" or $type eq "all"){

$date=`date`;
chomp ($date);
print STDERR "$date\tStart to find Helitron candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/Helitron";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
$date=`date`;
chomp ($date);
if ($overwrite eq 0 and -s "$genome.Helitron.raw.fa"){
	print STDERR "$date\tExisting result file $genome.Helitron.raw.fa found! Will keep this file without rerunning this module.\n\tPlease specify -overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify Helitron candidates from scratch.\n\n";

# run HelitronScanner
`sh $HelitronScanner $genome $threads`;

# filter candidates based on repeatness of flanking regions
`perl $format_helitronscanner -genome $genome -sitefilter 1 -minscore 12 -keepshorter 1 -extlen 30 -extout 1`;
`perl $flank_filter -genome $genome -query $genome.HelitronScanner.filtered.ext.fa -miniden 90 -mincov 0.9 -maxct 5 -blastplus $blastplus -t $threads`; #more relaxed
#`perl $flank_filter -genome $genome -query $genome.HelitronScanner.filtered.ext.fa -miniden 80 -mincov 0.8 -maxct 5 -blastplus $blastplus -t $threads`; #more stringent

# remove simple repeats and candidates with simple repeats at terminals
`${mdust}mdust $genome.HelitronScanner.filtered.ext.fa.pass.fa > $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -trf 0 -cleanN 1 -cleanT 1 -f $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted > $genome.Helitron.raw.fa`;
	}

# copy result files out
`cp $genome.Helitron.raw.fa ../$genome.Helitron.raw.fa`;
chdir '../..';

# check results
$date=`date`;
chomp ($date);
die "Error: Helitron results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.Helitron.raw.fa";
if (-s "$genome.EDTA.raw/$genome.Helitron.raw.fa"){
	print STDERR "$date\tFinish finding Helitron candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The Helitron result file has 0 bp!\n\n";
	}

}

$date=`date`;
chomp ($date);
print STDERR "$date\tExecution of EDTA_raw.pl is finished!\n\n";

