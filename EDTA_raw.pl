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
		-threads	[int]	Number of theads to run this script
		-help|-h	Display this help info
\n";

# pre-defined
my $genome = '';
my $species = 'others';
my $type = 'all';
my $threads = 4;
my $script_path = $FindBin::Bin;
my $genometools = "$script_path/bin/genometools-1.5.10/bin/gt";
my $LTR_FINDER = "$script_path/bin/LTR_FINDER_parallel/LTR_FINDER_parallel";
my $LTR_retriever = "$script_path/bin/LTR_retriever/LTR_retriever";
#my $TIR_Learner = "$script_path/bin/TIR-Learner1.13/TIR-Learner.sh";
my $TIR_Learner = "$script_path/bin/TIR-Learner1.15/TIR-Learner_1.15.sh";
my $MITE_Hunter = "$script_path/bin/MITE-Hunter2/MITE_Hunter_manager.pl";
my $HelitronScanner = "$script_path/util/run_helitron_scanner.sh";

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$species = $ARGV[$k+1] if /^-species$/i and $ARGV[$k+1] !~ /^-/;
	$type = lc $ARGV[$k+1] if /^-type$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
  }

print "Check files and dependencies, prepare working directories.\n\n";

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "The GenomeTools is not found in $genometools!\n" unless -s $genometools;
die "The LTR_FINDER_parallel is not found in $LTR_FINDER!\n" unless -s $LTR_FINDER;
die "The LTR_retriever is not found in $LTR_retriever!\n" unless -s $LTR_retriever;
die "The TIR_Learner is not found in $TIR_Learner!\n" unless -s $TIR_Learner;
die "The MITE_Hunter is not found in $MITE_Hunter!\n" unless -s $MITE_Hunter;
die "The HelitronScanner is not found in $HelitronScanner!\n" unless -s $HelitronScanner;

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

print "Start to find LTR candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/LTR";
`ln -s ../../$genome $genome` unless -s $genome;

# run LTRharvest
`$genometools suffixerator -db $genome -indexname $genome -tis -suf -lcp -des -ssp -sds -dna`;
`$genometools ltrharvest -index $genome -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $genome.harvest.scn`;
`rm $genome.esq $genome.lcp $genome.llv $genome.md5 $genome.prj $genome.sds $genome.suf`;

# run LTR_FINDER_parallel
`perl $LTR_FINDER -seq $genome -threads $threads -harvest_out -size 1000000 -time 300`;

# run LTR_retriever
`cat $genome.harvest.scn $genome.finder.combine.scn > $genome.rawLTR.scn`;
`perl $LTR_retriever -genome $genome -inharvest $genome.rawLTR.scn -threads $threads -noanno`;
`cp $genome.LTRlib.fa ../$genome.LTR.raw.fa`;
chdir '../..';

# check results
die "Error: LTR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.LTR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.LTR.raw.fa"){
	print "Finish finding LTR candidates.\n\n";
	} else {
	print "Warning: The LTR result file has 0 bp!\n\n";
	}

}

###########################
######  TIR-Learner  ######
###########################

if ($type eq "tir" or $type eq "all"){

print "Start to find TIR candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/TIR";
`ln -s ../../$genome $genome` unless -s $genome;

# run TIR-Learner
`sh $TIR_Learner -g $genome -s $species -t $threads`;
`cp TIR-Learner-Result/TIR-Learner_FinalAnn.fa ../$genome.TIR.raw.fa`;
chdir '../..';

# check results
die "Error: TIR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.TIR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.TIR.raw.fa"){
	print "Finish finding TIR candidates.\n\n";
	} else {
	`touch "$genome.EDTA.raw/$genome.TIR.raw.fa"`;
	print "Warning: The TIR result file has 0 bp!\n\n";
	}

}


###########################
######  MITE-Hunter  ######
###########################

if ($type eq "mite" or $type eq "all"){

print "Start to find MITE candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/MITE";
`rm -rf genome*`;
`ln -s ../../$genome $genome` unless -s $genome;

# run MITE-Hunter
`perl $MITE_Hunter -l 2 -w 1000 -L 80 -m 1 -S 12345678 -c $threads -i $genome`;
`cat *_Step8_* > $genome.MITE.raw.fa`;
`rm $genome $genome.index $genome.nhr $genome.nin $genome.nsq`;
`cp $genome.MITE.raw.fa ../$genome.MITE.raw.fa`;
chdir '../..';

# check results
die "Error: MITE results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.MITE.raw.fa";
if (-s "$genome.EDTA.raw/$genome.MITE.raw.fa"){
	print "Finish finding MITE candidates.\n\n";
	} else {
	print "Warning: The MITE result file has 0 bp!\n\n";
	}

}

#############################
###### HelitronScanner ######
#############################

if ($type eq "helitron" or $type eq "all"){

print "Start to find Helitron candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/Helitron";
`ln -s ../../$genome $genome` unless -s $genome;

# run HelitronScanner
`sh $HelitronScanner $genome $threads`;
`cp $genome.HelitronScanner.filtered.fa ../$genome.Helitron.raw.fa`;
chdir '../..';

# check results
die "Error: Helitron results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.Helitron.raw.fa";
if (-s "$genome.EDTA.raw/$genome.Helitron.raw.fa"){
	print "Finish finding Helitron candidates.\n\n";
	} else {
	print "Warning: The Helitron result file has 0 bp!\n\n";
	}

}

