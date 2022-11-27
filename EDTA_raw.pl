#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use Pod::Usage;

########################################################
##### Perform initial searches for TE candidates    ####
##### Shujun Ou (shujun.ou.1@gmail.com, 07/16/2020) ####
########################################################

## Input:
#	$genome

## Output:
#	$genome.LTR.raw.fa, $genome.LTR.intact.fa, $genome.LTR.intact.gff3
#	$genome.TIR.raw.fa, $genome.TIR.intact.fa, $genome.TIR.intact.gff3
#	$genome.Helitron.raw.fa, $genome.Helitron.intact.fa, $genome.Helitron.intact.gff3

my $usage = "\nObtain raw TE libraries using various structure-based programs

perl EDTA_raw.pl [options]
	--genome	[File]	The genome FASTA
	--species [rice|maize|others]	Specify the species for identification
					of TIR candidates. Default: others
	--type	[ltr|tir|helitron|all]	Specify which type of raw TE candidates
					you want to get. Default: all
	--overwrite	[0|1]	If previous results are found, decide to
				overwrite (1, rerun) or not (0, default).
	--convert_seq_name	[0|1]	Convert long sequence name to <= 15
					characters and remove annotations (1,
					default) or use the original (0)
	--u [float]	Neutral mutation rate to calculate the age of intact LTR elements.
			Intact LTR age is found in this file: *EDTA_raw/LTR/*.pass.list.
			Default: 1.3e-8 (per bp per year, from rice).
	--tesorter	[path]	Path to the TEsorter program. (default: find from ENV)
	--repeatmasker	[path]	Path to the RepeatMasker program. (default: find from ENV)
	--threads|-t	[int]	Number of theads to run this script. Default: 4
	--help|-h	Display this help info
\n";

# pre-defined
my $genome = '';
my $species = 'others';
my $type = 'all';
my $overwrite = 0; #0, no rerun. 1, rerun even old results exist.
my $convert_name = 1; #0, use original seq names; 1 shorten names.
my $maxint = 5000; #maximum interval length (bp) between TIRs (for GRF in TIR-Learner)
my $miu = 1.3e-8; #mutation rate, per bp per year, from rice
my $threads = 4;
my $script_path = $FindBin::Bin;
my $LTR_FINDER = "$script_path/bin/LTR_FINDER_parallel/LTR_FINDER_parallel";
my $LTR_HARVEST = "$script_path/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel";
my $TIR_Learner = "$script_path/bin/TIR-Learner2.5/TIR-Learner2.5.sh";
my $HelitronScanner = "$script_path/util/run_helitron_scanner.sh";
my $cleanup_misclas = "$script_path/util/cleanup_misclas.pl";
my $get_range = "$script_path/util/get_range.pl";
my $rename_LTR = "$script_path/util/rename_LTR_skim.pl";
my $filter_gff = "$script_path/util/filter_gff3.pl";
my $rename_tirlearner = "$script_path/util/rename_tirlearner.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $get_ext_seq = "$script_path/util/get_ext_seq.pl";
my $format_helitronscanner = "$script_path/util/format_helitronscanner_out.pl";
my $flank_filter = "$script_path/util/flanking_filter.pl";
my $make_bed = "$script_path/util/make_bed_with_intact.pl";
my $bed2gff = "$script_path/util/bed2gff.pl";
my $genometools = ''; #path to the genometools program
my $repeatmasker = ''; #path to the RepeatMasker program
my $LTR_retriever = ''; #path to the LTR_retriever program
my $TEsorter = ''; #path to the TEsorter program
my $blastplus = ''; #path to the blastn program
my $mdust = ''; #path to mdust
my $trf = ''; #path to trf
my $GRF = ''; #path to GRF
my $beta2 = 0; #0, beta2 is not ready. 1, try it out.
my $help = undef;

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^--genome$/i and $ARGV[$k+1] !~ /^-/;
	$species = $ARGV[$k+1] if /^--species$/i and $ARGV[$k+1] !~ /^-/;
	$type = lc $ARGV[$k+1] if /^--type$/i and $ARGV[$k+1] !~ /^-/;
	$overwrite = $ARGV[$k+1] if /^--overwrite$/i and $ARGV[$k+1] !~ /^-/;
	$convert_name = $ARGV[$k+1] if /^--convert_seq_name$/i and $ARGV[$k+1] !~ /^-/;
	$miu = $ARGV[$k+1] if /^--u$/i and $ARGV[$k+1] !~ /^-/;
	$genometools = $ARGV[$k+1] if /^--genometools/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^--repeatmasker$/i and $ARGV[$k+1] !~ /^-/;
	$LTR_retriever = $ARGV[$k+1] if /^--ltrretriever/i and $ARGV[$k+1] !~ /^-/;
	$TEsorter = $ARGV[$k+1] if /^--tesorter$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^--blastplus$/i and $ARGV[$k+1] !~ /^-/;
	$mdust = $ARGV[$k+1] if /^--mdust$/i and $ARGV[$k+1] !~ /^-/;
	$trf = $ARGV[$k+1] if /^--trf_path$/i and $ARGV[$k+1] !~ /^-/;
	$GRF = $ARGV[$k+1] if /^--GRF$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^--threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
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

die "The expected value for the type parameter is ltr or tir or helitron or all!\n" unless $type eq "ltr" or $type eq "tir" or $type eq "helitron" or $type eq "all";

# check bolean
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"};
if ($convert_name != 0 and $convert_name != 1){ die "The expected value for the convert_seq_name parameter is 0 or 1!\n"};
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"};
if ($miu !~ /[0-9\.e\-]+/){ die "The expected value for the u parameter is float value without units!\n"}

chomp (my $date = `date`);
print STDERR "$date\tEDTA_raw: Check dependencies, prepare working directories.\n\n";

# check files and dependencies
die "The LTR_FINDER_parallel is not found in $LTR_FINDER!\n" unless -s $LTR_FINDER;
die "The LTR_HARVEST_parallel is not found in $LTR_HARVEST!\n" unless -s $LTR_HARVEST;
die "The TIR_Learner is not found in $TIR_Learner!\n" unless -s $TIR_Learner;
die "The script get_range.pl is not found in $get_range!\n" unless -s $get_range;
die "The script rename_LTR.pl is not found in $rename_LTR!\n" unless -s $rename_LTR;
die "The script filter_gff3.pl is not found in $filter_gff!\n" unless -s $filter_gff;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script get_ext_seq.pl is not found in $get_ext_seq!\n" unless -s $get_ext_seq;
die "The HelitronScanner is not found in $HelitronScanner!\n" unless -s $HelitronScanner;
die "The script format_helitronscanner_out.pl is not found in $format_helitronscanner!\n" unless -s $format_helitronscanner;
die "The script flanking_filter.pl is not found in $flank_filter!\n" unless -s $flank_filter;
die "The script bed2gff.pl is not found in $bed2gff!\n" unless -s $bed2gff;
die "The script make_bed_with_intact.pl is not found in $make_bed!\n" unless -s $make_bed;

# GenomeTools
chomp ($genometools=`which gt 2>/dev/null`) if $genometools eq '';
$genometools =~ s/\s+$//;
$genometools = dirname($genometools) unless -d $genometools;
$genometools="$genometools/" if $genometools ne '' and $genometools !~ /\/$/;
die "Error: gt is not found in the genometools path $genometools!\n" unless -X "${genometools}gt";
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
`rm dummy060817.fa.$rand*`;
# LTR_retriever
chomp ($LTR_retriever=`which LTR_retriever 2>/dev/null`) if $LTR_retriever eq '';
$LTR_retriever =~ s/\s+$//;
$LTR_retriever = dirname($LTR_retriever) unless -d $LTR_retriever;
$LTR_retriever="$LTR_retriever/" if $LTR_retriever ne '' and $LTR_retriever !~ /\/$/;
die "Error: LTR_retriever is not found in the LTR_retriever path $LTR_retriever!\n" unless -X "${LTR_retriever}LTR_retriever";
# TEsorter
chomp ($TEsorter=`which TEsorter 2>/dev/null`) if $TEsorter eq '';
$TEsorter =~ s/\s+$//;
$TEsorter = dirname($TEsorter) unless -d $TEsorter;
$TEsorter="$TEsorter/" if $TEsorter ne '' and $TEsorter !~ /\/$/;
die "Error: TEsorter is not found in the TEsorter path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
# makeblastdb, blastn
chomp ($blastplus=`which makeblastdb 2>/dev/null`) if $blastplus eq '';
$blastplus =~ s/\s+$//;
$blastplus = dirname($blastplus) unless -d $blastplus;
$blastplus="$blastplus/" if $blastplus ne '' and $blastplus !~ /\/$/;
die "Error: makeblastdb is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}makeblastdb";
die "Error: blastn is not found in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";
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

# make a softlink to the genome
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

if ($convert_name == 1){
if (-s "$genome.mod"){
	$genome = "$genome.mod";
	} else {

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
			`rm $genome.temp`;
			die "$date\tERROR: Fail to convert seq IDs to <= $id_len_max characters! Please provide a genome with shorter seq IDs.\n\n";
			}
		}
	}
$genome = "$genome.mod";
}
}

# Make working directories
`mkdir $genome.EDTA.raw` unless -e "$genome.EDTA.raw" && -d "$genome.EDTA.raw";
`mkdir $genome.EDTA.raw/LTR` unless -e "$genome.EDTA.raw/LTR" && -d "$genome.EDTA.raw/LTR";
`mkdir $genome.EDTA.raw/TIR` unless -e "$genome.EDTA.raw/TIR" && -d "$genome.EDTA.raw/TIR";
`mkdir $genome.EDTA.raw/Helitron` unless -e "$genome.EDTA.raw/Helitron" && -d "$genome.EDTA.raw/Helitron";


###########################
###### LTR_retriever ######
###########################

if ($type eq "ltr" or $type eq "all"){

chomp ($date = `date`);
print STDERR "$date\tStart to find LTR candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/LTR";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
chomp ($date = `date`);
if ($overwrite eq 0 and -s "$genome.LTR.raw.fa"){
	print STDERR "$date\tExisting result file $genome.LTR.raw.fa found!\n\t\t\t\tWill keep this file without rerunning this module.\n\t\t\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify LTR retrotransposon candidates from scratch.\n\n";

# run LTRharvest
if ($overwrite eq 0 and -s "$genome.harvest.combine.scn"){
	print STDERR "$date\tExisting raw result $genome.harvest.scn found!\n\t\t\t\tWill use this for further analyses.\n\n";
	} else {
	`perl $LTR_HARVEST -seq $genome -threads $threads -gt $genometools -size 1000000 -time 300`;
	}

# run LTR_FINDER_parallel
if ($overwrite eq 0 and -s "$genome.finder.combine.scn"){
	print STDERR "$date\tExisting raw result $genome.finder.combine.scn found!\n\t\t\t\tWill use this for further analyses.\n\n";
	} else {
	`perl $LTR_FINDER -seq $genome -threads $threads -harvest_out -size 1000000 -time 300`;
	}

# run LTR_retriever
`cat $genome.harvest.combine.scn $genome.finder.combine.scn > $genome.rawLTR.scn`;
`${LTR_retriever}LTR_retriever -genome $genome -inharvest $genome.rawLTR.scn -u $miu -threads $threads -noanno -trf_path $trf -blastplus $blastplus -repeatmasker $repeatmasker`;

# get full-length LTR from pass.list
`awk '{if (\$1 !~ /#/) print \$1"\\t"\$1}' $genome.pass.list | perl $call_seq - -C $genome > $genome.LTR.intact.fa.ori`;
`perl -i -nle 's/\\|.*//; print \$_' $genome.LTR.intact.fa.ori`;
`perl $rename_LTR $genome.LTR.intact.fa.ori $genome.defalse > $genome.LTR.intact.fa.anno`;
`mv $genome.LTR.intact.fa.anno $genome.LTR.intact.fa.ori`;

# remove simple repeats and candidates with simple repeats at terminals
`${mdust}mdust $genome.LTR.intact.fa.ori > $genome.LTR.intact.fa.ori.dusted`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.LTR.intact.fa.ori.dusted > $genome.LTR.intact.fa.ori.dusted.cln`;

if ($beta2 == 1){
	# annotate and remove non-LTR candidates
	`${TEsorter}TEsorter $genome.LTR.intact.fa.ori.dusted.cln -p $threads`;
	`perl $cleanup_misclas $genome.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv`;
	`mv $genome.LTR.intact.fa.ori.dusted.cln.cln $genome.LTR.intact.fa`;
	`mv $genome.LTR.intact.fa.ori.dusted.cln.cln.list $genome.LTR.intact.fa.anno.list`;
	`cp $genome.LTR.intact.fa.anno.list ../`;
	} else {
	`mv $genome.LTR.intact.fa.ori.dusted.cln $genome.LTR.intact.fa`;
	}

# generate annotated output and gff
`perl $output_by_list 1 $genome.LTR.intact.fa.ori 1 $genome.LTR.intact.fa -FA -ex|grep \\>|perl -nle 's/>//; print "Name\\t\$_"' > $genome.LTR.intact.fa.ori.rmlist`;
`perl $filter_gff $genome.pass.list.gff3 $genome.LTR.intact.fa.ori.rmlist | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' > $genome.LTR.intact.gff3`;
`rm $genome`;
	}

# remove non-LTR sequence
#`${TEsorter}TEsorter $genome.LTRlib.fa -p $threads`;
#`awk '{if (\$2!="pararetrovirus" && \$2!="LTR")print \$0}' $genome.LTRlib.fa.rexdb.cls.tsv > $genome.LTRlib.fa.rexdb.cls.tsv.nonLTR`;
#`perl $output_by_list 1 $genome.LTRlib.fa 1 $genome.LTRlib.fa.rexdb.cls.tsv.nonLTR -ex -FA > $genome.LTRlib.fa.cln`;

# copy result files out
#`cp $genome.LTRlib.fa.cln $genome.LTR.raw.fa`;
#`cp $genome.LTRlib.fa.cln ../$genome.LTR.raw.fa`;
`touch $genome.LTRlib.fa` unless -e "$genome.LTRlib.fa";
`cp $genome.LTRlib.fa $genome.LTR.raw.fa`;
`cp $genome.LTRlib.fa ../$genome.LTR.raw.fa`;
`cp $genome.LTR.intact.fa $genome.LTR.intact.gff3 ../`;
chdir '../..';

# check results
chomp ($date = `date`);
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

chomp ($date = `date`);
print STDERR "$date\tStart to find TIR candidates.\n\n";

# enter the working directory, filter out short sequences and create genome softlink
chdir "$genome.EDTA.raw/TIR";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
chomp ($date = `date`);
if ($overwrite eq 0 and -s "$genome.TIR.raw.fa"){
	print STDERR "$date\tExisting result file $genome.TIR.raw.fa found!\n\t\t\t\tWill keep this file without rerunning this module.\n\t\t\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify TIR candidates from scratch.\n\n";
	print STDERR "Species: $species\n";

	# run TIR-Learner
	`bash $TIR_Learner -g $genome -s $species -t $threads -l $maxint`;
	`perl $rename_tirlearner ./TIR-Learner-Result/TIR-Learner_FinalAnn.fa | perl -nle 's/TIR-Learner_//g; print \$_' > $genome.TIR`;

	# clean raw predictions with flanking alignment
	`perl $get_ext_seq $genome $genome.TIR`;
	`perl $flank_filter -genome $genome -query $genome.TIR.ext30.fa -miniden 90 -mincov 0.9 -maxct 20 -blastplus $blastplus -t $threads`;

	# recover superfamily info
	`perl $output_by_list 1 $genome.TIR 1 $genome.TIR.ext30.fa.pass.fa -FA -MSU0 -MSU1 > $genome.TIR.ext30.fa.pass.fa.ori`;

	# remove simple repeats and candidates with simple repeats at terminals
	`${mdust}mdust $genome.TIR.ext30.fa.pass.fa.ori > $genome.TIR.ext30.fa.pass.fa.dusted`;
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.TIR.ext30.fa.pass.fa.dusted > $genome.TIR.ext30.fa.pass.fa.dusted.cln`;

	if ($beta2 == 1){
	# annotate and remove non-TIR candidates
	`${TEsorter}TEsorter $genome.TIR.ext30.fa.pass.fa.dusted.cln -p $threads`;
	`perl $cleanup_misclas $genome.TIR.ext30.fa.pass.fa.dusted.cln.rexdb.cls.tsv`;
	`mv $genome.TIR.ext30.fa.pass.fa.dusted.cln.cln $genome.TIR.raw.fa`;
	`cp $genome.TIR.ext30.fa.pass.fa.dusted.cln.cln.list $genome.TIR.intact.fa.anno.list`;
	`cp $genome.LTR.intact.fa.anno.list ../`;
	} else {
	`cp $genome.TIR.ext30.fa.pass.fa.dusted.cln $genome.TIR.raw.fa`;
	}

	# get gff3 of intact TIR elements
	`perl -nle 's/\\-\\+\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class='DNA'; \$class='MITE' if \$e-\$s+1 <= 600; my (\$tir, \$iden, \$tsd)=(\$1, \$2/100, \$3) if \$anno=~/TIR:(.*)_([0-9.]+)_TSD:([a-z0-9._]+)_LEN/i; print "\$chr \$s \$e \$chr:\$s..\$e \$class/\$supfam structural \$iden . . . TSD=\$tsd;TIR=\$tir"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3 | perl $output_by_list 4 - 1 $genome.TIR.raw.fa -MSU0 -MSU1 > $genome.TIR.intact.bed`;
	`perl $bed2gff $genome.TIR.intact.bed TIR > $genome.TIR.intact.gff3`;
	`cp $genome.TIR.raw.fa $genome.TIR.intact.fa`;
	}

# copy result files out
`touch $genome.TIR.raw.fa` unless -e "$genome.TIR.raw.fa";
`cp $genome.TIR.raw.fa $genome.TIR.intact.fa $genome.TIR.intact.gff3 $genome.TIR.intact.bed ../`;
chdir '../..';

# check results
chomp ($date = `date`);
die "Error: TIR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.TIR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.TIR.raw.fa"){
	print STDERR "$date\tFinish finding TIR candidates.\n\n";
	} else {
	print STDERR "Warning: The TIR result file has 0 bp!\n\n";
	}

}


#############################
###### HelitronScanner ######
#############################

if ($type eq "helitron" or $type eq "all"){

chomp ($date = `date`);
print STDERR "$date\tStart to find Helitron candidates.\n\n";

# enter the working directory and create genome softlink
chdir "$genome.EDTA.raw/Helitron";
`ln -s ../../$genome $genome` unless -s $genome;

# Try to recover existing results
chomp ($date = `date`);
if ($overwrite eq 0 and -s "$genome.Helitron.raw.fa"){
	print STDERR "$date\tExisting result file $genome.Helitron.raw.fa found!\n\t\t\t\tWill keep this file without rerunning this module.\n\t\t\t\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify Helitron candidates from scratch.\n\n";

# run HelitronScanner
`sh $HelitronScanner $genome $threads`;

# filter candidates based on repeatness of flanking regions
`perl $format_helitronscanner -genome $genome -sitefilter 1 -minscore 12 -keepshorter 1 -extlen 30 -extout 1`;
`perl $flank_filter -genome $genome -query $genome.HelitronScanner.filtered.ext.fa -miniden 90 -mincov 0.9 -maxct 5 -blastplus $blastplus -t $threads`; #more relaxed
#`perl $flank_filter -genome $genome -query $genome.HelitronScanner.filtered.ext.fa -miniden 80 -mincov 0.8 -maxct 5 -blastplus $blastplus -t $threads`; #more stringent

# remove simple repeats and candidates with simple repeats at terminals
`perl $output_by_list 1 $genome.HelitronScanner.filtered.fa 1 $genome.HelitronScanner.filtered.ext.fa.pass.fa -FA > $genome.HelitronScanner.filtered.fa.pass.fa`;
`${mdust}mdust $genome.HelitronScanner.filtered.fa.pass.fa > $genome.HelitronScanner.filtered.fa.pass.fa.dusted`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.HelitronScanner.filtered.fa.pass.fa.dusted | perl -nle 's/^(>.*)\\s+(.*)\$/\$1#DNA\\/Helitron\\t\$2/; print \$_' > $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln`;

if ($beta2 == 1){
# annotate and remove non-Helitron candidates
`${TEsorter}TEsorter $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln -p $threads`;
`perl $cleanup_misclas $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.rexdb.cls.tsv`;
`mv $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln $genome.Helitron.raw.fa`;
`cp $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln.list $genome.Helitron.intact.fa.anno.list`;
`cp $genome.Helitron.intact.fa.anno.list ../`;
} else {
`cp $genome.HelitronScanner.filtered.fa.pass.fa.dusted.cln $genome.Helitron.raw.fa`;
}

# get intact Helitrons and gff3
`cp $genome.Helitron.raw.fa $genome.Helitron.intact.fa`;
`perl $make_bed $genome.Helitron.intact.fa > $genome.Helitron.intact.bed`;
`perl $bed2gff $genome.Helitron.intact.bed HEL > $genome.Helitron.intact.gff3`;
	}

# copy result files out
`touch $genome.Helitron.raw.fa` unless -e "$genome.Helitron.raw.fa";
`cp $genome.Helitron.raw.fa $genome.Helitron.intact.fa $genome.Helitron.intact.gff3 $genome.Helitron.intact.bed ../`;
chdir '../..';

# check results
chomp ($date = `date`);
die "Error: Helitron results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.Helitron.raw.fa";
if (-s "$genome.EDTA.raw/$genome.Helitron.raw.fa"){
	print STDERR "$date\tFinish finding Helitron candidates.\n\n";
	} else {
	print STDERR "$date\tWarning: The Helitron result file has 0 bp!\n\n";
	}

}

chomp ($date = `date`);
print STDERR "$date\tExecution of EDTA_raw.pl is finished!\n\n";
