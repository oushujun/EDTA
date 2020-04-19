#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use File::Basename;
use Pod::Usage;

########################################################
##### Perform initial searches for TE candidates    ####
##### Shujun Ou (shujun.ou.1@gmail.com, 01/16/2020) ####
########################################################

## Input:
#	$genome

## Output:
#	$genome.LTR.raw.fa, $genome.LTR.intact.fa, $genome.LTR.intact.fa.gff3
#	$genome.TIR.raw.fa, $genome.TIR.intact.fa, $genome.TIR.intact.fa.gff
#	$genome.Helitron.raw.fa, $genome.Helitron.intact.fa, $genome.Helitron.intact.fa.gff

my $usage = "\nObtain raw TE libraries using various structure-based programs
	perl EDTA_raw.pl [options]
		--genome	[File]	The genome FASTA
		--species [rice|maize|others]	Specify the species for identification of TIR candidates. Default: others
		--type	[ltr|tir|helitron|all]	Specify which type of raw TE candidates you want to get. Default: all
		--overwrite	[0|1]	If previous results are found, decide to overwrite (1, rerun) or not (0, default).
		--convert_seq_name	[0|1]	Convert long sequence name to <= 15 characters and remove annotations (1, default) or use the original (0)
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
my $threads = 4;
my $script_path = $FindBin::Bin;
my $LTR_FINDER = "$script_path/bin/LTR_FINDER_parallel/LTR_FINDER_parallel";
my $TIR_Learner = "$script_path/bin/TIR-Learner2.5/TIR-Learner2.5.sh";
my $HelitronScanner = "$script_path/util/run_helitron_scanner.sh";
my $cleanup_misclas = "$script_path/util/cleanup_misclas.pl";
my $get_range = "$script_path/util/get_range.pl";
my $rename_LTR = "$script_path/util/rename_LTR.pl";
my $rename_tirlearner = "$script_path/util/rename_tirlearner.pl";
my $call_seq = "$script_path/util/call_seq_by_list.pl";
my $output_by_list = "$script_path/util/output_by_list.pl";
my $cleanup_tandem = "$script_path/util/cleanup_tandem.pl";
my $get_ext_seq = "$script_path/util/get_ext_seq.pl";
my $format_helitronscanner = "$script_path/util/format_helitronscanner_out.pl";
my $flank_filter = "$script_path/util/flanking_filter.pl";
my $make_gff = "$script_path/util/make_gff_with_intact.pl";
my $genometools = ''; #path to the genometools program
my $TEsorter = ''; #path to the TEsorter program
my $LTR_retriever = ''; #path to the LTR_retriever program
my $mdust = '';
my $blastplus = ''; #path to the blastn program
my $trf = ''; #path to trf
my $beta2 = 0; #0, beta2 is not ready. 1, try it out.
my $help = undef;

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^--genome$/i and $ARGV[$k+1] !~ /^-/;
	$species = $ARGV[$k+1] if /^--species$/i and $ARGV[$k+1] !~ /^-/;
	$type = lc $ARGV[$k+1] if /^--type$/i and $ARGV[$k+1] !~ /^-/;
	$TEsorter = $ARGV[$k+1] if /^--tesorter$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^--blastplus$/i and $ARGV[$k+1] !~ /^-/;
	$mdust = $ARGV[$k+1] if /^--mdust$/i and $ARGV[$k+1] !~ /^-/;
	$overwrite = $ARGV[$k+1] if /^--overwrite$/i and $ARGV[$k+1] !~ /^-/;
	$convert_name = $ARGV[$k+1] if /^--convert_seq_name$/i and $ARGV[$k+1] !~ /^-/;
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
		-message => "At least 1 parameter mandatory:\n1) Input fasta file: --genome\n".
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

if ($type){
	$type =~ s/LTR/ltr/i;
	$type =~ s/TIR/tir/i;
	$type =~ s/Helitron/helitron/i;
	$type =~ s/All/all/i;
	die "The expected value for the type parameter is ltr or tir or helitron or all!\n" unless $type eq "ltr" or $type eq "tir" or $type eq "helitron" or $type eq "all";
	}

die "The expected value for the overwrite parameter is 0 or 1!\n" if $overwrite != 0 and $overwrite != 1;
die "The expected value for the convert_seq_name parameter is 0 or 1!\n" if $convert_name != 0 and $convert_name != 1;
die "The expected value for the threads parameter is an integer!\n" if $threads !~ /^[0-9]+$/;

my $date=`date`;
chomp ($date);
print STDERR "$date\tEDTA_raw: Check dependencies, prepare working directories.\n\n";

# check files and dependencies
die "The LTR_FINDER_parallel is not found in $LTR_FINDER!\n" unless -s $LTR_FINDER;
die "The TIR_Learner is not found in $TIR_Learner!\n" unless -s $TIR_Learner;
die "The script get_range.pl is not found in $get_range!\n" unless -s $get_range;
die "The script rename_LTR.pl is not found in $rename_LTR!\n" unless -s $rename_LTR;
die "The script call_seq_by_list.pl is not found in $call_seq!\n" unless -s $call_seq;
die "The script output_by_list.pl is not found in $output_by_list!\n" unless -s $output_by_list;
die "The script rename_tirlearner.pl is not found in $rename_tirlearner!\n" unless -s $rename_tirlearner;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script get_ext_seq.pl is not found in $get_ext_seq!\n" unless -s $get_ext_seq;
die "The HelitronScanner is not found in $HelitronScanner!\n" unless -s $HelitronScanner;
die "The script format_helitronscanner_out.pl is not found in $format_helitronscanner!\n" unless -s $format_helitronscanner;
die "The script flanking_filter.pl is not found in $flank_filter!\n" unless -s $flank_filter;
die "The script make_gff_with_intact.pl is not found in $make_gff!\n" unless -s $make_gff;
# GenomeTools
$genometools=`which gt 2>/dev/null` if $genometools eq '';
$genometools=~s/gt\n//;
$genometools="$genometools/" if $genometools ne '' and $genometools !~ /\/$/;
die "gt is not exist in the genometools path $genometools!\n" unless -X "${genometools}gt";
# LTR_retriever
$LTR_retriever=`which LTR_retriever 2>/dev/null` if $LTR_retriever eq '';
$LTR_retriever=~s/LTR_retriever\n//;
$LTR_retriever="$LTR_retriever/" if $LTR_retriever ne '' and $LTR_retriever !~ /\/$/;
die "LTR_retriever is not exist in the LTR_retriever path $LTR_retriever!\n" unless -X "${LTR_retriever}LTR_retriever";
# TEsorter
$TEsorter=`which TEsorter 2>/dev/null` if $TEsorter eq '';
$TEsorter=~s/TEsorter\n//;
$TEsorter="$TEsorter/" if $TEsorter ne '' and $TEsorter !~ /\/$/;
die "TEsorter is not exist in the TEsorter path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
# makeblastdb, blastn
$blastplus=`which blastn 2>/dev/null` if $blastplus eq '';
$blastplus=~s/blastn\n//;
$blastplus="$blastplus/" if $blastplus ne '' and $blastplus !~ /\/$/;
die "makeblastdb is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}makeblastdb";
die "blastn is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";
# mdust
$mdust=`which mdust 2>/dev/null` if $mdust eq '';
$mdust=~s/mdust\n//;
die "mdust is not exist in the mdust path $mdust!\n" unless -X "${mdust}mdust";
# trf
$trf=`which trf 2>/dev/null` if $trf eq '';
$trf=~s/\n$//;
`$trf 2>/dev/null`;
die "Error: No Tandem Repeat Finder is working on the current system.
        Please report it to https://github.com/oushujun/EDTA/issues" if $?==32256;
die "\n\tTandem Repeat Finder not found!\n\n" unless $trf ne '';

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
	$date = `date`;
	chomp ($date);
	die "$date\tERROR: Identical sequence IDs found in the provided genome! Please resolve this issue and try again.\n";
	}

if ($convert_name == 1){
if (-s "$genome.mod"){
	$genome = "$genome.mod";
	} else {
# remove sequence annotations (content after the first space in sequence names)
`perl -nle 'my \$info=(split)[0]; print \$info' $genome > $genome.mod`;

# try to shortern sequences
if ($id_len > 15){
	$date = `date`;
	chomp ($date);
	print "$date\tThe longest sequence ID in the genome contains $id_len characters, which is longer than the limit (15)\n";
	print "\t\t\t\tTrying to reformat seq IDs...\n\t\t\t\tAttempt 1...\n";
	`perl -lne 'chomp; if (s/^>+//) {s/^\\s+//; \$_=(split)[0]; s/(.{1,15}).*/>\$1/g;} print "\$_"' $genome.mod > $genome.temp`;
	my $new_id = `grep \\> $genome.temp|sort -u|wc -l`;
	$date = `date`;
	chomp ($date);
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
			die "$date\tERROR: Fail to convert seq IDs to less than 15 characters! Please provide a genome with shorter seq IDs.\n\n";
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
	print STDERR "$date\tExisting result file $genome.LTR.raw.fa found! Will keep this file without rerunning this module.\n\tPlease specify --overwrite 1 if you want to rerun this module.\n\n";
	} else {
	print STDERR "$date\tIdentify LTR retrotransposon candidates from scratch.\n\n";

# run LTRharvest
if ($overwrite eq 0 and -s "$genome.harvest.scn"){
	print STDERR "$date\tExisting raw result $genome.harvest.scn found! Will use this for further analyses.\n\n";
	} else {
	`${genometools}gt suffixerator -db $genome -indexname $genome -tis -suf -lcp -des -ssp -sds -dna -mirrored 2>/dev/null`;
	`${genometools}gt ltrharvest -index $genome -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > $genome.harvest.scn 2>/dev/null`;
	`rm $genome.des $genome.esq $genome.lcp $genome.llv $genome.md5 $genome.prj $genome.sds $genome.ssp $genome.suf 2>/dev/null`;
}

# run LTR_FINDER_parallel
if ($overwrite eq 0 and -s "$genome.finder.combine.scn"){
	print STDERR "$date\tExisting raw result $genome.finder.combine.scn found! Will use this for further analyses.\n\n";
	} else {
	`perl $LTR_FINDER -seq $genome -threads $threads -harvest_out -size 1000000 -time 300`;
	}

# run LTR_retriever
`cat $genome.harvest.scn $genome.finder.combine.scn > $genome.rawLTR.scn`;
`${LTR_retriever}LTR_retriever -genome $genome -inharvest $genome.rawLTR.scn -threads $threads -noanno`;
#`for i in *.mod.*; do mv \$i \$(echo \$i|sed 's/\\.mod\\././'); done > /dev/null 2>&1`;
#`mv $genome.mod $genome` if -s "$genome.mod";

# get intact LTR elements
if (0){ #old, overly inclusive module
	`perl $get_range 1 $genome.rawLTR.scn $genome -f -g -max_ratio 50`;
	`cat $genome.rawLTR.scn.full | sort -fu > $genome.rawLTR.scn.full.uniq`;
	`perl $call_seq $genome.rawLTR.scn.full.uniq -C $genome > $genome.LTR`;
	`perl -i -nle 's/\\|.*//; print \$_' $genome.LTR`;
	`perl $get_ext_seq $genome $genome.LTR`;
	`perl $flank_filter -genome $genome -query $genome.LTR.ext30.fa -miniden 90 -mincov 0.9 -maxct 20 -blastplus $blastplus -t $threads`;

        # recover superfamily info
	`perl $output_by_list 1 $genome.LTR 1 $genome.LTR.ext30.fa.pass.fa -FA -MSU0 -MSU1 > $genome.LTR.ext30.fa.pass.fa.ori`;
	}

# get full-length LTR from pass.list
`awk '{if (\$1 !~ /#/) print \$1"\\t"\$1}' $genome.pass.list | perl $call_seq - -C $genome > $genome.LTR.intact.fa.ori`;
`perl -i -nle 's/\\|.*//; print \$_' $genome.LTR.intact.fa.ori`;

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
`perl $rename_LTR $genome $genome.LTR.intact.fa $genome.defalse > $genome.LTR.intact.fa.anno`;
`mv $genome.LTR.intact.fa.anno $genome.LTR.intact.fa`;
`cp $genome.LTR.intact.fa $genome.LTR.intact.fa.gff3 ../`;

`rm $genome`;
	}

# remove non-LTR sequence
#`${TEsorter}TEsorter $genome.LTRlib.fa -p $threads`;
#`awk '{if (\$2!="pararetrovirus" && \$2!="LTR")print \$0}' $genome.LTRlib.fa.rexdb.cls.tsv > $genome.LTRlib.fa.rexdb.cls.tsv.nonLTR`;
#`perl $output_by_list 1 $genome.LTRlib.fa 1 $genome.LTRlib.fa.rexdb.cls.tsv.nonLTR -ex -FA > $genome.LTRlib.fa.cln`;

# copy result files out
#`cp $genome.LTRlib.fa.cln $genome.LTR.raw.fa`;
#`cp $genome.LTRlib.fa.cln ../$genome.LTR.raw.fa`;
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
	print STDERR "Species: $species\n";

	# run TIR-Learner
	`sh $TIR_Learner -g $genome -s $species -t $threads -l $maxint`;
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
	} else {
	`cp $genome.TIR.ext30.fa.pass.fa.dusted.cln $genome.TIR.raw.fa`;
	}

	# get gff of intact TIR elements
	`perl -nle 's/\\-\\+\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class='DNA'; \$class='MITE' if \$e-\$s+1 <= 600; \$anno = "ID=\$chr:\$s..\$e#\$class/\$supfam;\$anno"; print "\$chr \$method \$class/\$supfam \$s \$e . . . \$anno \$chr:\$s..\$e"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3| perl $output_by_list 10 - 1 $genome.TIR.raw.fa -MSU0 -MSU1 | awk '{\$10=""; print \$0}' | perl -nle 's/\\s+/\\t/g; print \$_' > $genome.TIR.intact.fa.gff`;

	`cp $genome.TIR.raw.fa $genome.TIR.intact.fa`;
	`cp $genome.TIR.intact.fa $genome.TIR.intact.fa.gff ../`;
	if ($beta2 == 1){
		`cp $genome.TIR.ext30.fa.pass.fa.dusted.cln.cln.list $genome.TIR.intact.fa.anno.list`;
		`cp $genome.LTR.intact.fa.anno.list ../`;
		}

	}

# copy result files out
`cp $genome.TIR.raw.fa ../$genome.TIR.raw.fa`;
chdir '../..';

# check results
$date=`date`;
chomp ($date);
die "Error: TIR results not found!\n\n" unless -e "$genome.EDTA.raw/$genome.TIR.raw.fa";
if (-s "$genome.EDTA.raw/$genome.TIR.raw.fa"){
	print STDERR "$date\tFinish finding TIR candidates.\n\n";
	} else {
	`touch "$genome.EDTA.raw/$genome.TIR.raw.fa"`;
	print STDERR "Warning: The TIR result file has 0 bp!\n\n";
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
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted | perl -nle 's/^(>.*)\$/\$1#DNA\\/Helitron/; print \$_' > $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted.cln`;

if ($beta2 == 1){
# annotate and remove non-Helitron candidates
`${TEsorter}TEsorter $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted.cln -p $threads`;
`perl $cleanup_misclas $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted.cln.rexdb.cls.tsv`;
`mv $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted.cln.cln $genome.Helitron.raw.fa`;
} else {
`cp $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted.cln $genome.Helitron.raw.fa`;
}

# get intact Helitrons and gff
`cp $genome.Helitron.raw.fa $genome.Helitron.intact.fa`;
`perl $make_gff $genome.Helitron.intact.fa`;
`cp $genome.Helitron.intact.fa $genome.Helitron.intact.fa.gff ../`;
if ($beta2 == 1){
	`cp $genome.HelitronScanner.filtered.ext.fa.pass.fa.dusted.cln.cln.list $genome.Helitron.intact.fa.anno.list`;
	`cp $genome.Helitron.intact.fa.anno.list ../`;
	}
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
