#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

#To test the annotation performance of TE libraries by comparing to a reference annotation
#
#sens = P(test+|TE+) = TP/(TP+FN)
#spec = P(test-|TE-) = TN/(FP+TN)
#accu = (TP+TN)/(TP+TN+FP+FN)
#prec = TP/(TP+FP)
#FDR = 1-prec
#F1 = 2TP/(2TP+FP+FN)
#
#------------------------------  Whole genome
#+++++++++++++++++++	   	 Standard lib
#    ++++++++++++++++++++	 Test lib
# FN [    TP       ] FP  [ TN ]
#
#Author: Shujun Ou (oushujun@msu.edu), 03/08/2015
#Updates:
#	05/03/2023
#	05/09/2019
#
my $usage="\n\tTo test the annotation performance of TE libraries by comparing to a reference annotation

	perl lib-test.pl -genome genome.fasta -std genome.stdlib.RM.out -tst genome.testlib.RM.out -cat [options]
		-genome	[file]	FASTA format genome sequence
		-std	[file]	RepeatMasker .out file of the standard library
		-tst	[file]	RepeatMasker .out file of the test library
		-cat	[string]	Testing TE category. Use one of LTR|nonLTR|LINE|SINE|TIR|MITE|Helitron|Total|Classified
		-N	[0|1]	Include Ns in total length of the genome. Defaule: 0 (not include Ns).
		-unknown	[0|1]	Include unknown annotations to the testing category. This should be used when
					the test library has no classification and you assume they all belong to the
					target category specified by -cat. Default: 0 (not include unknowns)
		-rand	[int]	A randum number used to identify the current run. (default: generate automatically)
		-threads|-t	[int]	Number of threads to run this program. Default: 4
		-debug	[0|1]	Retain intermediate files for debugging. Default: 0 (delete intermediate files)
\n";

my $script_path = dirname(__FILE__);

my $genome='';
my $std_out='';
my $tst_out='';
my $category='total'; #TE category for testing, default total (all TEs)
my $unknown="0"; #0, not include unknown annotation in -tst, especially when -tst contains only 1 category (strict). 1 will include unknown..
my $rand=int(rand(1000000));
my $threads = 4;
my $debug = 0;

my $includeN=0; #0 will not include N in total length of the genome, 1 will include N

my $k=0;
foreach (@ARGV){
	$genome=$ARGV[$k+1] if /^-genome$/i;
	$std_out=$ARGV[$k+1] if /^-std$/i;
	$tst_out=$ARGV[$k+1] if /^-tst$/i;
	$includeN=$ARGV[$k+1] if /^-N$/i;
	$category=lc $ARGV[$k+1] if /^-cat$/i;
	$unknown=$ARGV[$k+1] if /^-unknown$/i;
	$rand=$ARGV[$k+1] if /^-rand$/i;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i;
	$debug=$ARGV[$k+1] if /^-debug$/i;
	$k++;
	}

## check files
die $usage unless -s $genome and -s $std_out and -s $tst_out;

## create softlinks to input files with unique ID in the current folder
my $genome_base = basename($genome);
my $stdout_base = basename($std_out);
my $tstout_base = basename($tst_out);
`ln -s $genome $genome_base.rand$rand`;
`ln -s $std_out $stdout_base.rand$rand`;
`ln -s $tst_out $tstout_base.rand$rand`;
$genome="$genome_base.rand$rand";
$std_out="$stdout_base.rand$rand";
$tst_out="$tstout_base.rand$rand";

## get whole genome range
open FASTA, "<$genome" or die $usage;
open All, ">$genome.list" or die $usage;

$/="\n>";
while (<FASTA>){
	s/>//g;
	next if /^\s+$/;
	my ($chr,$seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	my $len=length $seq;
	my $N=$seq=~tr/Nn//;
	$len=$len-$N unless $includeN;
	print All "$chr\t1\t$len\n";
	}
$/="\n";
close FASTA;
close All;

## Define TE categories. Categories are case insensitive, but hash keys of %category are all lowercases. Each category are specified in the format of "'cate1\\|cate2'". Don't forget ''.
my %category;
$category{'ltr'}="'RLG\\|RLC\\|RLB\\|RLR\\|RLE\\|\\\\s+LTR\\|RLX\\|Gypsy\\|Copia'";
$category{'nonltr'}="'SINE\\|LINE\\|Penelope\\|RIT\\|RIL\\|RST\\|RIX\\|RSX\\|nonLTR\\|\\\\s+YR'";
$category{'line'}="'LINE\\|RIL\\|RIT\\|RIX\\|Penelope'";
$category{'sine'}="'SINE\\|RST\\|RSX'";
$category{'tir'}="'TIR\\|MITE\\|hAT\\|hAT-Ac\\|MULE\\|MLE\\|MuDR\\|Tourist\\|CACT\\|PILE\\|POLE\\|Stowaway\\|TcMar-Stowaway\\|PIF\\|Harbinger\\|Tc1\\|En-Spm\\|EnSpm\\|CMC-EnSpm\\|PiggyBac\\|Mirage\\|P-element\\|Transib\\|DTA\\|DTH\\|DTT\\|DTM\\|DTC\\|DTA\\|TIR\\|DTX\\|DTR\\|DTE\\|Merlin\\|DTP\\|DTB\\|polinton'";
$category{'mite'}="MITE";
$category{'helitron'}="'Helitron\\|DHH\\|DHX\\|helitron'";
$category{'total'}="[0-9]"; #grep any line with numbers
$category{'classified'}="'Unknown\\|unknown\\/unknow\\|repeat_region\\|Unspecified'"; #unknown TEs of all kind

die "The specified catetory $category is not found in our database!\n" unless exists $category{$category};

## get all classified regions
if ($category eq "classified"){
	`grep -v -P '$category{$category}' $std_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - | sort -suV > $std_out.$category.cvg`;
	`grep -v -P '$category{$category}' $tst_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - | sort -suV > $tst_out.$category.cvg` if $category eq "classified";
	} else {

## get stdlib and testlib cover range based on $category{$category}
`grep -i -P '$category{$category}' $tst_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - | sort -suV > $tst_out.$category.cvg`;
`grep -i -P '$category{$category}' $std_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7"\t"\$11}' - | sort -suV > $std_out.$category.cvg`;
if ($category eq "mite" or $category eq "tir"){
	`grep -i -v -P 'noTIR' $std_out.$category.cvg | sort -suV > $std_out.$category.cvg.temp`;
	`mv $std_out.$category.cvg.temp $std_out.$category.cvg`;
	}

## if $unknown == 1, then include the "unknown" type annotation into the user specified $category{$category}.
`perl -nle 's/^\\s+//g; my \$cat=(split)[10]; \$cat = lc \$cat; print \$_ if \$cat eq "unknown" or \$cat eq "unspecified"' $tst_out |awk '{print \$5"\t"\$6"\t"\$7}' - | sort -suV >> $tst_out.$category.cvg` if $unknown == 1;
	}

## bed arithmetics
`perl $script_path/util/combine_overlap.pl $std_out.$category.cvg $std_out.$category.cvg.cbi`;
`perl $script_path/util/combine_overlap.pl  $tst_out.$category.cvg $tst_out.$category.cvg.cbi`;

`perl $script_path/util/substract_parallel.pl $tst_out.$category.cvg.cbi $std_out.$category.cvg.cbi $threads`;
`perl $script_path/util/substract_parallel.pl $std_out.$category.cvg.cbi $tst_out.$category.cvg.cbi $threads`;
`perl $script_path/util/substract_parallel.pl $std_out.$category.cvg.cbi $std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi $threads`;
`cat $std_out.$category.cvg.cbi $tst_out.$category.cvg.cbi > ${std_out}_$tst_out.$category-cmb`;
`perl $script_path/util/combine_overlap.pl  ${std_out}_$tst_out.$category-cmb ${std_out}_$tst_out.$category-cmb.cbi`;
`perl $script_path/util/substract_parallel.pl $genome.list ${std_out}_$tst_out.$category-cmb.cbi $threads`;

## FP - false positive
my ($FP, $FN, $TP, $TN, $sens, $spec, $accu, $prec, $FDR, $F1) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
$FP=`perl $script_path/util/count_mask.pl $tst_out.$category.cvg.cbi-$std_out.$category.cvg.cbi`;
chomp $FP;

## FN - false negative
$FN=`perl $script_path/util/count_mask.pl $std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi`;
chomp $FN;

## TP - true positive
$TP=`perl $script_path/util/count_mask.pl $std_out.$category.cvg.cbi-$std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi`;
chomp $TP;

## TN - true negative
$TN=`perl $script_path/util/count_mask.pl $genome.list-${std_out}_$tst_out.$category-cmb.cbi`;
chomp $TN;

$sens=$TP/($TP+$FN) if $TP+$FN > 0;
$spec=$TN/($FP+$TN) if $FP+$TN > 0;
$accu=($TP+$TN)/($TP+$TN+$FP+$FN) if $TP+$TN+$FP+$FN > 0;
$prec=$TP/($TP+$FP) if $TP+$FP > 0;
$FDR=1-$prec;
$F1=2*$TP/(2*$TP+$FP+$FN) if $TP+$FP+$FN > 0;

## output results
open Out, ">$tst_out.$category.lib.report" or die $usage;
print Out "
Genome:			$genome
Standard annotation:	$std_out
Testing annotation:	$tst_out

sens=TP/(TP+FN)
spec=TN/(FP+TN)
accu=(TP+TN)/(TP+TN+FP+FN)
prec=TP/(TP+FP)
FDR=1-prec=FP/(TP+FP)
F1=2TP/(2TP+FP+FN)

TP: $TP\n\nFN: $FN\n\nTN: $TN\n\nFP: $FP\n\nSensitivity: $sens\nSpecificity: $spec\nAccuracy: $accu\nPrecision: $prec\nFDR: $FDR\nF1 measure: $F1\n

#Metrics\tsens\tspec\taccu\tprec\tFDR\tF1\tTP\tTN\tFP\tFN
$tst_out.$category.lib.report\t$sens\t$spec\t$accu\t$prec\t$FDR\t$F1\t$TP\t$TN\t$FP\t$FN
";
close Out;

## remove temporary files
`rm $genome $std_out $tst_out $genome.list $std_out.$category.cvg $tst_out.$category.cvg $std_out.$category.cvg.cbi $tst_out.$category.cvg.cbi $tst_out.$category.cvg.cbi-$std_out.$category.cvg.cbi $std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi $std_out.$category.cvg.cbi-$std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi ${std_out}_$tst_out.$category-cmb ${std_out}_$tst_out.$category-cmb.cbi $genome.list-${std_out}_$tst_out.$category-cmb.cbi` unless $debug == 1;
