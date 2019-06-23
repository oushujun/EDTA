#!/usr/bin/perl -w
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
#Update 05/09/2019
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
\n";

my $script_path = dirname(__FILE__);

my $genome='';
my $std_out='';
my $tst_out='';
my $category='total'; #TE category for testing, default total (all TEs)
my $unknown="0"; #0, not include unknown annotation in -tst, especially when -tst contains only 1 category (strict). 1 will include unknown..

my $includeN=0; #0 will not include N in total length of the genome, 1 will include N

my $k=0;
foreach (@ARGV){
	$genome=$ARGV[$k+1] if /^-genome$/i;
	$std_out=$ARGV[$k+1] if /^-std$/i;
	$tst_out=$ARGV[$k+1] if /^-tst$/i;
	$includeN=$ARGV[$k+1] if /^-N$/i;
	$category=lc $ARGV[$k+1] if /^-cat$/i;
	$unknown=$ARGV[$k+1] if /^-unknown$/i;
	$k++;
	}

## check files
die $usage unless -s $genome and -s $std_out and -s $tst_out;

## create softlinks to input files with unique ID in the current folder
my $rand=int(rand(1000000));
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
open Out, ">$tst_out.$category.lib.report" or die $usage;

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
$category{'ltr'}="'RLG\\|RLC\\|RLB\\|RLR\\|RLE\\|LTR\\|RLX\\|Gypsy\\|Copia'";
$category{'nonltr'}="'SINE\\|LINE\\|Penelope'";
$category{'line'}="LINE";
$category{'sine'}="SINE";
$category{'tir'}="'MITE\\|hAT\\|MULE\\|Tourist\\|CACT\\|MLE\\|PILE\\|POLE\\|EnSpm\\|MuDR\\|Stowaway\\|PIF\\|Harbinger\\|Tc1\\|En-Spm\\|PiggyBac\\|Mirage\\|P-element\\|Transib\\|DTA\\|DTH\\|DTT\\|DTM\\|DTC'";
$category{'mite'}="MITE";
$category{'helitron'}="Helitron";
$category{'total'}="[0-9]"; #grep any line with numbers
$category{'classified'}="'Unknown\\|unknown\\/unknown'"; #unknown TEs of all kind

die "The specified catetory $category is not found in our database!\n" unless exists $category{$category};

## get all classified regions
if ($category eq "classified"){
	`grep -v -P '$category{$category}' $std_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - > $std_out.$category.cvg`;
	`grep -v -P '$category{$category}' $tst_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - > $tst_out.$category.cvg` if $category eq "classified";
	} else {

## get stdlib and testlib cover range based on $category{$category}
`grep -i -P '$category{$category}' $std_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - > $std_out.$category.cvg`;
`grep -i -P '$category{$category}' $tst_out | awk '{if (\$6~/[0-9]+/) print \$5"\t"\$6"\t"\$7}' - > $tst_out.$category.cvg`;

## if $unknown == 1, then include the "unknown" type annotation into the user specified $category{$category}.
`perl -nle 's/^\\s+//g; my \$cat=(split)[10]; \$cat = lc \$cat; print \$_ if \$cat eq "unknown" or \$cat eq "unspecified"' $tst_out |awk '{print \$5"\t"\$6"\t"\$7}' - >> $tst_out.$category.cvg` if $unknown == 1;
	}

## bed arithmetics
`perl $script_path/util/combine_overlap.pl $std_out.$category.cvg $std_out.$category.cvg.cbi`;
`perl $script_path/util/combine_overlap.pl  $tst_out.$category.cvg $tst_out.$category.cvg.cbi`;

`perl $script_path/util/substract.pl $tst_out.$category.cvg.cbi $std_out.$category.cvg.cbi`;
`perl $script_path/util/substract.pl $std_out.$category.cvg.cbi $tst_out.$category.cvg.cbi`;
`perl $script_path/util/substract.pl $std_out.$category.cvg.cbi $std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi`;
`cat $std_out.$category.cvg.cbi $tst_out.$category.cvg.cbi > ${std_out}_$tst_out.$category-cmb`;
`perl $script_path/util/combine_overlap.pl  ${std_out}_$tst_out.$category-cmb ${std_out}_$tst_out.$category-cmb.cbi`;
`perl $script_path/util/substract.pl $genome.list ${std_out}_$tst_out.$category-cmb.cbi`;

## FP - false positive
my $FP=`perl $script_path/util/count_mask.pl $tst_out.$category.cvg.cbi-$std_out.$category.cvg.cbi`;
chomp $FP;

## FN - false negative
my $FN=`perl $script_path/util/count_mask.pl $std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi`;
chomp $FN;

## TP - true positive
my $TP=`perl $script_path/util/count_mask.pl $std_out.$category.cvg.cbi-$std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi`;
chomp $TP;

## TN - true negative
my $TN=`perl $script_path/util/count_mask.pl $genome.list-${std_out}_$tst_out.$category-cmb.cbi`;
chomp $TN;

my $sens=$TP/($TP+$FN);
my $spec=$TN/($FP+$TN);
my $accu=($TP+$TN)/($TP+$TN+$FP+$FN);
my $prec=$TP/($TP+$FP);
my $FDR=1-$prec;
my $F1=2*$TP/(2*$TP+$FP+$FN);

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
`rm $genome $std_out $tst_out $genome.list $std_out.$category.cvg $tst_out.$category.cvg $std_out.$category.cvg.cbi $tst_out.$category.cvg.cbi $tst_out.$category.cvg.cbi-$std_out.$category.cvg.cbi $std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi $std_out.$category.cvg.cbi-$std_out.$category.cvg.cbi-$tst_out.$category.cvg.cbi ${std_out}_$tst_out.$category-cmb ${std_out}_$tst_out.$category-cmb.cbi $genome.list-${std_out}_$tst_out.$category-cmb.cbi`;

