#!/usr/bin/perl -w
use strict;

#To test the specificity and sensitivity of construct LTR library
#Warning: Run the program in the same directory at the same time may cause errors
#
#sens = P(test+|LTR+) = TP/(TP+FN)
#spec = P(test-|LTR-) = TN/(FP+TN)
#accu = (TP+TN)/(TP+TN+FP+FN)
#prec = TP/(TP+FP)
#
#-----------------------------  Whole genome
#    [=============]-----       Test lib
#----[=============]	 [    ] Standard lib
# FN      TP	     FP	   TN
#
#Author: Shujun Ou (oushujun@msu.edu), 03/08/2015
#
my $usage="	perl lib-test.pl all.fasta all.stdlib.out all.testlib.out\n";

my $genome=$ARGV[0];
my $std_out=$ARGV[1];
my $tst_out=$ARGV[2];
#my $script_path="/mnt/home/oushujun/git_bin/lib-testers";
my $script_path=`readlink -fn -- $0`;
$script_path=~s/(.+)\/.+$/$1/;

my $includeN=0; #0 will not include N in total length of the genome, 1 will include N

## get whole genome range
open FASTA, "<$genome" or die $usage;
open All, ">$genome.list" or die $usage;
open Out, ">$tst_out.lib.report" or die $usage;

if (1){ #test
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

## get stdlib and testlib cover range
`grep -i 'RLG\\|RLC\\|RLB\\|RLR\\|RLE\\|LTR\\|RLX\\|Gypsy\\|Copia' $std_out |awk '{if (\$6~/[0-9]+/)print \$5"\t"\$6"\t"\$7}' - > $std_out.LTR.cvg`;
`awk '{if (\$6~/[0-9]+/)print \$5"\t"\$6"\t"\$7}' $tst_out > $tst_out.cvg`;
#`awk '{if (\$6~/[0-9]+/)print \$5"\t"\$6"\t"\$7}' $tst_out |perl -ne 's/(.+)_([0-9]+)/\$1/;print \$_' > $tst_out.cvg`;
`perl $script_path/combine_overlap.pl $std_out.LTR.cvg $std_out.LTR.cvg.cbi`;
`perl $script_path/combine_overlap.pl  $tst_out.cvg $tst_out.cvg.cbi`;

`perl $script_path/substract.pl $tst_out.cvg.cbi $std_out.LTR.cvg.cbi`;
`perl $script_path/substract.pl $std_out.LTR.cvg.cbi $tst_out.cvg.cbi`;
`perl $script_path/substract.pl $std_out.LTR.cvg.cbi $std_out.LTR.cvg.cbi-$tst_out.cvg.cbi`;
`cat $std_out.LTR.cvg.cbi $tst_out.cvg.cbi > ${std_out}_$tst_out-cmb`;
`perl $script_path/combine_overlap.pl  ${std_out}_$tst_out-cmb ${std_out}_$tst_out-cmb.cbi`;
`perl $script_path/substract.pl $genome.list ${std_out}_$tst_out-cmb.cbi`;
} #test

#
## FP - false positive
my $FP=`perl $script_path/count_mask.pl $tst_out.cvg.cbi-$std_out.LTR.cvg.cbi`;

## FN - false negative
my $FN=`perl $script_path/count_mask.pl $std_out.LTR.cvg.cbi-$tst_out.cvg.cbi`;

## TP - true positive
my $TP=`perl $script_path/count_mask.pl $std_out.LTR.cvg.cbi-$std_out.LTR.cvg.cbi-$tst_out.cvg.cbi`;

## TN - true negative
my $TN=`perl $script_path/count_mask.pl $genome.list-${std_out}_$tst_out-cmb.cbi`;

my $sens=$TP/($TP+$FN);
my $spec=$TN/($FP+$TN);
my $accu=($TP+$TN)/($TP+$TN+$FP+$FN);
my $prec=$TP/($TP+$FP);

print Out "
sens=TP/(TP+FN)
spec=TN/(FP+TN)
accu=(TP+TN)/(TP+TN+FP+FN)
prec=TP/(TP+FP)

TP: $TP\nFN: $FN\nTN: $TN\nFP: $FP\nSensitivity: $sens\nSpecificity: $spec\nAccuracy: $accu\nPrecision: $prec\n";
close Out;

## remove temporary files
`rm $genome.list $std_out.LTR.cvg $tst_out.cvg $std_out.LTR.cvg.cbi $tst_out.cvg.cbi $tst_out.cvg.cbi-$std_out.LTR.cvg.cbi $std_out.LTR.cvg.cbi-$tst_out.cvg.cbi $std_out.LTR.cvg.cbi-$std_out.LTR.cvg.cbi-$tst_out.cvg.cbi ${std_out}_$tst_out-cmb ${std_out}_$tst_out-cmb.cbi $genome.list-${std_out}_$tst_out-cmb.cbi`;

