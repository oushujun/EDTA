##Generate gff3 file from RepeatMasker .out file of LTR_retriever
##Usage: perl make_gff3.pl genome.fa.out
##Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University
##Version: 1.0 01-27-2018

#!/usr/bin/env perl -w
use strict;

my $date=`date -u`;
chomp ($date);
open RMout, "sort -sV -k5,5 $ARGV[0]|" or die "ERROR: $!";
open GFF, ">$ARGV[0].gff" or die "ERROR: $!";
print GFF "##gff-version 2\n##date $date
##Chromosome Annotator Repeat_class/superfamily Start End Diversity(%) Strand SW_score Repeat_famliy\n";

my %seq_flag;
my $annotator="RepeatMasker";
while (<RMout>){
	s/^\s+//;
	my ($SW_score, $div, $chr, $chr_len, $element_start, $element_end, $element_length, $left_len, $strand, $TE_ID, $TE_class);
	($SW_score, $div, $chr, $element_start, $element_end, $left_len, $strand, $TE_ID, $TE_class)=(split)[0,1,4,5,6,7,8,9,10];
	next if $SW_score!~/[0-9]+/;
	$TE_ID=~s/_INT-int$/_INT/;
	$element_length=$element_end-$element_start+1;
	next if $SW_score<300 and $element_length<80;
	$strand="-" if $strand eq "C";
	unless (exists $seq_flag{$chr}){
		$seq_flag{$chr}=$chr;
		$left_len=~s/[\(\) ]+//g;
		$chr_len=$element_end+$left_len;
		print GFF "##sequence-region $chr 1 $chr_len\n";
		}
	print GFF "$chr\t$annotator\t$TE_class\t$element_start\t$element_end\t$div\t$strand\t$SW_score\t$TE_ID\n"
	}
close RMout;
close GFF;
