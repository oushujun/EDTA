#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

##Generate bed file from RepeatMasker .out file
##Usage: perl RMout2bed.pl genome.fa.out > genome.fa.out.bed
##Author: Shujun Ou (shujun.ou.1@gmail.com) 07/12/2020

my $usage = "\n\tperl RMout2bed.pl RM.out > RM.out.bed\n\n";

open RMout, "sort -sV -k5,5 $ARGV[0]|" or die $usage;
print "#seqid start end TE_id repeat_class/superfamily method identity sw_score strand phase extra_info\n";

my $method = "homology";
while (<RMout>){
	s/^\s+//;
	my ($SW_score, $div, $iden, $chr, $chr_len, $element_start, $element_end, $element_length, $left_len, $strand, $TE_ID, $TE_class);
	($SW_score, $div, undef, undef, $chr, $element_start, $element_end, $left_len, $strand, $TE_ID, $TE_class) = (split);
	next unless defined $SW_score and $SW_score =~ /[0-9]+/;
	$element_length = $element_end - $element_start + 1;
	next if $SW_score < 300 and $element_length < 80;

	$TE_ID =~ s/_INT-int$/_INT/;
	$strand = "-" if $strand eq "C";
	$strand = "." if $strand eq "NA";
	if ($div =~ /^([0-9\.ex\-\*]+)$/){
		$iden = 1 - $div/100;
		} else {
		$iden = "NA";
		}
	print "$chr\t$element_start\t$element_end\t$TE_ID\t$TE_class\t$method\t$iden\t$SW_score\t$strand\t.\t.\n";
	}
close RMout;
