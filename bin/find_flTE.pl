#!/usr/bin/env perl
use warnings;
use strict;
#Find full length (100%)  exact matches (100%) to the TE library.
#usage: perl find_flTE.pl genome.fa.out > genome.fa.out.fl
#Shujun Ou (shujun.ou.1@gmail.com)

my $RMout = $ARGV[0];

my $stringent = 0; #0, use flexible parameters below; 1, requires full length and exact match.
my $max_div = 20;
my $max_ins = 10;
my $max_del = 10;
my $min_cov = 0.8;

while (<>){
	s/[\(\)]+//g;
	next if /^\s+?$/;
	my ($SW, $div, $ins, $del, $chr, $start, $end, $strand, $id, $type, $TEs, $TEe, $TEleft) = (split)[0,1,2,3,4,5,6,8,9,10,11,12,13];
	next if $type eq "Simple_repeat";
	next unless $SW =~ /[0-9]+/;
	if ($stringent == 1){
		next unless $div == 0 and $ins == 0 and $del == 0;
		if ($strand eq "+"){
			next unless $TEs == 1 and $TEleft == 0;
			print $_;
			}
		if ($strand eq "C"){
			next unless $TEs == 0 and $TEleft == 1;
			print $_;
			}
		} else {
		next unless $div <= $max_div and $ins <= $max_ins and $del <= $max_del;
		my ($full_len, $len) = (0, 0);
		if ($strand eq "+"){
			$full_len = $TEe + $TEleft;
			$len = $TEe - $TEs + 1;
			}
		if ($strand eq "C"){
			$full_len = $TEs + $TEe;
			$len = $TEe - $TEleft + 1;
			}
		next unless $len/($full_len+1) >= $min_cov;
		print $_;
		}
	}

