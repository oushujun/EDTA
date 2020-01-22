#!/usr/bin/env perl -w
use strict;
#Find full length (100%)  exact matches (100%) to the TE library.
#usage: perl find_flTE.pl genome.fa.out > genome.fa.out.fl
#Shujun Ou (shujun.ou.1@gmail.com)

my $RMout = $ARGV[0];

while (<>){
	s/[\(\)]+//g;
	my ($SW, $div, $indel, $gap, $chr, $start, $end, $strand, $id, $type, $TEs, $TEe, $TEleft) = (split)[0,1,2,3,4,5,6,8,9,10,11,12,13];
	next if $type eq "Simple_repeat";
	next unless $SW =~ /[0-9]+/;
	next unless $div == 0 and $indel == 0 and $gap == 0;
	if ($strand eq "+"){
		next unless $TEs == 1 and $TEleft == 0;
		print $_;
		}
	if ($strand eq "C"){
		next unless $TEs == 0 and $TEleft == 1;
		print $_;
		}
	}

