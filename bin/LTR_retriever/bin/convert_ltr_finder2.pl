#!/usr/bin/env perl -w
use strict;
#convert LTR_FINDER -w 2 format to LTRharvest format
#Shujun Ou (shujun.ou.1@gmail.com) 11/13/2019
#Usage: perl this_script.pl LTR_FINDER.w2.scn > LTRharvest.scn

my $head ="## LTR_FINDER
## predictions are reported in the following way
## s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
## where:
## s = starting position
## e = ending position
## l = length
## ret = LTR-retrotransposon
## lLTR = left LTR
## rLTR = right LTR
## sim = similarity
## seq-nr = sequence number
";

print $head;

my $seq_id=-1;
while (<>){
	$seq_id++ if /^>/;
	next unless /^\[/;
	s/\[\s+/\[/g;
	my ($from, $to, $chr, $LTR_len, $lLTR_len, $rLTR_len, $lLTR_end, $rLTR_str, $similarity, $TSD, $motif, $direction, $loc, $len);
	$from=$to=$LTR_len=$from=$lLTR_end=$lLTR_len=$rLTR_str=$to=$rLTR_len=$similarity=$chr=$direction=$TSD=$motif=$loc=$len=' ';
	($chr, $loc, $len, $LTR_len, $direction, $similarity) = (split)[1,2,3,4,12,15];
	($from, $to) = ($1, $2) if $loc =~ /^([0-9]+)\-([0-9]+)$/;
	($lLTR_len, $rLTR_len) = ($1, $2) if $len =~ /^([0-9]+),([0-9]+)$/;
	($lLTR_end, $rLTR_str) = ($from+$lLTR_len-1, $to-$rLTR_len+1);

	print "$from $to $LTR_len $from $lLTR_end $lLTR_len $rLTR_str $to $rLTR_len $similarity $seq_id $chr\n";
	$from=$to=$LTR_len=$lLTR_len=$rLTR_len=$lLTR_end=$rLTR_str=$similarity=$TSD=$motif=$direction=$TSD=$motif=$loc=$len="NA";
	}


