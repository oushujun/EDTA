#!/usr/bin/env perl -w
use strict;

#Indentify solo LTRs in RM out files
#Shujun Ou (oushujun@msu.edu)
#07-25-2017

my $usage="perl this_script.pl RepeatMasker.out > solo_list";

my @arr;
while (<>){
	s/^\s+//g;
	next unless /^[0-9]+/;
	my ($score, $start, $end)=(split)[0,5,6];
	next if $score<=300 or $end-$start<=100; #require alignment score>300 or alignment length>100bp
	push @arr, $_;
	if (@arr>5){
		my ($div2, $chr2, $s2, $e2, $dir2, $loc2, $fam2)=(split /\s+/, $arr[2])[1,4,5,6,8,9,10];
		my $loc2_len=abs($2-$1) if $loc2=~/.*:([0-9]+)\.\.([0-9]+)_.*/;
		my $keep=1;
		$keep=0 if $loc2=~/int/i; #skip this line if encounter internal regions
		$keep=0 if abs($e2-$s2)<$loc2_len*0.8; #discard this entry if the length is less than 80% of the hit length
		foreach ($arr[0], $arr[1], $arr[3], $arr[4]){
			my ($div0, $chr0, $s0, $e0, $dir0, $loc0, $fam0)=(split /\s+/, $_)[1,4,5,6,8,9,10];
			$keep=0 if ($chr0 eq $chr2) and ($dir0 eq $dir2) and ($loc0 eq $loc2) and (abs($div2-$div0)<4); #require difference of diversity larger than 4% compared to other LTRs with the same ID matched to be a solo-LTR
			$keep=0 if ($chr0 eq $chr2) and ($loc0 =~ /int/i) and (abs($s2-$e0)<300) and (abs($s0-$e2)<300); #if there is int exists in this working window, require it has at least 300bp distance from a solo-LTR
			}
#		print "$loc2\n" if $keep==1; #print the library entry out
#		print "$chr2:$s2..$e2\t$loc2\n" if $keep==1;  #print the actual locus of the solo LTR
		print "$loc2\t$chr2:$s2..$e2\n" if $keep==1;  #print the actual locus of the solo LTR
		shift @arr;
		}
	}
