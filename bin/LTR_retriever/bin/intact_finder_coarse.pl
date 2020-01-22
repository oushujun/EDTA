#!/usr/bin/env perl -w
use strict;

#Indentify intact LTRs in RM out files
#Shujun Ou (oushujun@msu.edu)
#07-25-2017

my $usage="perl this_script.pl RepeatMasker.out > intact_list";

my @arr;
while (<>){
	s/^\s+//g;
	next unless /^[0-9]+/;
	my ($score, $start, $end)=(split)[0,5,6];
	next if $score<=300 or $end-$start<=100; #require alignment score>300 or alignment length>100bp
	push @arr, $_;
	if (@arr>4){
		my ($div0, $chr0, $s0, $e0, $dir0, $loc0, $fam0)=(split /\s+/, $arr[0])[1,4,5,6,8,9,10];
		my ($div, $chr, $s, $e, $dir, $loc, $fam)=(split /\s+/, $arr[1])[1,4,5,6,8,9,10];
		my ($div1, $chr1, $s1, $e1, $dir1, $loc1, $fam1)=(split /\s+/, $arr[2])[1,4,5,6,8,9,10];
		my ($div2, $chr2, $s2, $e2, $dir2, $loc2, $fam2)=(split /\s+/, $arr[3])[1,4,5,6,8,9,10];
#		print "$arr[0]\n$arr[1]\n$arr[2]\t$arr[3]\n" if ($chr0 eq $chr1) and ($dir0 eq $dir1) and ($loc0 =~ /LTR/i) and ($loc1 =~ /LTR/i) and ($loc =~ /int/i) and  ($loc0 !~ /int/i) and ($loc1 !~ /int/i) and ($s-$e0<300) and ($s1-$e<300) and (abs($div1-$div0)<4); #LTR-int-LTR, LTRs can be from different family
#		print "$arr[0]\n$arr[1]\n$arr[2]\t$arr[3]\n" if ($chr0 eq $chr2) and ($dir0 eq $dir2) and ($loc0  =~ /LTR/i) and ($loc2=~ /LTR/i) and ($loc0 !~ /int/i) and ($loc2 !~ /int/i) and ($loc =~ /int/i) and ($loc1 =~ /int/i) and ($s-$e0<10000) and ($s1-$e<300) and ($s2-$e1<300) and (abs($div2-$div0)<4); #LTR-int-int-LTR, LTRs can be from different family
		print "$loc0\n" if ($chr0 eq $chr1) and ($dir0 eq $dir1) and ($loc0 =~ /LTR/i) and ($loc1 =~ /LTR/i) and ($loc =~ /int/i) and  ($loc0 !~ /int/i) and ($loc1 !~ /int/i) and ($s-$e0<300) and ($s1-$e<300) and (abs($div1-$div0)<4); #LTR-int-LTR, LTRs can be from different family
		print "$loc0\n" if ($chr0 eq $chr2) and ($dir0 eq $dir2) and ($loc0  =~ /LTR/i) and ($loc2=~ /LTR/i) and ($loc0 !~ /int/i) and ($loc2 !~ /int/i) and ($loc =~ /int/i) and ($loc1 =~ /int/i) and ($s-$e0<10000) and ($s1-$e<300) and ($s2-$e1<300) and (abs($div2-$div0)<4); #LTR-int-int-LTR, LTRs can be from different family
		shift @arr;
		}
	}
