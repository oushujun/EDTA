##To annotate library with family and direction information
##Usage: perl annotate_lib.pl $index.scn.adj $index.LTRlib.fa
##Shujun Ou (oushujun@msu.edu) 06/23/2016


#!/usr/bin/env perl -w
use strict;

open List, "<$ARGV[0]" or die "ERROR: $!";
open Lib, "<$ARGV[1]" or die "ERROR: $!";

my %info;
while (<List>){
	next if /^#/;
	s/^\s+//;
	my ($start, $end, $lstart, $lend, $rstart, $rend, $direction, $fam)=(split)[0,1,3,4,6,7,12,18];
	my ($istart, $iend)=($lend+1, $rstart-1);
	next unless defined $fam;
	$info{"$start..$end"}=[$direction, $fam];
	$info{"$lstart..$lend"}=[$direction, $fam];
	$info{"$rstart..$rend"}=[$direction, $fam];
	$info{"$istart..$iend"}=[$direction, $fam];
	}

$/="\n>";
while (<Lib>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$id=~s/\s+//g;
	$seq=~s/\s+//g;
	my ($chr, $region)=("NA", "NA");
	($chr=$1, $id=$2, $region=$3) if $id=~/^(.*):([0-9]+\.\.[0-9]+)\|(.*)$/;
	$region="LTR" if $region=~/LTR_[1-2]/;
	$region="INT" if $region=~/LTR_IN/;
	my ($direction, $fam)=("?", "unknown");
	($direction, $fam)=@{$info{$id}}[0,1] if defined $info{$id};
	if ($direction eq "-"){
		$seq=~tr/ATGCatgc/TACGtacg/;
		$seq=reverse $seq;
		}

	print ">$chr:${id}_$region#LTR/$fam\n$seq\n";
	}
close Lib;
close List;
