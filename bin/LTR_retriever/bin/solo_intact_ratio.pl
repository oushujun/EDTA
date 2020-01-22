#!/usr/bin/env perl -w
use strict;

#compute solo-intact LTR ratio for given lists
#Shujun Ou (oushujun@msu.edu)
#07-25-2017
#update: 04/16/2019

my $usage="\n\tperl solo_intact_ratio.pl solo_list intact_list > solo_intact_ratio\n
	To make solo and intact lists:
		perl solo_finder.pl RepeatMasker.out > solo_list
		perl intact_finder_coarse.pl RepeatMasker.out > intact_list\n\n";

open Solo, "<$ARGV[0]" or die $usage;
open Intact, "<$ARGV[1]" or die $usage;

my %all;

while (<Intact>){
	chomp;
	s/^\s+//;
	my $id=(split)[0];
	if (exists $all{$id}){
		$all{$id}[1]++; #count intact number
		} else {
		$all{$id}=["0", "0"]; #initialize count of solo and intact number
		}
	}
close Intact;

while (<Solo>){
	chomp;
	s/^\s+//;
	my $id=(split)[0];
	if (exists $all{$id}){
		$all{$id}[0]++; #count solo number
		} else {
		$all{$id}=["0", "0"];
		}
	}
close Solo;

print "LTR_fam\tSolo_count\tIntact_count\tSolo-Intact-ratio\n";
foreach my $id (sort {$a cmp $b} keys %all){
	my $ratio="NA";
	my ($solo, $intact)=($all{$id}[0], $all{$id}[1]);
	if ($intact==0){
		$ratio="inf"; #no intact LTR found for this family
		}
	elsif ($solo==0){
		$ratio=0;
		}
	else {
		$ratio=sprintf ("%.1f", $solo/$intact);
		}
	print "$id\t$solo\t$intact\t$ratio\n";
	}

