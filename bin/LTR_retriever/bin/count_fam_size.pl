#!/usr/bin/env perl -w
use strict;

#Count the size of each LTR family. A family is determined by the LTR region.
#Usage: perl count_fam_size.pl LTR_retriever.pass.list RepeatMasker.out > LTR_retriever.family.size
#Author: Shujun Ou (oushujun@msu.edu) 08/16/2017

open List, "<$ARGV[0]" or die $!;
open Anno, "<$ARGV[1]" or die $!;

my %list;
my %translate;
while (<List>){
	next if /^#/;
	s/^\s+//g;
	my ($id, $in, $family, $age)=(split)[0,6,9,11];
	my ($chr, $start, $end)=($1, $2, $3) if $id=~/(.*):([0-9]+)\.\.([0-9]+)/;
	next unless defined $chr;
	$in=~s/IN://i;
	$age=4000000 if $age eq "NA"; #make NA the oldest element
	$age=$age/1000000; #convert to MY
	$translate{"$chr:$start"}=$id;
	$translate{"$chr:$end"}=$id;
	$translate{"$chr:$in"}=$id;
	$list{$id}=["0", $family, $age]; #size, family, age of the library entry
	}
close List;

while (<Anno>){
	s/^\s+//g;
	next unless /^[0-9]+/;
	my ($score, $start, $end, $id)=(split)[0,5,6,9];
	my $length=$end-$start+1;
	next if $score<300 or $length<100; #require SW alignment score>=300 or alignment length>=100bp
	my ($chr, $from, $to)=($1, $2, $3) if $id=~/(.*):([0-9]+)\.\.([0-9]+)/;
	my ($start_id, $end_id, $in_id)=("$chr:$from", "$chr:$to", "$chr:$from..$to"); #construct possible ids to search in %list
	foreach my $query ($start_id, $end_id, $in_id){
		next unless exists $translate{$query};
		my $id=$translate{$query};
		$list{$id}[0]+=$length;
		last;
		}
	}
close Anno;

print "#ID\tFam_size\tSuperfamily\tAge(MY)\n";
foreach my $id (sort { $list{$b}[0] <=> $list{$a}[0] } keys %list){ #sort value from large to small
	print "$id\t$list{$id}[0]\t$list{$id}[1]\t$list{$id}[2]\n" if $list{$id}[0]>0; #print out families that have size>0
	}


