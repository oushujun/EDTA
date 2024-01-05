#!/usr/bin/env perl
use strict;
use warnings;
#usage: This script helps to find TEs that are consistently annotated in a genome
#Shujun Ou (shujun.ou.1@gmail.com) 04/25/2021
#usage:
#	perl ~/bin/EDTA/util/evaluation.pl -anno genome.RM.out -maxcount 0 -threads 50 -genome genome.fa 
#	perl ~/bin/EDTA/util/find_flTE.pl genome.RM.out | perl get_consistent_TE.pl - genome.RM.out.TE.fa.stat|sort -k2,2 -nr|less -S

my $rmout = $ARGV[0]; #the original RepeatMasker out file can be preprocess by find_flTE.pl to find full-length TEs
my $stat = $ARGV[1];

open Out, "<$rmout" or die $!;
open STAT, "<$stat" or die $!;

my %anno;
while (<Out>){
	my ($chr, $str, $end, $id, $supfam) = (split)[4,5,6,9,10];
	$anno{"$chr:$str..$end"} = "$id#$supfam";
	}
close Out;

my %consistent;
while (<STAT>){
	next unless /Cleaned/;
	s/;//g;
	my ($id1, $id2) = (split)[0,7];
	$id1 =~ s/\|.*//;
	$id2 =~ s/\|.*//;
	next unless defined $anno{$id1} and defined $anno{$id2};
	#print "$anno{$id1}\t$id1\t$id2\n" if $anno{$id1} eq $anno{$id2}; #test
	push @{$consistent{$anno{$id1}}}, $id1;
	push @{$consistent{$anno{$id1}}}, $id2;
	}
close STAT;

foreach my $id (keys %consistent){
	my @ids = @{$consistent{$id}};
	my @ids_uniq = do { my %seen; grep { !$seen{$_}++ } @ids };
	my $count = @ids_uniq;
	print "$id\t$count\t@ids_uniq\n";
	}
