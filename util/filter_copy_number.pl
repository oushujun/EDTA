#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com; 05/18/2019)

my $usage = "\n\tperl filter_copy_number.pl file.fa\n\n";


my $seq = $ARGV[0]; #the seq file subjected to filter

my $mincopy = 2;
my $threads = 36;
my $set_cdhit="-c 0.8 -G 0.8 -s 0.8 -aL 0.8 -aS 0.8 -M 0"; #set parameters for cdhit
my $cdhitpath='';

#check cd-hit
$cdhitpath=`which cd-hit-est 2>/dev/null` if $cdhitpath eq '';
$cdhitpath=~s/cd-hit-est\n//;
die "cd-hit-est is not exist in the CDHIT path $cdhitpath!\n" unless -X "${cdhitpath}cd-hit-est";
die "\nThe fasta file is not found!\n$usage" unless -s $seq;

#run cd-hit for culstering
`${cdhitpath}cd-hit-est -i $seq -o $seq.clust $set_cdhit -T $threads`;
my $clust = "$seq.clust.clstr";

#filter clusters based on cluster size
open FA, "<$seq" or die $usage;
my %seq;
$/ = "\n>";
while (<FA>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my $key = ">$1" if $id =~ /^(.{19})/;
	$seq{$key} = ">$id\n$seq\n";
	}
close FA;

open Clust, "<$clust" or die $usage;
$/ = "\n>";
while (<Clust>){
	s/\n>$//;
	my @id = (split /\n/, $_);
	my $copy = @id - 1;
	next unless $copy >= $mincopy; 
	foreach (@id){
		next unless /^[0-9]+/;
		next unless /\*$/;
		my $key = (split)[2];
		$key =~ s/\.\.\.$//;
		my $pass = $seq{$key};
		print "$seq{$key}";
		}
	}
close Clust;


