#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com; 05/18/2019)

my $usage = "\n\tperl filter_copy_number.pl file.fa file.cd-hit.clstr\n\n";
my $mincopy = 2;

my $fasta = $ARGV[0]; #the seq file subjected to filter
my $clust = $ARGV[1]; #a cd-hit-est produced *.clstr file

die "\nThe fasta file is not found!\n$usage" unless -s $fasta;
die "\nThe clstr file is not found!\n$usage" unless -s $clust;

open FA, "<$fasta" or die $usage;
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


