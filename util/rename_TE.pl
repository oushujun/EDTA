#!/usr/bin/perl -w
use strict;

my $usage = "";
my $fasta = $ARGV[0];

open FA, "<$fasta" or die "\nInput not found!\n$usage";
$/ = "\n>";
my $num = 0;
while (<FA>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	my $name = (split /\s+/, $id)[0];
	$seq =~ s/\s+//g;
	my $class = $1 if $name =~ /#(.*)$/;
	$num = sprintf("%08d", $num);
	print ">TE_$num#$class\n$seq\n";
	$num++;
	}
close FA;


