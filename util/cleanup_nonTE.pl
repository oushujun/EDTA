#!/usr/bin/perl -w
use strict;

my $usage = "";

my $maxlen = 3000; #the max length allowed for non-repetitive sequence; longer than this will be removed
my $minlen = 100; #the min length to retain a seq, if shorter than this after removal it will be discarded
my $seq = $ARGV[0];

open FA, "<$seq" or die $usage;

$/ = "\n>";
while (<FA>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	$seq =~ s/([ATCGN]{$maxlen,})//g;
#	while ($seq =~ s/([ATCGN]{$maxlen,})//){my $len=length $1; print "$len\n";}
	$seq = uc $seq;
	my $len = length $seq;
	next if $len < $minlen;
	print ">$id\n$seq\n";
	}
close FA;


