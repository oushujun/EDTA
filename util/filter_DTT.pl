#!/usr/bin/perl -w
use strict;

#filter DTT candidates with CT***TC**T…A**GA***AG terminal motifs
$/ = "\n>";
while (<>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	if ($id =~ /DTT/){
		print ">$id\n$seq\n" if $seq =~ /^CT...TC..T.*A..GA...AG$/;
		} else {
		print ">$id\n$seq\n";
		}
	}
