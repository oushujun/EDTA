#!/usr/bin/env perl
use warnings;
use strict;

#filter DTT candidates with CT***TC**T…A**GA***AG terminal motifs
$/ = "\n>";
while (<>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	if ($id =~ /DTT/){
	#	print ">$id\n$seq\n" if $seq =~ /^CT...TC..T.*A..GA...AG$/;
	#	print ">$id\n$seq\n" if $seq =~ /^C[TA].[CT]CT.[CT].*[GA].AG[GA].[AT]G$/;
		print ">$id\n$seq\n" if $seq =~ /^C[TA].[CT]CT.[CT]/;
		} 
	elsif ($id =~ /DTM/){
		print ">$id\n$seq\n" if $seq =~ /^[GC][GA]/;
		}
	elsif ($id =~ /DTC/){
		print ">$id\n$seq\n" if $seq =~ /^CACT[AG]....A/;
		}
	else {
		print ">$id\n$seq\n";
		}
	}
