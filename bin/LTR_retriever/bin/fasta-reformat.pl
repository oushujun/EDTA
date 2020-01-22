#!/usr/bin/env perl -w
use strict;
#contributed by Ning Jiang (jiangn@msu.edu)

my $usage = "This script is written to convert fasta files into a prettier format. 
Usage: fasta-reformat.pl input-fasta-file number-of-positions-per-line\n";

if (@ARGV < 2) {die "ERROR: $usage";}
if ($ARGV[1] < 1) {die "ERROR: $usage";}

open(FA, "<$ARGV[0]") || die "ERROR: $usage";

my $seq = "";
while (<FA>) {
    if (/>\s*(.+)/) {
	if ($seq) {
	    my @sym = split(//, $seq);
	    my $ct = 0;
	    foreach my $sym (@sym) {
		print $sym;
		$ct ++;
		if ( !($ct%$ARGV[1]) ) {print "\n";}
	    }
	    if ($ct%$ARGV[1]) {print "\n";}
	}
	printf ">%s\n", $1;
	$seq = "";
    } else {
	chomp;
	$seq .= $_;
    }
}
close FA;

my @sym = split(//, $seq);
my $ct = 0;
foreach my $sym (@sym) {
    print $sym;
    $ct ++;
    if ( !($ct%$ARGV[1]) ) {print "\n";}
}
if ($ct%$ARGV[1]) {print "\n";}
