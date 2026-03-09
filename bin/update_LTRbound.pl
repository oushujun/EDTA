#!/usr/bin/env perl
use warnings;
use strict;

## Update LTR boundary file with renamed TE IDs from rename_TE.pl mapping
## Usage: perl update_LTRbound.pl <rename_map> <LTRbound> > <updated_LTRbound>
## Input rename_map format: TE_N#class\toriginal_name#class
## Input LTRbound format: original_name#class\ttotal_len\tlLTR_len\trLTR_len

die "Usage: perl update_LTRbound.pl <rename_map> <LTRbound>\n" unless @ARGV == 2;

open MAP, "<$ARGV[0]" or die "Cannot open $ARGV[0]: $!\n";
open BOUND, "<$ARGV[1]" or die "Cannot open $ARGV[1]: $!\n";

# Build reverse map: original_name => TE_name (strip #class for matching)
my %map;
while (<MAP>){
	chomp;
	my ($new_name, $ori_name) = split /\t/;
	next unless defined $ori_name;
	# Strip #class from original name to match boundary file entries
	my $ori_base = $ori_name;
	$ori_base =~ s/#.*//;
	$map{$ori_base} = $new_name;
}
close MAP;

# Read boundary file and replace names
while (<BOUND>){
	chomp;
	my ($name, @rest) = split /\t/;
	# Strip #class for matching
	my $base = $name;
	$base =~ s/#.*//;
	if (defined $map{$base}){
		print "$map{$base}\t" . join("\t", @rest) . "\n";
	}
}
close BOUND;