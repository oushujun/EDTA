#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com 04/22/2019)

my $usage = "\nFilter out GRF-main/mite candidates based on candidate, TIR, and TSD lengths. Takes FASTA and outputs a list of sequence names passed filtering.

	perl clean_GRF_TIR.pl GRF-mite.fasta > GRF-mite.fasta.clean.list

	You may change parameters in the script.
	You may use the script output_by_list.pl to extract sequences in the output list
	Input sequence names look like these:
		>10:484042:484071:5m1M4m:TAT (GRF-mite produced, \$TSDinfo = 1)
		>10:484042:484071:5m1M4m (GRF-main tir produced, \$TSDinfo = 0)
\n";



my $file = $ARGV[0];
my $minlen = 80; #min candidate seq length (bp)
my $mintir = 25; #combination of head and tail (bp)
my $TSDinfo = 1; #1 means sequence names carry TSD info; 0 means no TSD info

open File, "<$file" or die $usage;

my %cluster; #filter based on head (start) coordinates
while (<File>){
	next unless s/>//;
	chomp;
	my ($chr, $start, $end, $pairing, $tsd) = (split /:/, $_);

	#basic filtering
	my ($mutation, $iden) = (0, 0);
	$mutation += $1 while $pairing =~ s/([0-9]+)[MID]//;
	$iden += $1 while $pairing =~ s/([0-9]+)[m]//;
	next if $mutation + $iden < $mintir; #tir length filtering
	next if $end - $start + 1 < $minlen; #length filtering

	#deal with candidates with shared start, pick the one with less SNPs and Indels
	my $key = "$chr-$start";
	if (exists $cluster{$key}) {
		if ($cluster{$key}[4] = $mutation){
			$cluster{$key} = [$_, $chr, $start, $end, $mutation, $tsd] if $cluster{$key}[3] > $end; #keep the shorter candidate
			}
		elsif ($cluster{$key}[4] > $mutation){
			$cluster{$key} = [$_, $chr, $start, $end, $mutation, $tsd];
			}
		} else {
		$cluster{$key} = [$_, $chr, $start, $end, $mutation, $tsd];
		}

if ($TSDinfo == 1){
	#deal with candidates with +- 1bp shifted coordinates, pick the one with longer TSD
	my ($start_1p, $start_1m, $end_1p, $end_1m) = ($start+1, $start-1, $end+1, $end-1);
	#+1
	if (exists $cluster{"$chr-$start_1p"}){
		if (abs($end - $cluster{"$chr-$start_1p"}[3]) < 5){
			my $len = 0;
			$len = length $cluster{"$chr-$start_1p"}[5] if defined $cluster{"$chr-$start_1p"}[5];
			if ($len < length $tsd){
				delete $cluster{"$chr-$start_1p"}[5];
				} else {
				delete $cluster{$key};
				}
			}
		}
	#-1
	if (exists $cluster{"$chr-$start_1m"}){
		if (abs($end - $cluster{"$chr-$start_1m"}[3]) < 5){
			my $len = 0;
			$len = length $cluster{"$chr-$start_1m"}[5] if defined $cluster{"$chr-$start_1m"}[5];
			if ($len < length $tsd){
				delete $cluster{"$chr-$start_1m"}[5];
				} else {
				delete $cluster{$key};
				}
			}
		}
	}
	}
close File;

my %filter;  #filter based on taild (end) coordinates
foreach my $key (keys %cluster){
	my ($info, $chr, $start, $end, $mutation, $tsd) = @{$cluster{$key}};
	my $tail = "$chr-$end";
	if (exists $filter{$tail}) {
		if ($filter{$tail}[4] = $mutation){
			$filter{$tail} = [$info, $chr, $start, $end, $mutation, $tsd] if $filter{$tail}[2] < $start; #keep the shorter candidate
			}
		elsif ($filter{$tail}[4] > $mutation){ #keep the candidate with more identitcal TIR
			$filter{$tail} = [$info, $chr, $start, $end, $mutation, $tsd];
			}
		} else {
		$filter{$tail} = [$info, $chr, $start, $end, $mutation, $tsd];
		}
	}

#print out candidates
foreach my $tail (keys %filter){
	print "$filter{$tail}[0]\n";
	}

