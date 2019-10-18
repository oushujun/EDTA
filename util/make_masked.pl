#!/usr/bin/perl -w
use strict;

my $usage = "
Make hard/soft masked genome with filtering features.
Usage: perl make_masked.pl -genome unmasked_genome.fa -rmout repeatmasker.out [options]
		-maxdiv	[0-100]	Maximum divergence. Default: 40 (%)
		-minscore	[int]	Minimum SW score for a real match. Default: 300
		-minlen	[int]	Minimum alignment length. Default: 500 (bp)
		-hardmask	[0|1]	Output softmask (0) or hardmask (1, default) genome.
		-misschar	[chr]	Define the letter representing unknown sequences; default: N
\n";

# parameters
my $max_div = 40;
my $min_SW = 300;
my $min_len = 500;
my $hardmask = 1;
my $N = "N"; #the letter for hardmasked sequences
my $RMout = '';
my $genome = '';

# read parameters
my $k=0;
foreach (@ARGV){
	$max_div=$ARGV[$k+1] if /^-maxdiv$/i;
	$min_SW = $ARGV[$k+1] if /^-minscore$/i;
	$min_len=$ARGV[$k+1] if /^-minlen$/i;
	$hardmask=$ARGV[$k+1] if /^-hardmask$/i;
	$N=$ARGV[$k+1] if /^-misschar$/i;
	$RMout=$ARGV[$k+1] if /^-rmout$/i;
	$genome=$ARGV[$k+1] if /^-genome$/i;
	$k++;
	}

# read files
open Genome, "<$genome" or die $usage;
open RMout, "<$RMout" or die $usage;
open Masked, ">$genome.masked" or die $!;
open RMnew, ">$RMout.new" or die $!;

$/ = "\n>";
my %genome;
while (<Genome>){
	next if /^>\s?$/;
	chomp;
	s/>//g;
	s/^\s+//;
	my ($chr, $seq)=(split /\n/, $_, 2);
	$chr=~s/\s+$//;
	$seq=~s/\s+//g;
	$genome{$chr}=$seq;
	}
$/="\n";
close Genome;

while (<RMout>){
	s/^\s+//;
	next if /^\s+$/;
	next if /^$/;
	my ($SW_score, $div, $chr, $chr_len, $element_start, $element_end, $element_length, $left_len, $strand, $TE_ID, $TE_class);
	($SW_score, $div, $chr, $element_start, $element_end, $left_len, $strand, $TE_ID, $TE_class) = (split)[0,1,4,5,6,7,8,9,10];
	(print RMnew $_ and next) if $SW_score !~ /[0-9]+/;

	# filter short TEs
	$element_length = $element_end - $element_start + 1;
	next if $element_length < $min_len;

	# filter alignment scores and divergence
	next if $SW_score < $min_SW;
	next if $div > $max_div;
	$strand="-" if $strand eq "C";

	# mask the genome
	next unless exists $genome{$chr};
	substr($genome{$chr}, $element_start-1, $element_length) = $N x $element_length if $hardmask == 1;
	substr($genome{$chr}, $element_start-1, $element_length) = lc (substr($genome{$chr}, $element_start-1, $element_length)) if $hardmask == 0;

	# print new RMout line
	print RMnew $_;
	}
close RMnew;

foreach my $chr (sort {$a cmp $b} keys %genome){
	print Masked ">$chr\n$genome{$chr}\n";
	}


