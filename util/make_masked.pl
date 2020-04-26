#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;

my $usage = "
Make hard/soft masked genome with filtering features.
Usage: perl make_masked.pl -genome unmasked_genome.fa [options]
		-rmout	[file]	Required. The repeatmasker.out file
		-exclude	[file]	Optional. A list of regions in BED format to be excluded from TE annotation
		-maxdiv	[0-100]	Maximum divergence. Default: 40 (%)
		-minscore	[int]	Minimum SW score for a real match. Default: 300
		-minlen	[int]	Minimum alignment length. Default: 500 (bp)
		-hardmask	[0|1]	Output softmask (0) or hardmask (1, default) genome.
		-misschar	[chr]	Define the letter representing unknown sequences; default: N
		-threads|-t	[int]	Number of threads to run this program. Default: 4
\n";

# dependentcies
my $script_path = $FindBin::Bin;
my $substract = "$script_path/substract_parallel.pl";
my $combine = "$script_path/combine_overlap.pl";
die "The script substract_parallel.pl is not found in $substract!\n" unless -s $substract;
die "The script combine_overlap.pl is not found in $combine!\n" unless -s $combine;

# parameters
my $max_div = 40;
my $min_SW = 300;
my $min_len = 500;
my $hardmask = 1;
my $N = "N"; #the letter for hardmasked sequences
my $RMout = '';
my $genome = '';
my $exclude = '';
my $threads = 4;

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
	$exclude = $ARGV[$k+1] if /^-exclude$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i;
	$k++;
	}

# open file handles
open Genome, "<$genome" or die $usage;
open Masked, ">$genome.new.masked" or die $!;
open RMout, "<$RMout" or die $usage;
open RMnew, ">$RMout.new" or die $!;
die "The RepeatMasker.out file $RMout is empty or not exist!\n" unless -s $RMout;

# read genome files
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

# process the RM out file
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

	# print new RMout line
	print RMnew $_;
	}
close RMnew;

# get bed file from RMnew
`awk '{if (\$6~/[0-9]+/)print \$5"\t"\$6"\t"\$7"\t"\$10"#"\$11}' $RMout.new > $RMout.new.bed`;
`perl $combine $RMout.new.bed $RMout.new.bed.cbi`;
`cp $RMout.new.bed.cbi $RMout.target.bed`;

# remove exclude regions
if ($exclude ne ''){
	die "The exclude BED file $exclude is empty or not exist!\n" unless -s $exclude;
	`perl $combine $exclude $exclude.cbi`;
	`perl $substract $RMout.new.bed.cbi $exclude.cbi $threads`;
	`mv $RMout.new.bed.cbi-$exclude.cbi $RMout.target.bed`;
	`rm $RMout.new.bed $RMout.new.bed.cbi $exclude.cbi`;
	}

# mask the genome using $RMout.target.bed
open Mask, "<$RMout.target.bed" or die "The mask target is empty or not exist!\n";
while (<Mask>){
	my ($chr, $element_start, $element_end) = (split);
	my $element_length = $element_end - $element_start + 1;
	substr($genome{$chr}, $element_start-1, $element_length) = $N x $element_length if $hardmask == 1;
	substr($genome{$chr}, $element_start-1, $element_length) = lc (substr($genome{$chr}, $element_start-1, $element_length)) if $hardmask == 0;
	}
close Mask;

# print out masked genome
foreach my $chr (sort {$a cmp $b} keys %genome){
	print Masked ">$chr\n$genome{$chr}\n";
	}


