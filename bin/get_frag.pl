#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;

#function: For entries in the subtrahend.list that are overlapping with entries in the minuend.list, do;
#		1. Discard the minuend entries that are enclosed in subtrahend regions.
#		    This is handeled by the keep_nest.pl script.
#		2. Remove the region in the minuend entry that is overlapping with the subtrahend entry
#		3. Keep the minuend region that is not overlapping with the subtrahend region
#usage: Modified from substract_parallel.pl
#	perl get_frag.pl minuend.list subtrahend.list thread_num
#Author: Shujun Ou (shujun.ou.1@gmail.com), 12/19/2019

my $usage = "\n\tperl get_frag.pl EDTA.RM.bed EDTA.intact.bed.cmb thread_num\n\n";
my $minlen = 50; #fragments shorter than this will be discarded

## read thread number
my $threads = 4;
if (defined $ARGV[2]){
	$threads = $ARGV[2];
	}

## minuend − subtrahend = difference
open Minuend, "sort -suV $ARGV[0] |" or die $usage;
open Subtrahend, "sort -suV $ARGV[1] |" or die $usage;
open Diff, ">$ARGV[0]-$ARGV[1]" or die $!;

my %substr;
while (<Subtrahend>){
	next if /^\s+$/;
	next if /^#/;
	chomp;
	my ($chr, $from, $to, $type)=(split)[0,1,2,11];
	push @{$substr{$chr}}, [$from, $to, $type, $_];
	}

## multi-threading using queue, put candidate regions into queue for parallel computation
my %diff :shared;
my $queue = Thread::Queue -> new();
while (<Minuend>){
	next if /^\s+$/;
	next if /^#/;
	chomp;
	my ($chr, $from, $to, $type)=(split)[0,1,2,11];
	next unless defined $chr;
	$queue->enqueue([$chr, $from, $to, $type, $_]);
	}
$queue -> end();
close Minuend;

## initiate a number of worker threads and run
foreach (1..$threads){
	threads -> create(\&substract);
	}
foreach (threads -> list()){
	$_ -> join();
	}

## output results
foreach my $id (sort {$a cmp $b} keys %diff){
	chomp $diff{$id};
	print Diff "$diff{$id}\n";
	}
close Diff;

## subrotine to perform substraction
sub substract(){
	while (defined ($_ = $queue->dequeue())){
	my $keep=1;
	my ($chr, $from, $to, $type, $entry) = (@{$_}[0], @{$_}[1], @{$_}[2], @{$_}[3], @{$_}[4]);
	my $info = (split /\s+/, $entry, 4)[3]; #TE info without coordinates
	Run:
	foreach my $substr (@{$substr{$chr}}){
		my @range=@{$substr}; #[$from, $to, $type, $_]
		# skip this $substr range when its on the left side of $from, $to
		next if $range[1]<$from;

		# end the loop when $substr range is on the right side of $from, $to
		last if $range[0]>$to;

		# if the $substr range is covering the entire [$from, $to], discard this 
		# [$from, $to] region (this region is taken care of in the keep_nest.pl script)
		if ($range[0]<=$from and $range[1]>=$to){
			$keep=0;
			last;
			}

		# if the $substr range is enclosed in [$from, $to], keep the difference of [$from, $to]
		if ($range[0]>=$from and $range[1]<=$to){
			$range[0]--;
			$diff{"$chr:$from:$range[0]"} = "$chr\t$from\t$range[0]\t$info" if $range[0]-$from+1 >= $minlen;
			$from=$range[1]+1;
			$keep=1;
			goto Run;
			}

		# keep the region of [$from, $$range[0]-1] when this $substr range is overlapping 
		# on the right side of [$from, $to]
		if ($range[1]>$to){
			$keep=0;
			$range[0]--;
			$diff{"$chr:$from:$range[0]"} = "$chr\t$from\t$range[0]\t$info" if $range[0]-$from+1 >= $minlen;
			last;
			}

		# change [$from, $to] to  [$range[1]+1, $to] when this $substr range is overlapping
		# on the left side of [$from, $to]
		if ($range[0]<$from){
			$from=$range[1]+1;
			$keep=1;
			goto Run;
			}
		}
	$diff{"$chr:$from:$to"} = "$chr\t$from\t$to\t$info" if $keep==1 and $to-$from+1 >= $minlen;
	$keep=1;
	}
	}


