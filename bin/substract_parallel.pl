#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;

#usage: perl substract_parallel.pl minuend.list subtrahend.list thread_num
#Author: Shujun Ou (oushujun@msu.edu), 08/02/2019
# 05/30/2023, ChatGPT v3.5

my $usage = "\n\tperl substract_parallel.pl minuend.list subtrahend.list thread_num\n\n";
die $usage unless @ARGV >= 2;

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
	my ($chr, $from, $to)=(split)[0,1,2];
	push @{$substr{$chr}}, [$from, $to];
	}

## multi-threading using queue, put candidate regions into queue for parallel computation
my %diff :shared;
my $queue = Thread::Queue -> new();
while (<Minuend>){
	next if /^\s+$/;
	my ($chr, $from, $to)=(split)[0,1,2];
	next unless defined $chr;
	$queue->enqueue([$chr, $from, $to]);
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
	my ($chr, $from, $to) = (split /:/, $id);
	print Diff "$chr\t$from\t$to\n"
	}
close Diff;

## subrotine to perform substraction
sub substract(){
	while (defined ($_ = $queue->dequeue())){
	my $keep=1;
	my ($chr, $from, $to) = (@{$_}[0], @{$_}[1], @{$_}[2]);

	if (exists $substr{$chr}) {
	foreach my $info (@{$substr{$chr}}){
		my ($start, $end) = @{$info};
		last if $start > $to;
		next if $end < $from;
		$keep = 0 if ($start <= $from && $end >= $to);
		if ($start > $from) {
			$keep=0;
			$start--;
			$diff{"$chr:$from:$start"} = "$chr:$from:$start";
			} # if $range[0]>$from;
		if ($end < $to) {
			$from = $end + 1;
			$keep=1;
			}
		}
	}
	lock(%diff);
	$diff{"$chr:$from:$to"} = "$chr:$from:$to" if $keep==1;
	$keep=1;
	}
	}


