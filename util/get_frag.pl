#!/usr/bin/env perl -w
use strict;
use threads;
use Thread::Queue;
use threads::shared;

#usage: Modified from substract_parallel.pl
#Author: Shujun Ou (shujun.ou.1@gmail.com), 12/19/2019

my $usage = "\n\tperl get_frag.pl EDTA.TEanno.bed EDTA.intact.bed thread_num\n\n";
my $minlen = 50; #fragments shorter than this will be discarded

## read thread number
my $threads = 4;
if (defined $ARGV[2]){
	$threads = $ARGV[2];
	}

## minuend − subtrahend = difference
open Minuend, "<$ARGV[0]" or die $usage;
open Subtrahend, "<$ARGV[1]" or die $usage;
open Diff, ">$ARGV[0]-$ARGV[1]" or die $!;

my %substr;
while (<Subtrahend>){
	next if /^\s+$/;
	my ($chr, $from, $to, $type, $info)=(split /\s+/, $_, 5);
	push @{$substr{$chr}}, [$from, $to, $type, $info];
	}

## multi-threading using queue, put candidate regions into queue for parallel computation
my %diff :shared;
my $queue = Thread::Queue -> new();
while (<Minuend>){
	next if /^\s+$/;
	my ($chr, $from, $to, $type, $info)=(split /\s+/, $_, 5);
	next unless defined $chr;
	$queue->enqueue([$chr, $from, $to, $type, $info]);
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
	my ($chr, $from, $to, $type, $anno) = (@{$_}[0], @{$_}[1], @{$_}[2], @{$_}[3], @{$_}[4]);
	Run:
	foreach my $info (@{$substr{$chr}}){
		my @range=@{$info};
		last if $range[0]>$to;
		next if $range[1]<$from;
		$keep=0 if ($range[0]<=$from and $range[1]>=$to);
		if ($range[0]>$from){
			$keep=0;
			$range[0]--;
			$diff{"$chr:$from:$range[0]"} = "$chr\t$from\t$range[0]\t$type\t$anno" if $range[0]-$from+1 >= $minlen;
			}
		if ($range[1]<$to){
			$from=$range[1]+1;
			$keep=1;
			goto Run;
			}
		}
	$diff{"$chr:$from:$to"} = "$chr\t$from\t$to\t$type\t$anno" if $keep==1 and $to-$from+1 >= $minlen;
	$keep=1;
	}
	}


