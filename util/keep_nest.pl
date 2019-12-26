#!/usr/bin/perl -w
use strict;
use threads;
use Thread::Queue;
use threads::shared;

#usage: perl substract_parallel.pl minuend.list subtrahend.list thread_num
#Author: Shujun Ou (oushujun@msu.edu), 08/02/2019

my $usage = "\n\tperl substract_parallel.pl minuend.list subtrahend.list thread_num\n\n";

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
	my ($chr, $from, $to, $type, $info)=(split /\s+/, $_, 5);
#print "$chr, $from, $to, $type, $info\n";
	push @{$substr{$chr}}, [$from, $to, $type, $info];
	}

## multi-threading using queue, put candidate regions into queue for parallel computation
my %diff :shared;
my $queue = Thread::Queue -> new();
while (<Minuend>){
	next if /^\s+$/;
	my ($chr, $from, $to, $type, $info)=(split /\s+/, $_, 5);
	next unless defined $chr;
	$diff{"$chr:$from:$to"} = $_; #all minuend info are retained
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
	my ($chr, $from, $to, $type, $info) = (@{$_}[0], @{$_}[1], @{$_}[2], @{$_}[3], @{$_}[4]);
	foreach my $info (@{$substr{$chr}}){
		my @range=@{$info}; #[$from, $to, $type, $info]
		# skip this $substr range when its on the left side of $from, $to
		next if $range[1]<$from;
		# end the loop when $substr range is on the right side of $from, $to
		last if $range[0]>$to;
		# skip when $substr range is overlapping with the start of $from, $to, will let get_frag.pl deal with this
		next if ($range[0]<$from and $range[1]>=$from);
		# skip when $substr range is overlapping with the end of $from, $to, will let get_frag.pl deal with this
		next if ($range[0]>=$to and $range[1]>$to);
		# skip when $substr range is covering the entire $from, $to, will let get_frag.pl deal with this
		next if ($range[0]<$from and $range[1]>$to);
		# when $substr range is equal to or nested within $from, $to:
		if ($range[0]>=$from and $range[1]<=$to){
			# discard this range if it's the same type with $from, $to (a fragment)
			next if $range[2] eq $type;
			# discard this range if it's 80% covering the $from, $to but with different $type (misclassification)
			next if ($range[1]-$range[0]+1)/($to-$from+1) >= 0.8;
			# retain this range if the minuent entry is small and has a different $type
			$diff{"$chr:$range[0]:$range[1]"} = "$chr\t$range[0]\t$range[1]\t$range[2]\t$range[3]";
#print "$type\t$range[2]\t$range[3]\n" unless defined $range[3];
print $diff{"$chr:$range[0]:$range[1]"}."=$chr\t$range[0]\t$range[1]\t$range[2]\t$range[3]\n" unless defined $range[3];
			}
		}
	}
	}


