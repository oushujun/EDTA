#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;

#function: For regions in the subtrahend.list that are nested within or equal to regions in the minuend.list, do;
#		1. Discard this range if it's the same type with $from, $to (a fragment)
#		2. Discard this range if it's 80% covering the $from, $to but with different $type (misclassification)
#		3. Retain this range if the minuent entry is small and has a different $type
#usage: perl keep_nest.pl minuend.list subtrahend.list thread_num
#Author: Shujun Ou (shujun.ou.1@gmail.com), 08/02/2019

my $usage = "\n\tperl keep_nest.pl minuend.list subtrahend.list thread_num\n\n";

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
	chomp;
#	my ($chr, $from, $to, $type, $info)=(split /\s+/, $_, 5);
	my ($chr, $from, $to, $type)=(split)[0,1,2,11];
	push @{$substr{$chr}}, [$from||0, $to||0, $type, $_];
	}

# presort subtrahend intervals by start and build a prefix max-end array for fast nested scans
my %maxend;
foreach my $chr (keys %substr){
	@{$substr{$chr}} = sort { $a->[0] <=> $b->[0] || $a->[3] cmp $b->[3] } @{$substr{$chr}};
	my $max;
	$maxend{$chr} = [ map { $max = $_->[1] if !defined $max || $_->[1] > $max; $max } @{$substr{$chr}} ];
	}

## multi-threading using queue, put candidate regions into queue for parallel computation
my %diff :shared;
my $queue = Thread::Queue -> new();
my $nest_id = 0;
while (<Minuend>){
	next if /^\s+$/;
	chomp;
#	my ($chr, $from, $to, $type, $info)=(split /\s+/, $_, 5);
	my ($chr, $from, $to, $type)=(split)[0,1,2,11];
	next unless defined $chr;
	$diff{"$chr:$from:$to:$type#".++$nest_id} = $_; #all minuend info are retained
	$queue->enqueue([$chr, $from, $to, $type]);
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
	my ($chr, $from, $to, $type) = (@{$_}[0], @{$_}[1], @{$_}[2], @{$_}[3]);
	my $list = $substr{$chr};
	next unless defined $list and @$list;
	# binary search for the first subtrahend interval starting after $to
	my ($lo, $hi) = (0, scalar @$list);
	while ($lo < $hi){
		my $mid = int(($lo + $hi) / 2);
		if ($list->[$mid][0] <= $to){ $lo = $mid + 1 } else { $hi = $mid }
		}
	for (my $j = $lo - 1; $j >= 0; $j--){
		# no subtrahend interval at or left of $j has an end reaching $from, so none of them is relevant
		last if $maxend{$chr}[$j] < $from;
		my @range=@{$list->[$j]}; #[$from, $to, $type, $_]
		# skip this $substr range when its on the left side of $from, $to
		next if $range[1]<$from;
		# skip when $substr range is overlapping with the start of $from, $to, will let get_frag.pl deal with this
		next if ($range[0]<$from and $range[1]>=$from);
		# skip when $substr range is overlapping with the end of $from, $to, will let get_frag.pl deal with this
		next if ($range[0]<=$to and $range[1]>$to);
		# skip when $substr range is covering the entire $from, $to, will let get_frag.pl deal with this
		next if ($range[0]<$from and $range[1]>$to);
		# when $substr range is equal to or nested within $from, $to:
		if ($range[0]>=$from and $range[1]<=$to){
			# discard this range if it's the same type with $from, $to (a fragment)
			next if $range[2] eq $type;
			# discard this range if it's 80% covering the $from, $to but with different $type (misclassification)
			next if ($to-$from+1) > 0 and ($range[1]-$range[0]+1)/($to-$from+1) >= 0.8;
			# retain this range if the minuent entry is small and has a different $type
			$diff{"$chr:$range[0]:$range[1]:$range[2]#s$j"} = $range[3];
			#$diff{"$chr:$range[0]:$range[1]"} = "$chr\t$range[0]\t$range[1]\t$range[2]\t$range[3]";
			#print "$type\t$range[2]\t$range[3]\n" unless defined $range[3];
			#print $diff{"$chr:$range[0]:$range[1]"}."=$chr\t$range[0]\t$range[1]\t$range[2]\t$range[3]\n" unless defined $range[3];
			}
		}
	}
	}


