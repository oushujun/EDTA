#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
#Shujun Ou (shujun.ou.1@gmail.com) 03/26/2019
#Update: 07/26/2019
#Update: 10/26/2019
#Update: 11/04/2019

my $usage = "\n
Iteratively clean up nested TE insertions and remove redundancy.

Further info:
	Each sequence will be used as query to search the entire file.
	For a subject sequence containing >95% of the query sequence, the matching part in the subject will be removed.
	After removal, subject sequences shorter than the threadshold will be diacarded.
	The number of rounds of iterations is automatically decided (usually less than 8). User can also define this.

Usage:
perl cleanup_nested.pl -in file.fasta [options]
	-in	[file]	Input sequence file in FASTA format
	-cov	[float]	Minimum coverage of the query sequence to be considered as nesting. Default: 0.95
	-minlen	[int]	Minimum length of the clean sequence to retain. Default: 80 (bp)
	-miniden	[int]	Minimum identity of the clean sequence to retain. Default: 80 (%)
	-iter	[int]	Numbers of iteration to remove redundency. Default: automatic
	-blastplus [path]	Path to the blastn and makeblastdb program.
	-threads|-t	[int]	Threads to run this script. Default: 4
\n";

my $IN = "";
my $coverage = 0.95; #if a subject sequence covers >95% of a query sequence, the matching part in the subject sequence will be removed.
my $minlen = 80; #minimal length >=80bp, otherwise discard the sequence
my $min_iden = 80; #minimal identity >=80%, otherwise discard the sequence
my $iter = 1;
my $user_iter = 0;
my $blastplus = ""; #the path to blastn
my $threads = 4;

my $k=0;
foreach (@ARGV){
	$IN=$ARGV[$k+1] if /^-in$/i and $ARGV[$k+1] !~ /^-/;
	$coverage=$ARGV[$k+1] if /^-cov$/i and $ARGV[$k+1] !~ /^-/;
	$minlen=$ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$min_iden=$ARGV[$k+1] if /^-miniden$/i and $ARGV[$k+1] !~ /^-/;
	$user_iter=$ARGV[$k+1] if /^-iter$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus=$ARGV[$k+1] if /^-blastplus$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
	}

# checks
die "\nERROR: Input sequence file is not exist!\n$usage" unless -s $IN;
die "\nERROR: The -iter parameter receives non-integer input!\n$usage" unless $iter =~ /^[0-9]+$/;
$blastplus = "" unless defined $blastplus;
$blastplus=`which blastn 2>/dev/null` if $blastplus eq '';
$blastplus=~s/blastn\n//;
die "ERROR: blastn is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";

open IN, "<$IN" or die $!;
open STAT, ">$IN.stat" or die $!;

my %seq :shared;
$/ = "\n>";
while (<IN>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id =~ s/\s+.*//;
	$seq =~ s/\s+//g;
	$seq{$id} = $seq;
	}
$/ = "\n";
close IN;

# itreatively remove redundant sequences and nested insertions
my $queue;
my $num_stat = 0;
$iter = $user_iter if $user_iter != 0;
for (my $i=0; $i<$iter; $i++){
	my $date=`date`;
	chomp ($date);
	print "$date\tClean up nested insertions and redundancy. Working on iteration $i\n";
	# write seq to a file and make blast db
	open Seq, ">$IN.iter$i" or die $!;
	foreach my $id (sort {$a cmp $b} keys %seq){
		print Seq ">$id\n$seq{$id}\n";
		}
	close Seq;
	`${blastplus}makeblastdb -in $IN.iter$i -dbtype nucl`;

	# multi-threading using queue, put candidate regions into queue for parallel computation
	$queue = Thread::Queue->new();
	foreach my $id (keys %seq){
		last unless defined $seq{$id};
		$queue -> enqueue([$id, $i, "$IN.iter$i"]);
		}
	$queue -> end();

	# initiate a number of worker threads and run
	foreach (1..$threads){
		threads->create(\&condenser);
		}
	foreach (threads -> list()){
		$_->join();
		}
	`rm $IN.iter$i.nhr $IN.iter$i.nin $IN.iter$i.nsq`;

	# automatically increase iteration based on the stat result
	my $curr_stat = `wc -l "$IN.stat"`;
	$curr_stat = (split /\s+/, $curr_stat)[0];
	if ($num_stat == $curr_stat){
		print "Saturated at iter$i, automatically stop.\n\n";
		last;
		} else {
		$num_stat = $curr_stat;
		$iter++ if $user_iter == 0;
		}
	}

# output clean sequence
open CLN, ">$IN.cln" or die $!;
foreach my $id (sort {$a cmp $b} keys %seq){
	print CLN ">$id\n$seq{$id}\n";
	}
close CLN;
close STAT;


# subrotine for the condenser
sub condenser(){
        while (defined($_ = $queue->dequeue())){
                my ($id, $i, $db) = (@{$_}[0], @{$_}[1], @{$_}[2]);
                next unless exists $seq{$id};
		my $seq = ">$id\n$seq{$id}\n";
		my $length = length $seq{$id};
		my $exec="timeout 188s ${blastplus}blastn -query <(echo -e \"$seq\") -db $db -outfmt 6 -word_size 7 -evalue 1e-5 -dust no";
		my @Blast=();
		@Blast=qx(bash -c '$exec' 2> /dev/null);
		foreach (@Blast){
			my ($query, $subject, $iden, $len, $sbj_start, $sbj_end) = (split)[0,1,2,3,8,9];
			next unless exists $seq{$subject};
			next if $query eq $subject and $iden == 100;
			next unless defined $length and $length > 0;
			next if $iden < $min_iden;
			my $cov = $len/$length;
			next unless $cov >= $coverage;
			($sbj_start, $sbj_end) = ($sbj_end, $sbj_start) if $sbj_start > $sbj_end;
			my ($sbj_seq, $sbj_seq_p1, $sbj_seq_p2, $sbj_seq_new) = ($seq{$subject}, '', '', '');
			next unless defined $sbj_seq;
			my $sbj_len = length $sbj_seq;
			$sbj_seq_p1 = substr $sbj_seq, 1, $sbj_start if $sbj_start > 1;
			$sbj_seq_p2 = substr $sbj_seq, $sbj_end if $sbj_end < $sbj_len;
			$sbj_seq_new = "$sbj_seq_p1"."$sbj_seq_p2";
			my $sbj_len_new = length $sbj_seq_new;
			if ($sbj_len_new >= $minlen){
				print STAT "$subject\tIter$i\tCleaned. $sbj_start..$sbj_end covering $cov of $query for identity $iden%\n";
				$seq{$subject} = $sbj_seq_new; #update sequence, overwrite the current sequence
				} else {
				print STAT "$subject\tIter$i\tDiscarded. Has only $sbj_len_new bp after cleaning by $query\n";
				delete $seq{$subject}; #delete this sequence if new seq is too short
				}
			}
		}
	}


