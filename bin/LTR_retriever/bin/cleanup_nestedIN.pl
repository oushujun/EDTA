#!/usr/bin/env perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com) 03/26/2019

my $usage = "\n
Clean up LTR element nested insertions.

Further info:
	Each sequence will be used as query to search the entire file.
	For a subject sequence containing >95% of the query sequence, the matching part in the subject will be removed.
	After removal, subject sequences shorter than the threadshold will be diacarded.

Usage: cleanup_nestedIN.pl -in file.fasta > file.fasta.cln
Options:	-in		Input LTR sequence file in FASTA format
		-cov		Minimum coverage of the query sequence to be considered as nesting. Default: 0.95
		-minlen		Minimum length of the clean sequence to retain. Default: 100 (bp)
		-blastplus	Path to the blastn and makeblastdb program.
		-threads	Threads to run blastn

Free advise: Iterate this step 2-3 times for more thorough removal of nested insertions.
\n";

my $IN = "";
my $coverage = 0.95; #if a subject sequence covers >95% of a query sequence, the matching part in the subject sequence will be removed.
my $minlen = 100; #minimal length >=100bp, otherwise discard the sequence
my $blastplus = ""; #the path to blastn
my $threads = 4;

my $k=0;
foreach (@ARGV){
	$IN=$ARGV[$k+1] if /^-in$/i;
	$coverage=$ARGV[$k+1] if /^-cov$/i;
	$minlen=$ARGV[$k+1] if /^-minlen$/i;
	$blastplus=$ARGV[$k+1] if /^-blastplus$/i;
	$threads=$ARGV[$k+1] if /^-threads$/i;
	$k++;
	}

`${blastplus}makeblastdb -in $IN -dbtype nucl`;
open IN, "<$IN" or die "\n\nInput sequence file is not exist!\n\n$usage";

my %seq;
$/ = "\n>";
while (<IN>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	$seq{$id} = $seq;
	}
$/ = "\n";
close IN;

foreach my $id (keys %seq){
	next unless exists $seq{$id};
	my $seq = ">$id\n$seq{$id}\n";
	my $length = length $seq{$id};
	my $exec="${blastplus}blastn -outfmt 6 -num_threads $threads -query <(echo -e \"$seq\") -db $IN";
	my @Blast=();
	@Blast=qx(bash -c '$exec' 2> /dev/null);
	foreach (@Blast){
		my ($query, $subject, $iden, $len, $sbj_start, $sbj_end) = (split)[0,1,2,3,8,9];
		next if $query eq $subject and $iden == 100;
		next unless $len/$length >= $coverage;
		next unless exists $seq{$subject};
		($sbj_start, $sbj_end) = ($sbj_end, $sbj_start) if $sbj_start > $sbj_end;
		my ($sbj_seq, $sbj_seq_p1, $sbj_seq_p2, $sbj_seq_new) = ($seq{$subject}, '', '', '');
		my $sbj_len = length $sbj_seq;
		$sbj_seq_p1 = substr $sbj_seq, 1, $sbj_start if $sbj_start > 1;
		$sbj_seq_p2 = substr $sbj_seq, $sbj_end if $sbj_end < $sbj_len;
		$sbj_seq_new = "$sbj_seq_p1"."$sbj_seq_p2";
		my $sbj_len_new = length $sbj_seq_new;
		if ($sbj_len_new >= $minlen){
			$seq{$subject} = $sbj_seq_new; #update sequence
			} else {
			delete $seq{$subject}; #delete this sequence if new seq is too short
			}
		}
	}

foreach my $id (sort{$a cmp $b} (keys %seq)){
	print ">$id\n$seq{$id}\n";
	}
       
