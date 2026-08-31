#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com) 03/26/2019
#Update: 07/26/2019
#Update: 10/26/2019
#Update: 11/04/2019
#Update: 12/02/2020 by Sergei Ryazansky
#Update: 07/10/2025 minimap2 batch alignment replacing per-query blastn

my $usage = "\n
Iteratively clean up nested TE insertions and remove redundancy.
Uses minimap2 for batch all-vs-all alignment instead of per-query blastn.

Further info:
Each sequence will be used as query to search the entire file.
For a subject sequence containing >95% of the query sequence, the matching part in the subject will be removed.
After removal, subject sequences shorter than the threadshold will be discarded.
The number of rounds of iterations is automatically decided (usually less than 8). User can also define this.

Usage:
perl cleanup_nested_minimap2.pl -in file.fasta [options]
-in	[file]	Input sequence file in FASTA format
-cov	[float]	Minimum coverage of the query sequence to be considered as nesting. Default: 0.95
-minlen	[int]	Minimum length of the clean sequence to retain. Default: 80 (bp)
-miniden	[int]	Minimum identity of the clean sequence to retain. Default: 80 (%)
-clean	[int]	Clean nested sequences (1) or not (0). Default: 1
-maxcount	[int]	Specify the maximum number of stat lines you want to obtain. Default: 0 (no limit)
-iter	[int]	Numbers of iteration to remove redundency. Default: automatic
-minimap2 [path]	Path to the minimap2 program.
-kmer	[int]	Minimizer k-mer size for minimap2 seeding. Default: 7
		  Approximate blastn -word_size equivalence:
		    -kmer 7  ~ blastn -word_size 7  (sensitive, short/divergent seqs)
		    -kmer 10 ~ blastn -word_size 11 (faster, longer seqs)
		    -kmer 15 ~ blastn -word_size 16 (fast, near-identical seqs)
-mwindow [int]	Minimizer window size for minimap2. Default: 5
		  Smaller = more sensitive but slower. Recommended: ~2/3 of -kmer.
		    -mwindow 5  (sensitive, pairs with -kmer 7)
		    -mwindow 7  (moderate, pairs with -kmer 10)
		    -mwindow 10 (fast, pairs with -kmer 15)
-threads|-t	[int]	Threads to run minimap2. Default: 4
\n";

my $IN = "";
my $coverage = 0.95; #if a subject sequence covers >95% of a query sequence, the matching part in the subject sequence will be removed.
my $minlen = 80; #minimal length >=80bp, otherwise discard the sequence
my $min_iden = 80; #minimal identity >=80%, otherwise discard the sequence
my $offset = 7; #if two blast hits are less than $offset [default=7bp) away from each other, join them as one hit
my $clean = 1; #1, clean nested sequences; 0, will not clean nested, only discard highly overlapping (~100%) sequences
my $iter = 1;
my $user_iter = 0;
my $count_limit = 0; # the maximum number of stat lines you want to obtain. 0 = no limit.
my $minimap2 = ""; #the path to minimap2
my $mm_kmer = 7;   # minimap2 -k: minimizer k-mer size (~blastn -word_size)
my $mm_window = 5; # minimap2 -w: minimizer window size (smaller = more sensitive)
my $threads = 4;

my $k=0;
foreach (@ARGV){
	$IN=$ARGV[$k+1] if /^-in$/i and $ARGV[$k+1] !~ /^-/;
	$coverage=$ARGV[$k+1] if /^-cov$/i and $ARGV[$k+1] !~ /^-/;
	$minlen=$ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$min_iden=$ARGV[$k+1] if /^-miniden$/i and $ARGV[$k+1] !~ /^-/;
	$clean=$ARGV[$k+1] if /^-clean$/i and $ARGV[$k+1] !~ /^-/;
	$count_limit=$ARGV[$k+1] if /^-maxcount$/i and $ARGV[$k+1] !~ /^-/;
	$user_iter=$ARGV[$k+1] if /^-iter$/i and $ARGV[$k+1] !~ /^-/;
	$minimap2=$ARGV[$k+1] if /^-minimap2$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
	$mm_kmer=$ARGV[$k+1] if /^-kmer$/i and $ARGV[$k+1] !~ /^-/;
	$mm_window=$ARGV[$k+1] if /^-mwindow$/i and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
}

# checks
die "\nERROR: Input sequence file is not exist!\n$usage" unless -s $IN;
die "\nERROR: The -iter parameter receives non-integer input!\n$usage" unless $user_iter =~ /^[0-9]+$/;
$minimap2 = "" unless defined $minimap2;
$minimap2=`command -v minimap2 2>/dev/null` if $minimap2 eq '';
chomp $minimap2;
die "ERROR: minimap2 is not found in PATH or the provided path ($minimap2)!\n" unless -X $minimap2;

open IN, "<$IN" or die $!;
open STAT, ">$IN.stat" or die $!;

my %seq;
my %touched_seq;
$/ = "\n>";
while (<IN>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id =~ s/\s+.*//;
	$seq =~ s/\s+//g;
	$seq{$id} = $seq;
	$touched_seq{$id} = 0;
}
$/ = "\n";
close IN;

# iteratively remove redundant sequences and nested insertions
my $count_stat = 0; # count stat lines in realtime
my $num_stat = 0; # count stat lines at the end of each iteration
my $early_stop = 0;
$iter = $user_iter if $user_iter != 0;
for (my $i=0; $i<$iter; $i++){
	my $date=`date`;
	chomp ($date);
	print "$date\tClean up nested insertions and redundancy. Working on iteration $i\n";

	# write seq to a file for minimap2
	open Seq, ">$IN.iter$i" or die $!;
	%touched_seq = ();
	foreach my $id (sort {$a cmp $b} keys %seq){
		print Seq ">$id\n$seq{$id}\n";
		$touched_seq{$id}=0;
	}
	close Seq;

	# Phase 1: run minimap2 all-vs-all (single call, multi-threaded internally)
	my $paf_file = "$IN.iter$i.paf";
	system("$minimap2 -c -t $threads -k $mm_kmer -w $mm_window -p 0.01 -N 1000 $IN.iter$i $IN.iter$i > $paf_file 2>/dev/null");

	# Phase 2: stream PAF and process per-query batches
	open PAF, "<$paf_file" or die "Cannot open $paf_file: $!\n";
	my $prev_query = "";
	my %query_hsps;       # {subject} => [[start, end], ...]
	my %query_seq_len;    # {subject} => tlen
	my @query_aln_iden;   # [[alen, iden], ...]
	my $query_total_len = 0;
	my $query_length = 0; # length of current query

	while (<PAF>){
		chomp;
		my @f = split /\t/;
		# PAF: 0=qname 1=qlen 2=qstart 3=qend 4=strand 5=tname 6=tlen 7=tstart 8=tend 9=nmatch 10=alen 11=mapq
		my ($qname, $qlen, $tname, $tlen, $tstart, $tend, $nmatch, $alen) =
		   ($f[0], $f[1], $f[5], $f[6], $f[7], $f[8], $f[9], $f[10]);

		# detect query boundary -> process previous query's batch
		if ($qname ne $prev_query && $prev_query ne ""){
			process_query_batch($prev_query, $query_length, $i,
				\%query_hsps, \%query_seq_len, \@query_aln_iden, $query_total_len);
			last if $early_stop;
			%query_hsps = ();
			%query_seq_len = ();
			@query_aln_iden = ();
			$query_total_len = 0;
		}
		$prev_query = $qname;
		$query_length = $qlen;

		# filter
		next unless defined $tname && defined $tlen && defined $nmatch && defined $alen;
		next unless $alen > 0;
		next unless exists $seq{$tname};
		next if $qname eq $tname;              # skip self-hit
		next unless exists $touched_seq{$tname};
		next if $touched_seq{$tname} == 1;     # skip already-modified subjects

		# compute gap-excluded identity from PAF CIGAR
		# minimap2's nmatch/alen includes gap columns in denominator, producing
		# systematically lower identity than blastn for gappy alignments.
		# Parse CIGAR to count gap bases (I/D ops) and exclude them.
		my $gap_bases = 0;
		if (/cg:Z:(\S+)/) {
			my $cigar = $1;
			while ($cigar =~ /(\d+)[ID]/g) {
				$gap_bases += $1;
			}
		}
		my $non_gap_cols = $alen - $gap_bases; # aligned bases only (matches + mismatches)
		my $iden = ($non_gap_cols > 0) ? ($nmatch / $non_gap_cols) * 100 : 0;
		next if $iden < $min_iden;

		# convert PAF 0-based half-open to 1-based inclusive
		my $sbj_start = $tstart + 1;
		my $sbj_end   = $tend;  # half-open end = 1-based inclusive end
		my $aln_len = $sbj_end - $sbj_start + 1;
		next if $aln_len < $minlen;

		push @{$query_hsps{$tname}}, [$sbj_start, $sbj_end];
		$query_seq_len{$tname} = $tlen;
		$query_total_len += $alen;
		push @query_aln_iden, [$alen, $iden];
	}
	# process final query batch
	if ($prev_query ne "" && !$early_stop){
		process_query_batch($prev_query, $query_length, $i,
			\%query_hsps, \%query_seq_len, \@query_aln_iden, $query_total_len);
	}
	close PAF;

	# cleanup iteration files
	unlink "$IN.iter$i", $paf_file;

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
	last if $early_stop;
}

# output clean sequence
open CLN, ">$IN.cln" or die $!;
foreach my $id (sort {$a cmp $b} keys %seq){
	print CLN ">$id\n$seq{$id}\n";
}
close CLN;
close STAT;


# subroutine to process one query's batch of HSPs
sub process_query_batch {
	my ($id, $length, $i, $hsps_ref, $seq_len_ref, $aln_iden_ref, $total_len) = @_;

	return unless exists $seq{$id};
	return if $touched_seq{$id} == 1;
	return unless defined $length and $length > 0;
	return unless keys %$hsps_ref;

	# calculate weighted identity
	my $scaled_iden = 0;
	foreach my $pair (@$aln_iden_ref) {
		my ($len, $iden) = @$pair;
		$scaled_iden += sprintf("%.3f", $iden * $len / $total_len);
	}

	# merge all overlapped HSPs and calculate total coverage per subject
	my $merged = 0;
	my %merged_hsps_size;
	foreach my $sbj (keys %$hsps_ref) {
		my ($ref1, $ref2) = merger(@{$hsps_ref->{$sbj}});
		$hsps_ref->{$sbj} = $ref1;
		$merged = $$ref2;
		foreach my $interval (@$ref1) {
			my ($start, $end) = ($interval->[0], $interval->[1]);
			$merged_hsps_size{$sbj} += $end - $start + 1;
		}
	}

	# removing the regions from the subject that are inserted into the query
	foreach my $sbj (keys %merged_hsps_size) {
		my $sbj_len = $seq_len_ref->{$sbj};
		my $seq_new = $seq{$sbj};
		next unless defined $seq_new;
		next if length $seq_new ne $sbj_len; #if the subject length changes, it has been modified. Skip to avoid mismodification.
		my $poss = ''; # positions of the non-overlapped HSP regions that will be removed from the subject
		my ($qcov, $scov) = ($merged_hsps_size{$sbj}/$length, $merged_hsps_size{$sbj}/$sbj_len);
		$qcov = sprintf("%.3f", $qcov);
		$scov = sprintf("%.3f", $scov);

		if ($qcov >= $coverage or $scov >= $coverage) {
			# replace bases of HSP regions to R (aka Remove)
			for my $hsp (@{$hsps_ref->{$sbj}}) {
				my ($start, $end) = ($hsp->[0], $hsp->[1]);
				$poss = $poss . $start . ".." . $end . ",";
				my $len = $end - $start + 1;
				substr($seq_new, $start-1, $len) = "R" x $len if length $seq_new >= $start + $len - 1;
			}
			$seq_new =~ s/R//g;
			my $sbj_len_new = length $seq_new;
			if ($sbj_len_new >= $minlen and $sbj_len_new < length $seq{$sbj} and $clean == 1){
				print STAT "$sbj\tIter$i\tCleaned. $poss covering $qcov of $id; scov: $scov; identity: $scaled_iden; merged $merged\n";
				$seq{$sbj} = $seq_new;
				$touched_seq{$sbj} = 1;
				$count_stat++;
			} elsif ($sbj_len_new < $minlen) {
				print STAT "$sbj\tIter$i\tDiscarded. Has only $sbj_len_new bp after cleaning by $id; qcov: $qcov; scov: $scov; identity: $scaled_iden; merged $merged\n";
				delete $seq{$sbj};
				$touched_seq{$sbj} = 1;
				$count_stat++;
			}
			# When $count_stat reaches the user defined stat count, stop
			if ($count_stat >= $count_limit and $count_limit > 0){
				print STAT "Reached user defined $count_limit of processed sequences at iter$i, stopping...\n\n";
				$early_stop = 1;
				last;
			}
		}
	}
}

sub merger() {
	my @hsps = @_;
	my $merged = 0;
	my @intervals = sort {
		$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]
	} @hsps;
	my @merged;
	my $current = $intervals[0];
	for my $i (1..$#intervals) {
		if ($intervals[$i][0] > $current->[1] + $offset) { # allow 7bp offset
			push @merged, $current;
			$current = $intervals[$i];
		} else {
			next unless $intervals[$i][1] > $current->[1];
			$current->[1] = $intervals[$i][1];
			$merged++;
		}
	}
	push @merged, $current;
	return (\@merged, \$merged);
}
