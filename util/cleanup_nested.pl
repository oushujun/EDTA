#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
use Data::Dumper;
#Shujun Ou (shujun.ou.1@gmail.com) 03/26/2019
#Update: 07/26/2019
#Update: 10/26/2019
#Update: 11/04/2019
#Update: 12/02/2020 by Sergei Ryazansky

my $usage = "\n
Iteratively clean up nested TE insertions and remove redundancy.

Further info:
Each sequence will be used as query to search the entire file.
For a subject sequence containing >95% of the query sequence, the matching part in the subject will be removed.
After removal, subject sequences shorter than the threadshold will be discarded.
The number of rounds of iterations is automatically decided (usually less than 8). User can also define this.

Usage:
perl cleanup_nested.pl -in file.fasta [options]
-in	[file]	Input sequence file in FASTA format
-cov	[float]	Minimum coverage of the query sequence to be considered as nesting. Default: 0.95
-minlen	[int]	Minimum length of the clean sequence to retain. Default: 80 (bp)
-miniden	[int]	Minimum identity of the clean sequence to retain. Default: 80 (%)
-clean	[int]	Clean nested sequences (1) or not (0). Default: 1
-maxcount	[int]	Specify the maximum number of stat lines you want to obtain. Default: 0 (no limit)
-iter	[int]	Numbers of iteration to remove redundency. Default: automatic
-blastplus [path]	Path to the blastn and makeblastdb program.
-threads|-t	[int]	Threads to run this script. Default: 4
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
my $blastplus = ""; #the path to blastn
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
	$blastplus=$ARGV[$k+1] if /^-blastplus$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
}

# checks
die "\nERROR: Input sequence file is not exist!\n$usage" unless -s $IN;
die "\nERROR: The -iter parameter receives non-integer input!\n$usage" unless $user_iter =~ /^[0-9]+$/;
$blastplus = "" unless defined $blastplus;
$blastplus=`command -v blastn 2>/dev/null` if $blastplus eq '';
$blastplus=~s/blastn\n//;
die "ERROR: blastn is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";

open IN, "<$IN" or die $!;
open STAT, ">$IN.stat" or die $!;

my %seq :shared;
my %touched_seq :shared; # here we will save all subject sequences that were modified after removing nested sequence
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

# itreatively remove redundant sequences and nested insertions
my $queue;
my $count_stat :shared = 0; # count stat lines in realtime
my $num_stat = 0; # count stat lines at the end of each iteration
$iter = $user_iter if $user_iter != 0;
for (my $i=0; $i<$iter; $i++){
	my $date=`date`;
	chomp ($date);
	print "$date\tClean up nested insertions and redundancy. Working on iteration $i\n";
	# write seq to a file and make blast db
	open Seq, ">$IN.iter$i" or die $!;
	%touched_seq = ();
	foreach my $id (sort {$a cmp $b} keys %seq){
		print Seq ">$id\n$seq{$id}\n";
		$touched_seq{$id}=0;
	}
	close Seq;
	`${blastplus}makeblastdb -in $IN.iter$i -dbtype nucl`;

	# multi-threading using queue, put candidate regions into queue for parallel computation
	$queue = Thread::Queue->new();
	foreach my $id (keys %seq){
		last unless defined $seq{$id};
		$queue -> enqueue([$id, $i, "$IN.iter$i"]);
	}
	$queue -> end(); # signal that no more items will be added to the queue

	# initiate a number of worker threads and run
	my @threads;
	foreach (1..$threads){
		push @threads, threads->create(\&condenser);
		#threads->create(\&condenser);
	}
	foreach my $thread (@threads){
		$thread->join();
	}
	`rm $IN.iter$i.nhr $IN.iter$i.nin $IN.iter$i.nsq $IN.iter$i.ndb $IN.iter$i.not $IN.iter$i.ntf $IN.iter$i.nto 2>/dev/null`;

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
		next if $touched_seq{$id} == 1;
		my $seq = ">$id\n$seq{$id}\n";
		my $length = length $seq{$id}; #query length
		next unless defined $length and $length > 0;
		my $exec="timeout 188s ${blastplus}blastn -query <(echo -e \"$seq\") -db $db -word_size 7 -evalue 1e-5 -dust no -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen\"";
		my @Blast=();
		my %seq_len; #store subject seq length info
		my %merged_hsps;
		my %merged_hsps_size; # collection of summing size of non-overlapped HSPs for checking the coverage
		my @aln_iden; # store alignment length and iden of this $id
		my $scaled_iden; # overall identity given all blast hits that pass the filter
		my $total_len; # total length of all alignments to $id
		@Blast=qx(bash -c '$exec' 2> /dev/null);
		# collect BLAST HSPs into the hash
		foreach (@Blast){
			my ($query, $subject, $iden, $len, $sbj_start, $sbj_end, $sbj_len) = (split)[0,1,2,3,8,9,11];
			my @vars = ($query, $subject, $iden, $len, $sbj_start, $sbj_end, $sbj_len);
			my @undefined_vars = grep { !defined($_) } @vars;
			next if @undefined_vars;
			next unless exists $seq{$subject};
			next if $query eq $subject;
			next if $touched_seq{$subject} == 1; # skip the iteration if the subject sequence was already modified (including discarded)
			next if $iden < $min_iden;
			next if $len < $minlen; # the length of HSPs should be more than 80 bp
			($sbj_start, $sbj_end) = ($sbj_end, $sbj_start) if $sbj_start > $sbj_end;
			push @{$merged_hsps{$subject}}, [$sbj_start, $sbj_end];
			$seq_len{$subject} = $sbj_len;
			$total_len += $len;
			push @aln_iden, [$len, $iden];
			}

		# calculate weighted identity
		foreach (@aln_iden) {
			my ($len, $iden) = @{$_};
			$scaled_iden += sprintf("%.3f", $iden * $len / $total_len);
			}

		# merge all overlapped HSPs and calculating the total covering by HSPs of subjects on the query
		my $merged = 0; # number of overlapping HSPs
		map {
			my $sbj = $_;
			# merging
			my ($ref1, $ref2) = &merger(@{$merged_hsps{$sbj}});
			@{$merged_hsps{$sbj}} = @$ref1;
			$merged = $$ref2;
			# total coverage by HSPs
			map {
				my ($start,$end) = ($_->[0],$_->[1]);
				$merged_hsps_size{$sbj} += $end - $start + 1;
			} @{$merged_hsps{$sbj}};
		} keys %merged_hsps;

		# removing the regions from the subject that are inserted into the query
		map {
			my ($sbj, $sbj_len) = ($_, $seq_len{$_});
			my $seq_new = $seq{$sbj};
			next unless defined $seq_new;
			next if length $seq_new ne $sbj_len; #if the subject length changes, it has been modified. Skip to avoid mismodification.
			my $poss = ''; # positions of the non-overlapped rHSPs egions that will be removed from the subject
			my ($qcov, $scov) = ($merged_hsps_size{$sbj}/$length, $merged_hsps_size{$sbj}/$sbj_len);
			$qcov = sprintf("%.3f", $qcov);

			if ($qcov >= $coverage or $scov >= $coverage) {
				# replace bases of HSPs regions to R (aka Remove); this masking is nessary since the subject sequence 
				# may be cleaned several times, for each non-overlapping merged HSPs regions.
				for my $hsp (@{$merged_hsps{$sbj}}) {
					my ($start, $end) = ($hsp->[0], $hsp->[1]);
					$poss = $poss . $start . ".." . $end . ",";
					my $len = $end - $start + 1;
					substr($seq_new, $start-1, $len) = "R" x $len if length $seq_new >= $start + $len - 1;
				}
				$seq_new =~ s/R//g;
				my $sbj_len_new = length $seq_new;
				if ($sbj_len_new >= $minlen and $sbj_len_new < length $seq{$sbj} and $clean == 1){
					print STAT "$sbj\tIter$i\tCleaned. $poss covering $qcov of $id; identity: $scaled_iden; merged $merged\n";
					$seq{$sbj} = $seq_new; #overwrite the sbj sequence if the new one is shorter
					$touched_seq{$sbj} = 1; # this subject sequence was modifed, and we will not deal with it any more in the current iteration
					$count_stat++;
				} elsif ($sbj_len_new < $minlen) {
					print STAT "$sbj\tIter$i\tDiscarded. Has only $sbj_len_new bp after cleaning by $id; identity: $scaled_iden; merged $merged\n";
					delete $seq{$sbj}; #delete this sequence if new seq is too short
					$touched_seq{$sbj} = 1; # this subject sequence was modifed (removed), and we will not deal with it any more in the current iteration
					$count_stat++;
				}
				# When $count_stat reaches the user defined stat count, we end the entire queue.
				if ($count_stat >= $count_limit and $count_limit > 0){
					print STAT "Reached user defined $count_limit of processed sequences at iter$i, stopping...\n\n";
					$queue->end();
					last;
				}
			}
		} keys %merged_hsps_size;
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
