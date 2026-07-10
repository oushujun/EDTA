#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
use IO::Handle;
use Data::Dumper;
#Shujun Ou (shujun.ou.1@gmail.com) 03/26/2019
#Update: 07/26/2019
#Update: 10/26/2019
#Update: 11/04/2019
#Update: 12/02/2020 by Sergei Ryazansky
#Update: 07/07/2026 query-level checkpoint/resume: auto-continues a killed run from its
#	partial .stat + per-iteration snapshots, skipping already-processed queries. Use
#	-overwrite 1 to force a clean restart.

my $usage = "\n
Iteratively clean up nested TE insertions and remove redundancy.

Further info:
Each sequence will be used as query to search the entire file.
For a subject sequence containing >95% of the query sequence, the matching part in the subject will be removed.
After removal, subject sequences shorter than the threadshold will be discarded.
The number of rounds of iterations is automatically decided (usually less than 8). User can also define this.

This script is check-pointed: if a previous run was interrupted, rerunning on the same input
auto-resumes from the partial .stat file and the per-iteration snapshots, re-BLASTing only the
queries that had not been processed yet. Pass -overwrite 1 to ignore any partial results and
restart cleanly.

Usage:
perl cleanup_nested.pl -in file.fasta [options]
-in	[file]	Input sequence file in FASTA format
-cov	[float]	Minimum coverage of the query sequence to be considered as nesting. Default: 0.95
-minlen	[int]	Minimum length of the clean sequence to retain. Default: 80 (bp)
-miniden	[int]	Minimum identity of the clean sequence to retain. Default: 80 (%)
-clean	[int]	Clean nested sequences (1) or not (0). Default: 1
-maxcount	[int]	Specify the maximum number of stat lines you want to obtain. Default: 0 (no limit)
-iter	[int]	Numbers of iteration to remove redundency. Default: automatic
-overwrite	[0|1]	Ignore any partial results and restart cleanly (1), or auto-resume if found (0, default).
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
my $overwrite = 0; # 1, ignore partial results and restart; 0, auto-resume from partial results if found.
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
	$overwrite=$ARGV[$k+1] if /^-overwrite$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus=$ARGV[$k+1] if /^-blastplus$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
}

# checks
die "\nERROR: Input sequence file is not exist!\n$usage" unless -s $IN;
die "\nERROR: The -iter parameter receives non-integer input!\n$usage" unless $user_iter =~ /^[0-9]+$/;
die "\nERROR: The -overwrite parameter expects 0 or 1!\n$usage" unless $overwrite =~ /^[01]$/;
$blastplus = "" unless defined $blastplus;
$blastplus=`command -v blastn 2>/dev/null` if $blastplus eq '';
$blastplus=~s/blastn\n//i;
die "ERROR: blastn is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";

# shared state
my %seq :shared;
my %touched_seq :shared; # here we will save all subject sequences that were modified after removing nested sequence
my $count_stat :shared = 0; # count stat lines in realtime (drives -maxcount and the saturation check)
my $stat_lock :shared;      # serializes writes to STAT so lines never interleave (also keeps .stat replayable)
my $queue;
my $num_stat = 0; # stat count at the end of the previous iteration; drives the saturation check

# Decide whether to resume from a previous interrupted run.
# A resumable checkpoint = a partial .stat plus the snapshot of the iteration that was in
# progress ($IN.iter<L>). The snapshots are written atomically (tmp + rename), so the highest
# $IN.iter<N> present is guaranteed complete and reflects the start-of-iteration-N state.
my $resume = 0;
my $start_iter = 0;
my %done_q; # query ids already processed in the resumed iteration (skip re-BLASTing them)
my $L = ($overwrite == 1) ? -1 : &latest_iter($IN);
if ($L >= 0 and -s "$IN.stat" and -s "$IN.iter$L"){
	$resume = 1;
	$start_iter = $L;
} else {
	# fresh run (no valid checkpoint, or -overwrite 1): clear any stale checkpoint artifacts
	unlink glob "$IN.iter*";
}

# load sequences and open the stat file
if ($resume){
	&load_fasta("$IN.iter$start_iter");      # start-of-iteration-L state
	&replay_stat("$IN.stat", $start_iter);   # re-apply this iteration's Cleaned/Discarded events -> mid-iteration state
	my ($total, $by_iter) = &count_stat_lines("$IN.stat");
	$count_stat = $total;                    # total events so far (for -maxcount)
	$num_stat = 0;
	$num_stat += $by_iter->{$_} foreach grep { $_ < $start_iter } keys %$by_iter; # events through iter L-1
	&load_done(\%done_q, $start_iter);
	open STAT, ">>$IN.stat" or die $!;       # append; keep the partial results
	my $date=`date`; chomp $date;
	print "$date\tResuming cleanup_nested at iteration $start_iter: ", scalar(keys %seq),
		" sequences restored, $count_stat stat lines kept, ", scalar(keys %done_q),
		" queries already processed (will be skipped).\n";
} else {
	&load_fasta($IN);
	open STAT, ">$IN.stat" or die $!;
}
STAT->autoflush(1); # flush each line promptly so a kill truncates at most the final line

# ensure the loop actually runs the resumed iteration under automatic iteration counting
$iter = $user_iter if $user_iter != 0;
$iter = $start_iter + 1 if $resume and $user_iter == 0 and $iter <= $start_iter;

# iteratively remove redundant sequences and nested insertions
for (my $i=$start_iter; $i<$iter; $i++){
	my $date=`date`;
	chomp ($date);
	print "$date\tClean up nested insertions and redundancy. Working on iteration $i\n";

	# write seq to a file and make blast db.
	# On the resumed iteration we reuse the existing snapshot and the replayed %seq/%touched_seq state.
	unless ($resume and $i == $start_iter){
		open Seq, ">$IN.iter$i.tmp" or die $!;
		%touched_seq = ();
		foreach my $id (sort {$a cmp $b} keys %seq){
			print Seq ">$id\n$seq{$id}\n";
			$touched_seq{$id}=0;
		}
		close Seq;
		rename "$IN.iter$i.tmp", "$IN.iter$i" or die $!; # atomic: snapshot appears only once complete
	}
	`${blastplus}makeblastdb -in $IN.iter$i -dbtype nucl`;

	# multi-threading using queue, put candidate regions into queue for parallel computation
	$queue = Thread::Queue->new();
	foreach my $id (keys %seq){
		last unless defined $seq{$id};
		next if $resume and $i == $start_iter and exists $done_q{$id}; # already processed before the interruption
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
	`rm $IN.iter$i.nhr $IN.iter$i.nin $IN.iter$i.nsq $IN.iter$i.ndb $IN.iter$i.not $IN.iter$i.ntf $IN.iter$i.nto $IN.iter$i.njs 2>/dev/null`;

	# stop cleanly once the -maxcount cap is reached (no need to spin through remaining iterations)
	last if $count_limit > 0 and $count_stat >= $count_limit;

	# automatically increase iteration based on the stat result (saturation = an iteration adds no new events)
	my $curr_stat = $count_stat;
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

# successful completion: remove checkpoint artifacts (snapshots, done-logs, tmp) so a future
# run on this input starts fresh rather than trying to resume a finished job.
unlink glob "$IN.iter*";


# subrotine for the condenser
sub condenser(){
	my $done_fh; # per-thread log of processed queries (opened lazily; no cross-thread contention)
	while (defined($_ = $queue->dequeue())){
		my ($id, $i, $db) = (@{$_}[0], @{$_}[1], @{$_}[2]);
		unless (defined $done_fh){
			open($done_fh, ">>", "$IN.iter$i.done." . threads->tid()) or warn "cannot open done-log: $!\n";
			$done_fh->autoflush(1) if defined $done_fh;
		}
		next unless exists $seq{$id};
		next if $touched_seq{$id} == 1;
		my $seq = ">$id\n$seq{$id}\n";
		my $length = length $seq{$id}; #query length
		next unless defined $length and $length > 0;
		my $exec="timeout -s KILL 120s ${blastplus}blastn -query <(echo -e \"$seq\") -db $db -word_size 7 -evalue 1e-5 -dust no -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen\"";
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
			$scov = sprintf("%.3f", $scov);

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
					{ lock($stat_lock); print STAT "$sbj\tIter$i\tCleaned. $poss covering $qcov of $id; scov: $scov; identity: $scaled_iden; merged $merged\n"; $count_stat++; }
					$seq{$sbj} = $seq_new; #overwrite the sbj sequence if the new one is shorter
					$touched_seq{$sbj} = 1; # this subject sequence was modifed, and we will not deal with it any more in the current iteration
				} elsif ($sbj_len_new < $minlen) {
					{ lock($stat_lock); print STAT "$sbj\tIter$i\tDiscarded. Has only $sbj_len_new bp after cleaning by $id; qcov: $qcov; scov: $scov; identity: $scaled_iden; merged $merged\n"; $count_stat++; }
					delete $seq{$sbj}; #delete this sequence if new seq is too short
					$touched_seq{$sbj} = 1; # this subject sequence was modifed (removed), and we will not deal with it any more in the current iteration
				}
				# When $count_stat reaches the user defined stat count, we end the entire queue.
				if ($count_stat >= $count_limit and $count_limit > 0){
					{ lock($stat_lock); print STAT "Reached user defined $count_limit of processed sequences at iter$i, stopping...\n\n"; }
					$queue->end();
					last;
				}
			}
		} keys %merged_hsps_size;

		# this query is fully processed; record it so a resumed run will not re-BLAST it
		print $done_fh "$id\n" if defined $done_fh;
	}
	close $done_fh if defined $done_fh;
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

# ---- checkpoint / resume helpers ----

# load a FASTA into the shared %seq / %touched_seq (touched reset to 0)
sub load_fasta {
	my ($file) = @_;
	open(my $fh, "<", $file) or die "ERROR: cannot read $file: $!\n";
	local $/ = "\n>";
	while (<$fh>){
		s/>//g;
		my ($id, $s) = (split /\n/, $_, 2);
		next unless defined $id and defined $s;
		$id =~ s/\s+.*//;
		$s =~ s/\s+//g;
		next if $id eq '';
		$seq{$id} = $s;
		$touched_seq{$id} = 0;
	}
	close $fh;
}

# highest N for which the snapshot $IN.iter<N> exists (atomic writes guarantee it is complete); -1 if none
sub latest_iter {
	my ($in) = @_;
	my $max = -1;
	foreach my $f (glob "$in.iter*"){
		if ($f =~ /\Q$in\E\.iter(\d+)$/){
			$max = $1 if $1 > $max;
		}
	}
	return $max;
}

# re-apply a single recorded "Cleaned." event to reconstruct the subject's current sequence
sub apply_clean {
	my ($sbj, $poss) = @_;
	return unless exists $seq{$sbj};
	my $seq_new = $seq{$sbj};
	foreach my $range (split /,/, $poss){
		next unless $range =~ /^(\d+)\.\.(\d+)$/;
		my ($start, $end) = ($1, $2);
		my $len = $end - $start + 1;
		substr($seq_new, $start-1, $len) = "R" x $len if length $seq_new >= $start + $len - 1;
	}
	$seq_new =~ s/R//g;
	$seq{$sbj} = $seq_new;
	$touched_seq{$sbj} = 1;
}

# replay the Cleaned/Discarded events of iteration $L from the partial .stat to restore
# %seq/%touched_seq to the mid-iteration state (each subject is modified at most once per iteration)
sub replay_stat {
	my ($file, $L) = @_;
	open(my $fh, "<", $file) or return;
	while (<$fh>){
		if (/^(\S+)\tIter(\d+)\tCleaned\. ([\d.,]+) covering/){
			&apply_clean($1, $3) if $2 == $L;
		} elsif (/^(\S+)\tIter(\d+)\tDiscarded\./){
			if ($2 == $L){ delete $seq{$1}; $touched_seq{$1} = 1; }
		}
	}
	close $fh;
}

# count Cleaned/Discarded events in .stat: returns (total, {iter => count})
sub count_stat_lines {
	my ($file) = @_;
	my %by_iter;
	my $total = 0;
	open(my $fh, "<", $file) or return (0, {});
	while (<$fh>){
		if (/^\S+\tIter(\d+)\t(?:Cleaned|Discarded)\./){
			$by_iter{$1}++;
			$total++;
		}
	}
	close $fh;
	return ($total, \%by_iter);
}

# union the per-thread processed-query logs for iteration $i into %$set
sub load_done {
	my ($set, $i) = @_;
	foreach my $df (glob "$IN.iter$i.done.*"){
		open(my $fh, "<", $df) or next;
		while (my $q = <$fh>){
			chomp $q;
			$set->{$q} = 1 if length $q;
		}
		close $fh;
	}
}
