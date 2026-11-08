#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
use IO::Handle;
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
#Shujun Ou (shujun.ou.1@gmail.com) 03/26/2019
#Update: 07/26/2019
#Update: 10/26/2019
#Update: 11/04/2019
#Update: 12/02/2020 by Sergei Ryazansky
#Update: 07/07/2026 query-level checkpoint/resume: auto-continues a killed run from its
#	partial .stat + per-iteration snapshots, skipping already-processed queries. Use
#	-overwrite 1 to force a clean restart.
#Update: 11/06/2026 the same output on every run: threads only run BLAST, and the main thread
#	applies the results one query at a time in a fixed order. A fragment of an entry of its own
#	class no longer cuts that entry; the entry discards it instead (-sameclass).

my $usage = "\n
Iteratively clean up nested TE insertions and remove redundancy.

Further info:
Each sequence will be used as query to search the entire file.
For a subject sequence containing >95% of the query sequence, the matching part in the subject will be removed.
If the two IDs carry the same class (name#class, e.g. #LINE/L1; #Unknown does not count), the query is a
fragment of the subject rather than an insertion in it: the subject is kept whole, and discards the fragment
when it is the query (if it is still there and unchanged by then).
After removal, subject sequences shorter than the threadshold will be discarded.
The number of rounds of iterations is automatically decided (usually less than 8). User can also define this.
Threads only run BLAST. The results are applied one query at a time in a fixed order (by the MD5 of the ID),
so the output does not depend on the number of threads or on timing.

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
-sameclass	[0|1]	A fragment of an entry of its own class: discarded by that entry (1, default), or cut out of it (0).
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
my $mb_wordsize = 20; #word size to salvage a query that detonates word_size 7 (times out): re-blast with megablast at this word size. Must stay >7 to avoid re-detonating on high-copy queries against a repetitive library.
my $offset = 7; #if two blast hits are less than $offset [default=7bp) away from each other, join them as one hit
my $clean = 1; #1, clean nested sequences; 0, will not clean nested, only discard highly overlapping (~100%) sequences
my $sameclass = 1; #1, a fragment of an entry of its own class is discarded by that entry instead of cut out of it; 0, it is cut out
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
	$mb_wordsize=$ARGV[$k+1] if /^-mb_wordsize$/i and $ARGV[$k+1] !~ /^-/;
	$clean=$ARGV[$k+1] if /^-clean$/i and $ARGV[$k+1] !~ /^-/;
	$sameclass=$ARGV[$k+1] if /^-sameclass$/i and $ARGV[$k+1] !~ /^-/;
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
die "\nERROR: The -sameclass parameter expects 0 or 1!\n$usage" unless $sameclass =~ /^[01]$/;
$blastplus = "" unless defined $blastplus;
$blastplus=`command -v blastn 2>/dev/null` if $blastplus eq '';
$blastplus=~s/blastn\n//i;
die "ERROR: blastn is not exist in the BLAST+ path $blastplus!\n" unless -X "${blastplus}blastn";

# shared state. Only the main thread changes %seq and %touched_seq (in apply_hits); the BLAST threads read them.
my %seq :shared;
my %touched_seq :shared; # here we will save all subject sequences that were modified after removing nested sequence
my $stop :shared = 0;    # set once -maxcount is reached: the BLAST threads skip what is left of the queue
my $count_stat = 0;      # stat lines so far (drives -maxcount and the saturation check)
my ($queue, $results);   # queries to BLAST [index, id, db]; BLAST results [index, id, hits]
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

	# The queries in a fixed order: by the MD5 of the ID (a fixed shuffle, so no group of IDs goes first),
	# then by ID. Which of two nested entries is the query first decides the outcome (the fragment cuts the
	# entry, or the entry discards the fragment), so this order, not thread timing, must decide it.
	my $N = 0;
	{
	my %key = map { $_ => md5_hex($_) } keys %seq;
	$queue = Thread::Queue->new();
	foreach my $id (sort { $key{$a} cmp $key{$b} or $a cmp $b } keys %key){
		next if $resume and $i == $start_iter and exists $done_q{$id}; # already processed before the interruption
		$queue -> enqueue([$N++, $id, "$IN.iter$i"]);
	}
	$queue -> end(); # signal that no more items will be added to the queue
	} # %key is freed here, before the threads copy the interpreter
	%done_q = () if $resume and $i == $start_iter;
	$results = Thread::Queue->new();
	$stop = 0;

	# initiate a number of worker threads to BLAST the queries
	my @threads;
	foreach (1..$threads){
		push @threads, threads->create(\&blaster);
	}

	# apply the BLAST results in queue order, whatever order they arrive in. A query modified or discarded
	# by an earlier one is skipped, as is a subject modified earlier in this iteration.
	open my $done_fh, ">>", "$IN.iter$i.done.main" or die "ERROR: cannot write $IN.iter$i.done.main: $!\n";
	$done_fh->autoflush(1); # per-iteration log of applied queries; a resumed run skips them
	my %ready;
	for my $n (0..$N-1){
		until (exists $ready{$n}){
			my $r = $results->dequeue();
			$ready{$r->[0]} = $r;
		}
		my ($id, $hits, $err) = @{delete $ready{$n}}[1, 2, 3];
		die "ERROR: BLAST thread failed on $id: $err\n" if defined $err;
		next unless defined $hits; # discarded or modified before its BLAST
		next unless exists $seq{$id} and !$touched_seq{$id};
		&apply_hits($id, $i, $hits);
		print $done_fh "$id\n";
		# When $count_stat reaches the user defined stat count, we end the entire queue.
		if ($count_limit > 0 and $count_stat >= $count_limit){
			print STAT "Reached user defined $count_limit of processed sequences at iter$i, stopping...\n\n";
			$stop = 1;
			last;
		}
	}
	close $done_fh;
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


# BLAST thread: search each query against the iteration's snapshot and pass the HSPs that pass the
# identity and length filters to the main thread. A query already discarded or modified is skipped.
sub blaster(){
	while (defined(my $job = $queue->dequeue())){
		my ($n, $id, $db) = @$job;
		my $s = $stop ? undef : $seq{$id};
		if (!defined $s or $touched_seq{$id} or length $s == 0){
			$results->enqueue([$n, $id, undef]);
			next;
		}
		my $hits = eval {
			my $seq = ">$id\n$s\n";
			my $exec="timeout -s KILL 120s ${blastplus}blastn -query <(echo -e \"$seq\") -db $db -word_size 7 -evalue 1e-5 -dust no -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen\"";
			my @Blast=qx(bash -c '$exec' 2> /dev/null);
			if (($? >> 8) == 137){ # word_size 7 detonated (timeout/OOM) on a high-copy query -> salvage with megablast (word_size $mb_wordsize)
				(my $exec_mb = $exec) =~ s/-word_size 7 -evalue 1e-5 -dust no/-word_size $mb_wordsize -evalue 1e-5/;
				@Blast=qx(bash -c '$exec_mb' 2> /dev/null);
			}
			my @hits;
			foreach (@Blast){
				my ($query, $subject, $iden, $len, $sbj_start, $sbj_end, $sbj_len) = (split)[0,1,2,3,8,9,11];
				next if grep { !defined } ($query, $subject, $iden, $len, $sbj_start, $sbj_end, $sbj_len);
				next if $query eq $subject;
				next if $iden < $min_iden;
				next if $len < $minlen; # the length of HSPs should be more than 80 bp
				push @hits, join("\t", $subject, $iden, $len, $sbj_start, $sbj_end, $sbj_len);
			}
			join("\n", @hits);
		};
		$results->enqueue(defined $hits ? [$n, $id, $hits] : [$n, $id, undef, $@ || 'unknown error']);
	}
}

# Main thread: apply one query's HSPs to the current sequences. A subject covering >= $coverage of the
# query, or covered >= $coverage by it, loses the matching regions and is discarded when shorter than
# $minlen. Between entries of the same class, the one that holds the other is kept whole and discards it.
sub apply_hits {
	my ($id, $i, $hits) = @_;
	my $length = length $seq{$id}; #query length
	my %merged_hsps;
	my %seq_len; #store subject seq length info
	my @aln_iden; # store alignment length and iden of this $id
	my $scaled_iden; # overall identity given all blast hits that pass the filter
	my $total_len; # total length of all alignments to $id
	foreach (split /\n/, $hits){
		my ($subject, $iden, $len, $sbj_start, $sbj_end, $sbj_len) = split /\t/;
		next unless exists $seq{$subject};
		next if $touched_seq{$subject}; # skip the iteration if the subject sequence was already modified (including discarded)
		($sbj_start, $sbj_end) = ($sbj_end, $sbj_start) if $sbj_start > $sbj_end;
		push @{$merged_hsps{$subject}}, [$sbj_start, $sbj_end];
		$seq_len{$subject} = $sbj_len;
		$total_len += $len;
		push @aln_iden, [$len, $iden];
	}
	return unless %merged_hsps;

	# calculate weighted identity
	foreach (@aln_iden) {
		my ($len, $iden) = @{$_};
		$scaled_iden += sprintf("%.3f", $iden * $len / $total_len);
	}

	# merge the overlapping HSPs of each subject, then remove the regions from the subject that are inserted
	# into the query
	foreach my $sbj (sort keys %merged_hsps){
		my ($ref1, $ref2) = &merger(@{$merged_hsps{$sbj}});
		my @hsps = @$ref1;
		my $merged = $$ref2; # number of overlapping HSPs
		my $size = 0; # total coverage by HSPs
		$size += $_->[1] - $_->[0] + 1 foreach @hsps;
		my $sbj_len = $seq_len{$sbj};
		my $seq_new = $seq{$sbj};
		next unless defined $seq_new;
		next if length $seq_new != $sbj_len; #if the subject length changes, it has been modified. Skip to avoid mismodification.
		my ($qcov, $scov) = (sprintf("%.3f", $size/$length), sprintf("%.3f", $size/$sbj_len));
		next unless $qcov >= $coverage or $scov >= $coverage;

		# One entry is a fragment of the other, of the same class, not an insertion in it: the entry stays
		# whole and discards the fragment when the entry is the query. Waiting for the entry's own turn means
		# a fragment goes only while the entry that holds it is still there unchanged.
		if ($sameclass and ($qcov >= $coverage xor $scov >= $coverage) and &same_class($id, $sbj)){
			next if $qcov >= $coverage; # the query is the fragment: leave the subject whole
			print STAT "$sbj\tIter$i\tDiscarded. Contained in $id of the same class; qcov: $qcov; scov: $scov; identity: $scaled_iden; merged $merged\n";
			$count_stat++;
			delete $seq{$sbj};
			$touched_seq{$sbj} = 1;
			return if $count_limit > 0 and $count_stat >= $count_limit;
			next;
		}

		# replace bases of HSPs regions to \0 (aka Remove); this masking is nessary since the subject sequence
		# may be cleaned several times, for each non-overlapping merged HSPs regions.
		my $poss = ''; # positions of the non-overlapped rHSPs egions that will be removed from the subject
		for my $hsp (@hsps) {
			my ($start, $end) = ($hsp->[0], $hsp->[1]);
			$poss = $poss . $start . ".." . $end . ",";
			my $len = $end - $start + 1;
			substr($seq_new, $start-1, $len) = "\0" x $len if length $seq_new >= $start + $len - 1;
		}
		$seq_new =~ tr/\0//d;
		my $sbj_len_new = length $seq_new;
		if ($sbj_len_new >= $minlen and $sbj_len_new < length $seq{$sbj} and $clean == 1){
			print STAT "$sbj\tIter$i\tCleaned. $poss covering $qcov of $id; scov: $scov; identity: $scaled_iden; merged $merged\n";
			$count_stat++;
			$seq{$sbj} = $seq_new; #overwrite the sbj sequence if the new one is shorter
			$touched_seq{$sbj} = 1; # this subject sequence was modifed, and we will not deal with it any more in the current iteration
		} elsif ($sbj_len_new < $minlen) {
			print STAT "$sbj\tIter$i\tDiscarded. Has only $sbj_len_new bp after cleaning by $id; qcov: $qcov; scov: $scov; identity: $scaled_iden; merged $merged\n";
			$count_stat++;
			delete $seq{$sbj}; #delete this sequence if new seq is too short
			$touched_seq{$sbj} = 1; # this subject sequence was modifed (removed), and we will not deal with it any more in the current iteration
		}
		return if $count_limit > 0 and $count_stat >= $count_limit;
	}
}

# true if both IDs carry the same class after the last # (case-insensitive); Unknown is no class
sub same_class {
	my ($c1, $c2) = map { /#([^#]+)$/ ? lc $1 : '' } @_;
	return ($c1 ne '' and $c1 ne 'unknown' and $c1 eq $c2);
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
		substr($seq_new, $start-1, $len) = "\0" x $len if length $seq_new >= $start + $len - 1;
	}
	$seq_new =~ tr/\0//d;
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

# union the processed-query logs for iteration $i into %$set: the main thread's ($IN.iter$i.done.main),
# and the per-thread ones of a run interrupted before 11/06/2026
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
