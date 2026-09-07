#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use threads::shared;

my $usage = "
	A script to purify a TE library based on another TE file containing the target contaminant.
	This is to use the richness difference between TE1 and TE2. Real contaminants is less abundant in TE1 but rich in TE2.
	Identified contaminated sequences will be converted into lowercases in the TE1-TE2.fa output.
		Usage: perl TE_purifier.pl -TE1 [fasta] -TE2 [fasta]
		options:	-TE1	[fasta]	The file to be purified.
				-TE2	[fasta]	The file that mainly consists of TE1 contaminants.
				-lower	[0|1]	Mask contaminants in TE1 with lowercase letters (1, default) or Ns (0).
				-minlen	[int]	The shortest length (bp) of sequence matches to be considered. Default: 50
				-miniden	[int]	The minimum identity (%) to be considered a real match. Default: 60
				-mindiff	[float]	The minimum fold difference in richness between TE1 and TE2 for a 
							sequence to be considered as real to TE1.
				-reprocess	[0|1]	Skip (1) RepeatMasking and use existing *stat file to regenerate the 
							TE1-TE2.fa file. Useful to test different -mindiff settings. default 0.
				-repeatmasker	[path]	The directory containing RepeatMasker (default: read from ENV)
				-blastplus	[path]	The directory containing Blastn (default: read from ENV)
				-threads	[int]	Number of theads to run this script
				-help|-h	Display this help info
\n";

# user input
my $TE1 = ""; #the file to be purified
my $TE2 = ""; #the file that has lots of $TE1 contaminants

# pre-defined
my $lower = 1; #use lower case (1, default) or Ns (0) to mask qualified contaminants
my $minlen = 50; #shortest length of match to be considered. I choose half the size of the shortest element (100bp) here.
my $miniden = 60; #minimum identity (%) to be considered a real match
my $mb_wordsize = 20; #word size used to salvage a query that detonates word_size 7 (times out): re-blast with megablast at this word size. Smaller = more sensitive (closer to word_size 7) but slower; must stay >7 to avoid re-detonating on high-copy queries against a repetitive library.
my $mindiff = "0.4"; #minimum richness difference between $TE1 and $TE2 for a sequence to be considered as real to $TE1
my $reprocess = 0; #skip (1) RepeatMasking and use existing *stat file to regenerate the TE1-TE2.fa file. default 0.
my $script_path = $FindBin::Bin;
my $call_seq = "$script_path/call_seq_by_list.pl";
my $repeatmasker = "";
my $blastplus = "";
my $threads = 8;
my $batch_size = 1000; #masked regions per batched blastn invocation

# read parameters
my $k=0;
foreach (@ARGV){
	$TE1 = $ARGV[$k+1] if /^-TE1$/i and $ARGV[$k+1] !~ /^-/;
	$TE2 = $ARGV[$k+1] if /^-TE2$/i and $ARGV[$k+1] !~ /^-/;
	$lower = $ARGV[$k+1] if /^-lower$/i and $ARGV[$k+1] !~ /^-/;
	$minlen = $ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$miniden = $ARGV[$k+1] if /^-miniden/i and $ARGV[$k+1] !~ /^-/;
	$mb_wordsize = $ARGV[$k+1] if /^-mb_wordsize/i and $ARGV[$k+1] !~ /^-/;
	$mindiff = $ARGV[$k+1] if /^-mindiff/i and $ARGV[$k+1] !~ /^-/;
	$reprocess = $ARGV[$k+1] if /^-reprocess/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^-blastplus/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
        }

# some checks
die "The TE1 file $TE1 is not found!\n$usage" unless -e $TE1;
die "The TE2 file $TE2 is not found!\n$usage" unless -e $TE2;

# read $TE1 into memory and count total length
my ($TE1_len, $TE2_len) = (0, 0);
open TE1, "<$TE1" or die $!;
open TE2, "<$TE2" or die $!;
my %TE1;
$/ = "\n>";
while (<TE1>){
	chomp;
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id = (split)[0];
	next if length $id > 80;
	$seq =~ s/\s+//g;
	$seq = uc $seq; # convert all TE1 sequences into uppercase
	my $len = length $seq;
	$TE1_len += $len;
	$TE1{$id} = $seq;
	}

# count $TE2 total length
while (<TE2>){
	chomp;
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my $len = length $seq;
	$TE2_len += $len;
	}
close TE1;
close TE2;
$/ = "\n";

# empty libraries are not fatal: an empty TE1 (query) library produces empty outputs; an empty
# TE2 (masking) library has nothing to mask with, so TE1 passes through unmodified
if (!%TE1 or $TE2_len == 0){
	warn "WARNING: The TE1 library $TE1 contains no sequence. The output will be empty.\n" unless %TE1;
	warn "WARNING: The TE2 library $TE2 contains no sequence. Nothing to mask with: TE1 passes through unmodified.\n" if %TE1;
	open my $out, ">", "$TE1-$TE2.fa" or die $!;
	foreach my $id (sort {$a cmp $b} keys %TE1){
		print $out ">$id\n$TE1{$id}\n";
		}
	close $out;
	open my $stat, ">", "$TE1-$TE2.stat" or die $!;
	print $stat "TE1_id\tTE1_len\tTE2_len\tTE1_richness\tTE2_richness\tFold_diff\n";
	close $stat;
	exit 0;
	}


if ($reprocess == 0){

# Repeatmask TE1 with TE2; make blast db for $TE1 and $TE2
my $div = 100 - $miniden;
my $err = '';
unlink "$TE1.out"; # drop a stale .out from a killed earlier run so it can't be mistaken for this run's result
$err = `${repeatmasker}RepeatMasker -e ncbi -pa $threads -qq -no_is -nolow -div $div -lib $TE2 $TE1 > ${TE2}-${TE1}.RM.status` // 'no output captured';
die "ERROR: RepeatMasker failed ($?): $err\n" if $? != 0;
my $mbdb1 = `${blastplus}makeblastdb -in $TE1 -out $TE1 -dbtype nucl 2>&1` // 'no output captured';
die "ERROR: makeblastdb failed on $TE1 ($?): $mbdb1\n" if $? != 0;
my $mbdb2 = `${blastplus}makeblastdb -in $TE2 -out $TE2 -dbtype nucl 2>&1` // 'no output captured';
die "ERROR: makeblastdb failed on $TE2 ($?): $mbdb2\n" if $? != 0;
print STDERR "$err\n" if $err ne '';

# get masked regions of TE1
if (`grep "No repetitive sequences were detected" ${TE2}-${TE1}.RM.status`){
	print STDERR "RepeatMasker ran correctly. No repetitive sequences were detected.\n\n";
	`touch $TE1.out`
	}
open RM, "<$TE1.out" or die $!;
my %RM;
while (<RM>){
	s/^\s+//;
	my ($SW_score, $id, $from, $to) = (split)[0,4,5,6];
	next unless defined $SW_score;
	next unless $SW_score =~ /^[0-9]+$/;
	next if $SW_score < 300;
	my $len = abs ($to - $from) + 1;
	($from, $to) = ($to, $from) if $to < $from;
	next if $len < $minlen;
	$RM{$id} .= "$from-$to ";
	}
close RM;

###############
# Identify fold difference between TE1 and TE2
my $stat_lock :shared; #serializes writes to STAT so lines never interleave
open STAT, ">$TE1-$TE2.stat" or die $!;
print STAT "TE1_id\tTE1_len\tTE2_len\tTE1_richness\tTE2_richness\tFold_diff\n";

# collect every masked region as one query; the plain serial id "r<N>" cannot be altered by
# BLAST id parsing (unlike ids containing "|" or ":"), %region maps it back to the region
my (%region, @queries);
foreach my $id (sort {$a cmp $b} keys %RM){
	while ($RM{$id} =~ s/([0-9]+)-([0-9]+)//){
		my ($from, $to, $seqlen) = ($1, $2, $2-$1+1);
		next unless exists $TE1{$id};
		my $seq = substr $TE1{$id}, $from-1, $seqlen;
		next unless defined $seq;
		my $qid = "r".scalar(@queries);
		$region{$qid} = [$id, $from, $to];
		push @queries, [$qid, $seq];
		}
	}

# blast the regions in batches of $batch_size: ONE blastn against TE1 and ONE against TE2 per
# batch (each -num_threads $threads) replaces the two single-threaded blastn processes per
# region of the original implementation
while (@queries){
	my @batch = splice @queries, 0, $batch_size;
	my $rich_te1 = &batch_richness(\@batch, $TE1);
	my $rich_te2 = &batch_richness(\@batch, $TE2);

	#calculate the fold difference in richness. the smaller the more likely it belongs to $TE2 (contaminant of $TE1)
	foreach my $query (@batch){
		my ($qid) = @$query;
		my ($id, $from, $to) = @{$region{$qid}};
		my ($seq_te1_len, $seq_te2_len) = ($rich_te1->{$qid}, $rich_te2->{$qid});
		my ($seq_te1_percent, $seq_te2_percent) = ($seq_te1_len/$TE1_len, $seq_te2_len/$TE2_len);
		my $diff = 1000;
		$diff = $seq_te1_percent/$seq_te2_percent if $seq_te2_percent > 0;
		{ lock($stat_lock); print STAT "$id:$from..$to\t$seq_te1_len\t$seq_te2_len\t$seq_te1_percent\t$seq_te2_percent\t$diff\n"; }
		}
	}
close STAT;
###############
}


# filter sequences based on min_diff
open STAT, "<$TE1-$TE2.stat" or die $!;
while (<STAT>){
	next if /^TE1_id\s+/;
	my ($info, $diff) = (split)[0,5];
	unless (defined $info and defined $diff and $info =~ /(.*):([0-9]+)\.\.([0-9]+)/){
		warn "WARNING: skip malformed stat line $. in $TE1-$TE2.stat: $_";
		next;
		}
	my ($id, $from, $to, $seqlen) = ($1, $2, $3, $3-$2+1);
	unless (exists $TE1{$id}){
		warn "WARNING: skip stat line for unknown TE1 id $id\n";
		next;
		}
	my $ori_seq = $TE1{$id};
	my $seq = substr $ori_seq, $from-1, $seqlen;

	# convert contaminated TE1 sequences into lowercase
	substr ($ori_seq, $from-1, $seqlen) = lc $seq if $diff < $mindiff;
	$TE1{$id} = $ori_seq;
	}


# output unprocessed and clean sequence
open Seq, ">$TE1-$TE2.fa" or die $!;
foreach my $id (sort {$a cmp $b} keys %TE1){
	print Seq ">$id\n$TE1{$id}\n";
	}
close STAT;
close Seq;


# batch_richness(\@batch, $db): the richness (see richness()) of every query of the batch
# against $db, computed with ONE multi-threaded blastn instead of one blastn process per query.
# If some query of the batch detonates word_size 7 and the batch blastn is KILLed, the whole
# batch falls back to the original per-query richness() (188s timeout + megablast salvage).
sub batch_richness {
	my ($batch, $db) = @_;
	my %sum;
	$sum{$_->[0]} = 0 foreach @$batch;
	my $query_file = "$TE1-$TE2.query.tmp";
	open Q, ">$query_file" or die $!;
	foreach my $query (@$batch){
		print Q ">$query->[0]\n$query->[1]\n";
		}
	close Q;
	my $budget = 188 * int((scalar(@$batch) + $threads - 1)/$threads); #the old per-query timeout, scaled to the batch
	my $exec = "timeout -s KILL ${budget}s ${blastplus}blastn -db $db -query $query_file -outfmt 6 -word_size 7 -evalue 1e-5 -dust no -num_threads $threads";
	my $out = qx($exec 2> /dev/null);
	if (($? >> 8) == 137){
		foreach my $query (@$batch){
			$sum{$query->[0]} = &richness($query->[1], $db);
			}
	} else {
		warn "WARNING: batched blastn vs $db exited with status $?, results may be incomplete\n" if $? != 0;
		foreach my $row (split /\n/, $out){
			my ($qid, $iden, $len) = (split ' ', $row)[0,2,3];
			next unless defined $len and $iden =~ /^[0-9]/ and $len =~ /^[0-9]/;
			$qid =~ s/^lcl\|//; #some BLAST+ versions prefix local ids with lcl|
			next if $iden < $miniden or $len < $minlen;
			$sum{$qid} += $len;
			}
		}
	unlink $query_file;
	return \%sum;
	}



# richness($seq, $db): total length of qualified matches of $seq in $db (its "richness").
# Default to the sensitive word_size 7, which reproduces the original count exactly. If word_size 7
# detonates on a high-copy query and is KILLed (times out), salvage the count with megablast
# (word_size $mb_wordsize), which completes fast. Non-detonating queries are therefore counted
# identically to the original; only the few detonators -- which the original could not finish at
# all -- differ.
sub richness {
	my ($seq, $db) = @_;
	my ($out, $rc) = &run_blast($seq, $db, "-word_size 7 -evalue 1e-5 -dust no");
	($out, $rc) = &run_blast($seq, $db, "-word_size $mb_wordsize -evalue 1e-5") if ($rc >> 8) == 137;
	my $sum = 0;
	foreach (@$out){
		my ($iden, $len) = (split)[2,3];
		next unless defined $len and $iden =~ /^[0-9]/ and $len =~ /^[0-9]/;
		next if $iden < $miniden or $len < $minlen;
		$sum += $len;
	}
	return $sum;
}

# run_blast($seq, $db, $opts): one blastn with a hard timeout, returning (\@rows, $rc). Retry only
# transient failures, never a KILL (timeout/OOM) -- a detonating query then costs one timeout
# instead of a 10x retry storm.
sub run_blast {
	my ($seq, $db, $opts) = @_;
	my $exec = "timeout -s KILL 188s ${blastplus}blastn -db $db -query <(echo -e \"$seq\") -outfmt 6 $opts";
	my (@out, $rc);
	for (my $try = 0; $try < 3; $try++){
		@out = qx(bash -c '$exec' 2> /dev/null);
		$rc = $?;
		last if $rc == 0;
		last if ($rc >> 8) == 137;
	}
	return (\@out, $rc);
}
