#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;
use threads;
use Thread::Queue;
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
my $mindiff = "0.4"; #minimum richness difference between $TE1 and $TE2 for a sequence to be considered as real to $TE1
my $reprocess = 0; #skip (1) RepeatMasking and use existing *stat file to regenerate the TE1-TE2.fa file. default 0.
my $script_path = $FindBin::Bin;
my $call_seq = "$script_path/call_seq_by_list.pl";
my $repeatmasker = "";
my $blastplus = "";
my $threads = 36;
my $queue = Thread::Queue->new();

# read parameters
my $k=0;
foreach (@ARGV){
	$TE1 = $ARGV[$k+1] if /^-TE1$/i and $ARGV[$k+1] !~ /^-/;
	$TE2 = $ARGV[$k+1] if /^-TE2$/i and $ARGV[$k+1] !~ /^-/;
	$lower = $ARGV[$k+1] if /^-lower$/i and $ARGV[$k+1] !~ /^-/;
	$minlen = $ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$miniden = $ARGV[$k+1] if /^-miniden/i and $ARGV[$k+1] !~ /^-/;
	$mindiff = $ARGV[$k+1] if /^-mindiff/i and $ARGV[$k+1] !~ /^-/;
	$reprocess = $ARGV[$k+1] if /^-reprocess/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^-blastplus/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
        }

# some checks
die "The TE1 file $TE1 is not found or it's empty!\n$usage" unless -s $TE1;
die "The TE2 file $TE2 is not found or it's empty!\n$usage" unless -s $TE2;

# define RepeatMasker -pa parameter
my $rm_threads = int($threads/4);

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


if ($reprocess == 0){

# Repeatmask TE1 with TE2; make blast db for $TE1 and $TE2
my $div = 100 - $miniden;
my $err = '';
$err = `${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -qq -no_is -nolow -div $div -lib $TE2 $TE1 >/dev/null`;
`${blastplus}makeblastdb -in $TE1 -out $TE1 -dbtype nucl 2> /dev/null`;
`${blastplus}makeblastdb -in $TE2 -out $TE2 -dbtype nucl 2> /dev/null`;
print STDERR "$err\n" if $err ne '';

# get masked regions of TE1
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
open STAT, ">$TE1-$TE2.stat" or die $!;
print STAT "TE1_id\tTE1_len\tTE2_len\tTE1_richness\tTE2_richness\tFold_diff\n";

# multi-threading using queue, put candidate regions into queue for parallel computation
$queue = Thread::Queue->new();
foreach my $id (keys %RM){
	last unless defined $RM{$id};
	$queue -> enqueue([$id, $RM{$id}]);
	}
$queue -> end();

# initiate a number of worker threads and run
foreach (1..$threads){
	threads->create(\&purifier);
	}
foreach (threads -> list()){
	$_->join();
	}
close STAT;
###############
}


# filter sequences based on min_diff
open STAT, "<$TE1-$TE2.stat" or die $!;
while (<STAT>){
	next if /^TE1_id\s+/;
	my ($info, $diff) = (split)[0,5];
	my ($id, $from, $to, $seqlen) = ($1, $2, $3, $3-$2+1) if $info =~ /(.*):([0-9]+)\.\.([0-9]+)/;
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

# fixing the formatting error created by simutaniously writing the same file
`perl -i -nle 's/>/\\n>/g unless /^>/; print \$_' $TE1-$TE2.fa`;


# subrotine for the purifier
sub purifier(){
	while (defined($_ = $queue->dequeue())){
		my ($id, $coor) = (@{$_}[0], @{$_}[1]);
		next unless exists $TE1{$id};
		my $ori_seq = $TE1{$id};

		while ($coor =~ s/([0-9]+)-([0-9]+)//){
			my ($from, $to, $seqlen) = ($1, $2, $2-$1+1);
			my $seq = substr $ori_seq, $from-1, $seqlen;
			next unless defined $seq;
			my $target = ">$id:$from..$to\\n$seq";

			# count the size of $target in $TE1
			my $exec = "timeout 188s ${blastplus}blastn -db $TE1 -query <(echo -e \"$seq\") -outfmt 6 -word_size 7 -evalue 1e-5 -dust no";
			my @blast_te1 = ();
			my $try = 0;
			while ($try < 10){ #try 10 times to guarantee the blast is run correctly
				@blast_te1 = qx(bash -c '$exec' 2> /dev/null) if defined $seq;
				last if $? == 0;
				$try++;
				}
			my $seq_te1_len = 0;
			foreach (@blast_te1){
				my ($id, $iden, $len) = (split)[1,2,3];
				next unless (defined $id && defined $iden && defined $len);
				next if $iden < $miniden or $len < $minlen;
				$seq_te1_len += $len;
				}

			# count the size of $target in $TE2
			$exec = "timeout 188s ${blastplus}blastn -db $TE2 -query <(echo -e \"$seq\") -outfmt 6 -word_size 7 -evalue 1e-5 -dust no";
			my @blast_te2 = ();
			$try = 0;
			while ($try < 10){ #try 10 times to guarantee the blast is run correctly
				@blast_te2 = qx(bash -c '$exec' 2> /dev/null) if defined $seq;
				last if $? == 0;
				$try++;
				}
			my $seq_te2_len = 0;
			foreach (@blast_te2){
				my ($id, $iden, $len) = (split)[1,2,3];
				next if $iden < $miniden or $len < $minlen;
				$seq_te2_len += $len;
				}

			#calculate the fold difference in richness. the smaller the more likely it belongs to $TE2 (contaminant of $TE1)
			my ($seq_te1_percent, $seq_te2_percent) = ($seq_te1_len/$TE1_len, $seq_te2_len/$TE2_len);
			my $diff = 1000;
			$diff = $seq_te1_percent/$seq_te2_percent if $seq_te2_percent > 0;
			print STAT "$id:$from..$to\t$seq_te1_len\t$seq_te2_len\t$seq_te1_percent\t$seq_te2_percent\t$diff\n";
			}
		}
	}

