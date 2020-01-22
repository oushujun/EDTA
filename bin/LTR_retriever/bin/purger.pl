#!/usr/bin/env perl -w
use strict;
use File::Basename;

my $usage="
Clean up sequence using blast resuls

perl purger.pl -blast blast_outfmt6 -seq seq.fa [options]

Options:
	-eval	[0-1]	e-value cutoff; discard the hit if >= this number; default 0.001
	-len	[int]	length cutoff; discard the hit if < this number; default 90 (bp)
	-iden	[0-100]	identity cutoff; discard the hit if <= tis number; default 35 (%)
	-cov	[0-1]	coverage cutoff; discard the entire sequence if >= this number; default 1
	-purge	[0|1]	purge switch; switch on=1(default)/off=0 to clean up aligned region and joint unaligned sequences

Dependency: combine_overlap.pl, call_seq_by_list.pl

BLAST example:
blastn -subject seq.fa -query removal.lib.fa -outfmt=6 > blast_outfmt6

Shujun Ou (oushujun\@msu.edu)
04/17/2017
\n";


#take blast outfmt=6 output, with each column means:
#query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

#use Cwd 'abs_path';
my $script_path = dirname(__FILE__);

my $seq; #provide the sequence to be purged
my $blast; #provide the blast outfmt=6 result
my $evalue=0.001; #evalue cutoff for blast entries. Evalues lower than this cutoff is considered a real alignment.
my $length=90; #length cutoff (bp, default 90) to be considered as a real alignment (= alignment length - mismatch)
my $identity=35; #identity cutoff (%, default 35) to be considered as a read alignment
my $purge=1; #switch on=1(default)/off=0 to clean up aligned region and joint unaligned sequences
my $coverage=1; #if the excluded portion is too long (default 1, [0-1]), discard the entire sequence

my $k=0;
foreach my $para (@ARGV){
	$seq=$ARGV[$k+1] if $para=~/^-seq$/i;
	$blast=$ARGV[$k+1] if $para=~/^-blast$/i;
	$evalue=$ARGV[$k+1] if $para=~/^-eval$/i;
	$length=$ARGV[$k+1] if $para=~/^-len$/i;
	$coverage=$ARGV[$k+1] if $para=~/^-cov$/i;
	$identity=$ARGV[$k+1] if $para=~/^-iden$/i;
	$purge=$ARGV[$k+1] if $para=~/^-purge$/i;
	$k++;
	}

open File, "<$blast" or die "ERROR: Please specify the BLAST result!\n$usage";
my %query; #store query information
my $info='';
while (<File>){
	s/^\s+//;
	my ($query, $iden, $len, $mismatch, $qstart, $qend, $eval)=(split)[0,2,3,4,6,7,10];
	($qstart, $qend)=($qend, $qstart) if $qstart>$qend;
	$info.="$query\t$qstart\t$qend\n" if ($eval<=$evalue and $len-$mismatch>=$length and $iden>=$identity);
}
$info="Good news! No sequence is needed to be purged.\n" if $info=~/^(\s+)?$/;
open Out, ">$seq.exclude.temp";
print Out "$info";
close Out;

if ($info=~/Good news!/i){
	`mv $seq.exclude.temp $seq.exclude.list`;
	`cp $seq $seq.clean`;
	} else {
	`perl $script_path/combine_overlap.pl $seq.exclude.temp $seq.exclude.list`;
	`rm $seq.exclude.temp`;
	`awk '{print \$1\"\\t\"\$1\":\"\$2\"..\"\$3}' $seq.exclude.list | perl $script_path/call_seq_by_list.pl - -C $seq -ex -cov $coverage -purge $purge > $seq.clean`;
	}
