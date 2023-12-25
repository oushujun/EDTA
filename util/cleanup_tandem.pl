#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

my $usage="
        Usage: perl cleanup_tandem.pl -f sample.fa [options] > sample.cln.fa 
	Options:
		-misschar	[n|l]	Define the letter representing unknown sequences; default: n. l: recognize lower case letters
		-Nscreen	[0|1]	Enable (1) or disable (0) the -nc parameter; default: 1
		-nc		[int]	Ambuguous sequence len cutoff; discard the entire sequence if > this number; default: 0
		-nr		[0-1]	Ambuguous sequence percentage cutoff; discard the entire sequence if > this number; default: 1
		-minlen		[int]	Minimum sequence length filter after clean up; default: 100 (bp)
		-maxlen		[int]	Maximum sequence length filter after clean up; default: 25000 (bp)
		-cleanN		[0|1]	Retain (0) or remove (1) the -misschar taget in output sequence; default: 0
		-cleanT		[0|1]	Remove entire seq. if any terminal seq (default: 20bp, set by -Tlen) has 75% of N (-cleanT 1); disabled by default (-cleanT 0).
		-Tlen		[int]	The length of terminal sequence for the -cleanT parameter; default: 20 (bp)
		-minrm		[int]	The minimum length of -misschar to be removed if -cleanN 1; default: 1.
		-trf		[0|1]	Enable (1) or disable (0) tandem repeat finder (trf); default: 1
		-trf_path	path	Path to the trf program
        \n";
my $version="
cleanup_tandem.pl
Clean up tandem repeat sequence and N gaps in the fasta file
Author: Shujun Ou (shujun.ou.1 @ gmail.com)
Version:	1.7 Add the Tlen parameter. 2023/12/23
		1.6 Add alignment score control (-e [int]) and missing count control (-c [int])
		1.5 Add missing rate control (-r [0, 1]) and modify missing count control (-nc not control)
		1.0 2014/06/02
\n";

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

my $target="n";
my $func_nc=1; #1, do $n_count screening
my $n_count=0; #count the $target in each sequence, if it exceeds this number, it will be discarted.
my $n_rate=1; #count the $target in each sequence, if it exceeds this percentage, it will be discarted.
my $minlen=100; #Minimum sequence length filter after clean up; default: 100 (bp)
my $maxlen=25000; #Maximum sequence length filter after clean up; default: 25000 (bp)
my $minrm=1; #The minimum length of -misschar to be removed if -cleanN 1
my $align_score=1000; #-e para, dft:1000 
my $max_seed=2000; #maximum period size to report
my $cleanN=0; #1 will remove $target="n" in output sequence
my $trf=1; #1 will enable tandem repeat finder (default), 0 will not
my $cleanT=0; #1 will distard the entire sequence if any terminal sequence (head/tail) contains >= 75% of $target="n"; #0 will not do so.
my $Tlen=20; #The length of terminal sequence for the -cleanT parameter; default: 20 (bp)
my $trf_path='';
my $file;

my $k=0;
foreach (@ARGV){
	$target=$ARGV[$k+1] if /^-misschar$/i;
	$func_nc=0 if /^-Nscreen$/i;
	$n_count=int($ARGV[$k+1]) if /^-nc$/i;
	$n_rate=$ARGV[$k+1] if /^-nr$/i;
	$minlen=int($ARGV[$k+1]) if /^-minlen$/i;
	$maxlen=int($ARGV[$k+1]) if /^-maxlen$/i;
	$minrm=$ARGV[$k+1] if /^-minrm$/i;
	$align_score=$ARGV[$k+1] if /^-minscore$/i;
	$file=$ARGV[$k+1] if /^-f$/i;
	$cleanN=$ARGV[$k+1] if /^-cleanN$/i;
	$cleanT=$ARGV[$k+1] if /^-cleanT$/i;
	$Tlen=int($ARGV[$k+1]) if /^-Tlen$/i;
	$trf=$ARGV[$k+1] if /^-trf$/i;
	$trf_path=$ARGV[$k+1] if /^-trf_path$/i;
	$k++;
	}

# check TRF
$trf_path=`which trf 2>/dev/null` if $trf_path eq '';
$trf_path=~s/\n$//;
`$trf_path 2>/dev/null`;
die "Error: No Tandem Repeat Finder is working on the current system.
	Please report it to https://github.com/oushujun/EDTA/issues" if $?==32256;
die "\n\tTandem Repeat Finder not found!\n\n$usage" unless $trf_path ne '';
die "\n\tInput file \"$file\" not found!\n\n$usage" unless -e $file;
die "\n\tTlen must be > 0!\n\n$usage" unless $Tlen > 0;

my %tandem;
my $tandem='';
$tandem=`$trf_path $file 2 7 7 80 10 $align_score $max_seed -ngs -h -l 6` if $trf==1;
while ($tandem=~s/\@(.*)\n?// and $tandem ne ''){
	$tandem{$1}=$1;
	}

open Info, ">$file.cleanup";
open File, "<$file" or die "ERROR: $!";
$/="\n>";
while (<File>){
	if (/^$/){next}
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	my $length=length $seq;
	my $mark=0;
	my $count=0;
	if ($target =~ /n/i){ #hardmask N
		$count++ while $seq=~/$target/gi;
		}
	if ($target =~ /^l$/i){ #softmask lowercase
		$count++ while $seq=~/[atcgnN]/g;
		}
	my $count_rate=$count/$length;

#missing control
	if ($count_rate>=$n_rate){
		print Info "$id\t$count_rate missing\n";
		$mark=1;
		}
	elsif ($count>=$n_count && $func_nc==1 && $n_count>0){
		print Info "$id\t$count sequence gap\n";
		$mark=1;
		}

#tandem sequence control
	if (exists $tandem{$id}){
		print Info "$id\ttandem sequence\n";
		$mark=1;
		}

#terminal missing control
	if ($cleanT==1){
		my $start_20=substr $seq, 0, $Tlen;
		my $end_20=substr $seq, -$Tlen;
		my ($count_s, $count_e) = (0,0);
		if ($target =~ /n/i){ #hardmask
			$count_s++ while $start_20=~/$target/gi;
			$count_e++ while $end_20=~/$target/gi;
			}
		elsif ($target =~ /^l$/i){ #softmask
			$count_s++ while $start_20=~/[atcgnN]/g;
			$count_e++ while $end_20=~/[atcgnN]/g;
			}
		(print Info "$id\tSequence head has $count_s/$Tlen bp missing\n" and $mark=1) if $count_s/$Tlen>=0.75;
		(print Info "$id\tSequence tail has $count_e/$Tlen bp missing\n" and $mark=1) if $count_e/$Tlen>=0.75;
		}

#remove missing seq and length control
	if ($cleanN==1){
#		$seq=~s/$target//gi if $target =~ /n/i;
		$seq=~s/([nN]{$minrm,})//gi if $target =~ /n/i; #hardmask
		$seq =~ s/([atcgnN]{$minrm,})//g if $target =~ /^l$/i; #softmask
		my $len=length $seq;
		(print Info "$id\tOnly $len bp left after cleanup\n" and $mark=1) if $len < $minlen;
		(print Info "$id\tSequence $len bp > the limit $maxlen bp after cleanup\n" and $mark=1) if $len > $maxlen;
		}

	$seq = uc $seq;
	print ">$id\n$seq\n" unless $mark==1;
	}
close Info;
close File;
