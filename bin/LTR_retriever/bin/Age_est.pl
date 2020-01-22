#!/usr/bin/env perl -w
use strict;
use File::Basename;

my $usage="
Estimate the mean identity of LTR regions in the genome via all-versus-all BLAST
Usage: perl Age_est.pl -blast [path to the blastn program] -RMout [Repeatmasker .out file] -genome [genome file] [*options]
Options:
	-q	quick estimation of LTR identity (much faster for large genomes, may sacrifice ~0.5% of accuracy)
	-t	number of threads to run blastn
	-h	print this usage message
Credit: Shujun Ou (oushujun\@msu.edu) 02/07/2018
\n";

my $genome=''; #genome file
my $RMout=''; #Repeatmasker .out file generated using the LTR_retriever library
my $quick=0; #quick estimation of LTR identity (may sacrifice ~0.5% of accuracy)
my $read_coverage=0.8; #as an LTR-derived read, at least 80% of the read should be mapped to LTR
my $evalue=0.0001; #maximum e-value for a real hit
my $iden_cut=100; #alignment identity cutoff, hits higher than this value will be filtered out. Default: 100 (%)
my $threads=4;
my $sample=$ARGV[0];

#paths to dependent programs
my $blast=''; #blastn and makeblastdb
my $script_path = dirname(__FILE__);

my $k=0;
foreach (@ARGV){
	$blast=$ARGV[$k+1] if /^-blast$/i;
	$RMout=$ARGV[$k+1] if /^-RMout$/i;
	$genome=$ARGV[$k+1] if /^-genome$/i;
	$quick=1 if /^-q$/i;
	$threads=$ARGV[$k+1] if /^-t$/i;
	$iden_cut=$ARGV[$k+1] if /^-iden_cut$/i;
	die $usage if /^-?-h|help$/i;
        $k++;
        }

$blast=~s/blastn$//;
$blast.='/' unless $blast=~/\/$/ or $blast eq '';
die $usage unless -s "$RMout" and -s "$genome";

open RMout, "<$RMout" or die $!;
my @list;
while (<RMout>){
	next unless /_LTR/;
	s/^\s+//;
	my ($SW, $chr, $from, $to, $id)=(split)[0,4,5,6,9];
	next unless $SW =~ /^[0-9]+$/;
	next if $SW < 300 or $to-$from+1 < 100;
	push @list, "$chr:$from..$to\t$chr:$from..$to";
	}
close RMout;

my $n_all = $#list +1;
my $iden_all;
if ($quick==1 and $n_all > 60000){ #at least exceed 60K sequences to execute the quick estimation module, otherwise, no need.
#for quick estimation, randomly choose 20K, 40K, and 60K sequences
	my ($n1, $n2, $n3)=(20000, 40000, 60000);

#work on the first set
	my $out="$RMout.$n1";
	open List10, ">$out.LAI.LTRlist" or die $!;
	for (my $i=1; $i<=$n1; $i++){
		my $k=int(rand($#list));
		print List10 "$list[$k]\n";
		}
	close List10;
	my $iden_10K=&Age_est($out);
	`rm $out.LAI.LTRlist $out.LAI.LTR.fa*`;

#work on the second set
	$out="$RMout.$n2";
	open List20, ">$out.LAI.LTRlist" or die $!;
	for (my $i=1; $i<=$n2; $i++){
		my $k=int(rand($#list));
		print List20 "$list[$k]\n";
		}
	close List20;
	my $iden_20K=&Age_est($out);
	`rm $out.LAI.LTRlist $out.LAI.LTR.fa*`;

#work on the third set
	$out="$RMout.$n3";
	open List30, ">$out.LAI.LTRlist" or die $!;
	for (my $i=1; $i<=$n3; $i++){
		my $k=int(rand($#list));
		print List30 "$list[$k]\n";
		}
	close List30;
	my $iden_30K=&Age_est($out);
	`rm $out.LAI.LTRlist $out.LAI.LTR.fa*`;

#use the Lease Squares Estimate method to estimate the slope k in iden = iden0 + k * lg(n)
	my $n_mean=(log($n1)+log($n2)+log($n3))/3;
	my $i_mean=($iden_10K+$iden_20K+$iden_30K)/3;
	my $k = ((log($n1)-$n_mean)*($iden_10K-$i_mean) + (log($n2)-$n_mean)*($iden_20K-$i_mean) + (log($n3)-$n_mean)*($iden_30K-$i_mean))/((log($n1)-$n_mean)**2 + (log($n2)-$n_mean)**2 + (log($n3)-$n_mean)**2);
	$iden_all = $iden_30K + $k * (log($n_all) - log($n3));

	} else {

	open List, ">$RMout.LAI.LTRlist" or die $!;
	print List "$_\n" foreach @list;
	close List;
	$iden_all=&Age_est($RMout);
	`rm $RMout.LAI.LTRlist $RMout.LAI.LTR.fa*`;
	}

#print out the result and write on a file
open Age, ">$RMout.LAI.LTR.ava.age" or die $! if $quick != 1;
open Age, ">$RMout.q.LAI.LTR.ava.age" or die $! if $quick == 1;
print "Input:$RMout\tSeq_num:$n_all\tMean_identity:$iden_all\n";
print Age "Input:$RMout\tSeq_num:$n_all\tMean_identity:$iden_all\n";
close Age;



#the subroutine to estimate mean identity of a given sequence set
sub Age_est {
my $RMout=$_[0];
`perl $script_path/call_seq_by_list.pl $RMout.LAI.LTRlist -C $genome > $RMout.LAI.LTR.fa`;

die "$RMout.LAI.LTR.fa is empty, please check the $genome file and the $RMout.LAI.LTRlist file\n" unless -s "$RMout.LAI.LTR.fa";

`${blast}makeblastdb -in $RMout.LAI.LTR.fa -dbtype nucl`;
`${blast}blastn -word_size 20 -outfmt 6 -evalue 0.0001 -num_alignments 10 -num_threads $threads -query $RMout.LAI.LTR.fa -db $RMout.LAI.LTR.fa -out $RMout.LAI.LTR.ava.out`;

my ($read1, $read2)=('', '');
my ($identity1, $identity2, $total_iden, $read_num)=(0,0,0,0);

open Blast, "<$RMout.LAI.LTR.ava.out" or die "$RMout.LAI.LTR.ava.out is empty, the blastn program may not run correctly.\n";
while (<Blast>){
	my ($chr, $alignment_length, $gap, $eval);
	s/^\s+//;
	($read1, $chr, $identity1, $alignment_length, $gap, $eval)=(split)[0,1,2,3,5,10];
#query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	next if $identity1 > $iden_cut; #filter out hits with identity higher than $iden_cut
	next if $read1 eq $chr; #exclude self alignments
	my $read_length=$2-$1+1 if $read1=~/:([0-9]+)\.\.([0-9]+)\|/; #obtain the length of the query from the query ID
	next if $eval>$evalue;
	next if $alignment_length/$read_length < $read_coverage or $alignment_length < 100; #filter out low-quality hits
	if ($read2 eq $read1){
		$identity2=$identity1 if $identity1>$identity2;
		} else {
		$total_iden+=$identity2;
		$read_num++;
		$read2=$read1;
		$identity2=$identity1;
		}
	}
my $average_iden=$total_iden/$read_num;
return $average_iden;
}
