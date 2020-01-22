#!/usr/bin/env perl -w
use strict;

#Description: This is the script to calculate the LTR-RT distribution of different superfamilies
#Author: Shujun Ou (oushujun@msu.edu)
#Last updated: 05/10/2018
#Note: Memory consumption of this scrip is approx. 4X the size of the input genome

my $usage= "perl LTR_sum.pl -genome genome.fa -all genome.fa.RM.out [options]
	-window [int]	bp size of the sliding window, default 3,000,000
	-step [int]	bp size of the moving step, defalut 300,000
	-intact		indicate the -all file is an LTR_retriever .pass.list instead of a RepeatMasker .out file";

my $window="3000000"; #3Mb/window
my $step="300000"; #300Kb/step
my $total="";
my $genome="";
my $intact="0"; #a switch to control input file. 0 for RepeatMasker .out file; 1 for LTR_retriever .pass.list file

#my ($totLTR, $iden, $iden_slope);

my $k=0;
foreach (@ARGV){
	$genome=$ARGV[$k+1] if /^-genome$/i;
	$total=$ARGV[$k+1] if /^-all$/i;
	$window=$ARGV[$k+1] if /^-window$/i;
	$step=$ARGV[$k+1] if /^-step$/i;
	$intact=1 if /^-intact$/i;
	$k++;
	}

$window=~s/,//g;
$step=~s/,//g;
die "$usage\n" unless -s "$total";

#extract LTR info from the RM.out file into bed format
open TOTAL, "awk '{if (\$6~/[0-9]/ && \$7-\$6+1>=80)print \$5\"\t\"\$6\"\\t\"\$7\"\\t\"\$11}' $total|sort -suV -k1,3|" or die "ERROR: $!\n$usage\n" if $intact == 0;
open TOTAL, "perl -nle 'my (\$id, \$superfam)=(split)[0,9]; my (\$chr, \$from, \$to)=(\$1, \$2, \$3) if \$id=~/^(.*)\:([0-9]+)\.\.([0-9]+)\$/; next unless defined \$chr; print \"\$chr\t\$from\t\$to\t\$superfam\"' $total|sort -suV -k1,3|" or die "ERROR: $!\n$usage\n" if $intact == 1;
open Genome, "<$genome" or die $!;

my $genome_len=0; #length of the genome
my %length; #store sequence length
my %total; #store all LTR sequence info
my %copia; #store copia LTR-RT info
my %gypsy; #store gypsy LTR-RT info
my %unknown; #store unknown LTR-RT info
my @seqID; #store chr ID names in input order
my $output=''; #stores output info

$/="\n>";
while (<Genome>){
	next if /^>\s?$/;
	chomp;
	s/>//g;
	s/^\s+//;
	my ($chr, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	$chr=~s/\s+$//; #remove space at the end of the seq ID
	push @seqID, $chr;
	$seq=length($seq);
	$genome_len+=$seq;
	$length{$chr}=$seq;
	$total{$chr}='0' x $seq; #create another '0' string that has the same length as the chr
	$copia{$chr}='0' x $seq; #create a '0' string that has the same length as the chr
	$gypsy{$chr}='0' x $seq; #create a '0' string that has the same length as the chr
	$unknown{$chr}='0' x $seq; #create a '0' string that has the same length as the chr
	}
$/="\n";
close Genome;

while (<TOTAL>){
	s/^\s+//;
	my ($chr, $from, $to, $superfam)=(split)[0,1,2,3];
	next unless exists $total{$chr};
	my $len=$to-$from+1;
	substr($total{$chr}, $from-1, $len)="a" x $len; #substitute '0' with 'a' where LTR sequence is occurred
	substr($copia{$chr}, $from-1, $len)="c" x $len if $superfam =~ /Copia/i; #substitute '0' with 'a' where LTR sequence is occurred
	substr($gypsy{$chr}, $from-1, $len)="g" x $len if $superfam =~ /Gypsy/i; #substitute '0' with 'a' where LTR sequence is occurred
	substr($unknown{$chr}, $from-1, $len)="u" x $len if $superfam =~ /unknown/i; #substitute '0' with 'a' where LTR sequence is occurred
	}
close TOTAL;

my ($tot_cop_count, $tot_gyp_count, $tot_ukn_count, $tot_all_count, $tot_cop_per, $tot_gyp_per, $tot_ukn_per, $tot_all_per)=(0, 0, 0, 0, 0, 0, 0, 0);
foreach my $chr (@seqID){
	my $all_count = $total{$chr} =~ tr/a/a/;
	my $copia_count = $copia{$chr} =~ tr/c/c/;
	my $gypsy_count = $gypsy{$chr} =~ tr/g/g/;
	my $unknown_count = $unknown{$chr} =~ tr/u/u/;
	$tot_cop_count += $copia_count;
	$tot_gyp_count += $gypsy_count;
	$tot_ukn_count += $unknown_count;
	$tot_all_count += $all_count;
	}

foreach my $chr (@seqID){
#estimate LAI based on windows and steps
	my $win_len = $window;
	for (my $start=1; $win_len == $window; $start += $step){ #update $win_len everytime, exit when != $window (occurs at chromosome end)
		my $end = $start+$window-1;
		$end = $length{$chr} if $end > $length{$chr};
		$win_len = $end-$start+1; #the actual size of the window (chromosome end may have win_len < window)
		my $win_seq_all = substr($total{$chr}, $start, $win_len); #total LTR sequence information in the win
		my $win_seq_cop = substr($copia{$chr}, $start, $win_len); #copia LTR-RT information in the window
		my $win_seq_gyp = substr($gypsy{$chr}, $start, $win_len); #gypsy LTR-RT information in the window
		my $win_seq_ukn = substr($unknown{$chr}, $start, $win_len); #unknown LTR-RT information in the window
		my $win_all_count = $win_seq_all =~ tr/a/a/; #total LTR sequence length
		my $win_cop_count = $win_seq_cop =~ tr/c/c/; #copia LTR-RT length
		my $win_gyp_count = $win_seq_gyp =~ tr/g/g/; #gypsy LTR-RT length
		my $win_ukn_count = $win_seq_ukn =~ tr/u/u/; #unknown LTR-RT length
		my $win_all_per = sprintf("%.4f", $win_all_count/$win_len); #propotion of total LTR sequence in the window
		my $win_cop_per = sprintf("%.4f", $win_cop_count/$win_len); #propotion of intact LTR-RT in the window
		my $win_gyp_per = sprintf("%.4f", $win_gyp_count/$win_len); #propotion of intact LTR-RT in the window
		my $win_ukn_per = sprintf("%.4f", $win_ukn_count/$win_len); #propotion of intact LTR-RT in the window
		$output .= "$chr\t$start\t$end\t$win_cop_per\t$win_gyp_per\t$win_ukn_per\t$win_all_per\n";
		}
	}
#estimate genome-wide LAI
$tot_all_per = sprintf("%.4f", $tot_all_count/$genome_len);
$tot_cop_per = sprintf("%.4f", $tot_cop_count/$genome_len);
$tot_gyp_per = sprintf("%.4f", $tot_gyp_count/$genome_len);
$tot_ukn_per = sprintf("%.4f", $tot_ukn_count/$genome_len);
print "Chr\tFrom\tTo\tCopia\tGypsy\tunknown\tTotal_LTR\n";
print "whole_genome\t1\t$genome_len\t$tot_cop_per\t$tot_gyp_per\t$tot_ukn_per\t$tot_all_per\n$output"; #print out all LTR info

