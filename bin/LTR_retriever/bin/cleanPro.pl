#Usage: perl cleanPro.pl fastafile blastxout(fmt6) -o [1|0] -l [int] -c [0,1] > sequence_name_out
#Optional parameters:
#-o [1|0] output option, default 1, to output the entries with protein hits, 0 for protein-free hits
#-l [int] #bp, cumulated alignment length more than this will be treated as protein contained
#-c [0,1] #propotional, cumulated alignment coverage more than this will be treated as protein contained

#Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University.
#Version: v1.0 12/17/2014


#!usr/bin/perl -w
use strict;

##filters for the sequence alignment
my $len_cut=1000; #bp, cumulated alignment length more than this will be treated as protein contained
my $cvg_cut=0.3; #propotional, cumulated alignment coverage more than this will be treated as protein contained

##filters for the hits
my $evalue=0.001; #0.001; #hits with evalue more than this number will not be counted.
my $identity=30; #%, hits with identity less than this number will not be counted.
my $length=90; #bp, hits with alignment length less than this number will not be counted.

##define output format
my $out=1; #1 to output the entries with protein hits, 0 for protein-free hits

##parameter customization
my $i=0;
foreach (@ARGV){
	$len_cut=$ARGV[$i+1] if /^-l|length|-prolen$/i;
	$cvg_cut=$ARGV[$i+1] if /^-c|coverage|-procvg$/i;
	$length=$ARGV[$i+1] if /^-prolensig$/i;
	$out=$ARGV[$i+1] if (/^-o|prolist$/i);
	$i++;
	}

open FA, "<$ARGV[0]" or die "ERROR: $!";
$/="\n>";
my %fasta;
my %blast;
while (<FA>){
	s/>//;
	chomp;
	s/^\s+//;
	my ($head, $seq)=(split /\n/, $_, 2);
	$head=(split)[0];
	$seq=~s/\s//g;
	my $length=length $seq;
	$fasta{$head}=[$seq, $length];
	$blast{$head}=();
	}
$/="\n";
close FA;

open BLAST, "<$ARGV[1]" or die "ERROR: $!";
while (<BLAST>){
	next if /^#/;
	next unless /[0-9]+/;
	s/^\s+//;
#query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	my ($qid, $sid, $idty, $len, $mis, $gap, $qs, $qe, $ss, $se, $eval)=split;
	next if ($idty<$identity or $eval>$evalue or $len<$length);
	($qs, $qe)=($qe, $qs) if $qe<$qs;
	foreach my $i ($qs..$qe){
		$blast{$qid}->[$i]=1;
		}
	}

foreach my $qid (keys %blast){
	my $count=0;
	foreach my $i (0..$#{$blast{$qid}}){
		$count++ if defined ${$blast{$qid}}[$i];
		}
	next unless defined $fasta{$qid}->[1];
	my $cvg=$count/$fasta{$qid}->[1];
	print "$qid\t$fasta{$qid}->[1]\t$count\t$cvg\n" if (($count>=$len_cut or $cvg>=$cvg_cut) and $out==1);
	print "$qid\t$fasta{$qid}->[1]\t$count\t$cvg\n" if (($count<$len_cut and $cvg<$cvg_cut) and $out==0);
	}
