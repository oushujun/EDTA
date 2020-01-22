#To run MGEScan-ltr on genomes containing more than one sequence/scaffolds
#Usage: perl run_MGEScan.pl genome
#Author: Shujun Ou (oushujun@msu.edu)	08/14/2016



#!/usr/bin/env perl -w
use strict;

my $genome=$ARGV[0];
my $DAWGPAWS_path="/mnt/home/oushujun/git_bin/MGEScan_LTR"; #path to the DAWGPAWS program find_ltr_DAWGPAWS.pl
open Genome, "<$genome" or die "ERROR: $!";
die "ERROR: $!" unless (-e "$DAWGPAWS_path/find_ltr_DAWGPAWS.pl");

my $date=`date +"%m-%d-%y_%H%M"`;
chomp ($date);
`mv $genome.ltrpos $genome.pre$date.ltrpos` if -s "$genome.ltrpos";

$/="\n>";
while (<Genome>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_,  2);
	$seq=~s/\s+//g;
	$id=~s/\s+//g;
	open Seq, ">$id.temp" or die "ERROR: $!";
	print Seq ">$id\n$seq\n";
	close Seq;
	`perl $DAWGPAWS_path/find_ltr_DAWGPAWS.pl -seq=$id.temp`;
	`rm $id.temp`;
	`awk '{print "$id\t"\$0}' $id.temp.ltrpos >> $genome.ltrpos`;
	}
close Genome;
