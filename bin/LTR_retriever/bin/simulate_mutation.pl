#!/usr/bin/env perl -w
use strict;

my $mut_rate="0.02"; #%2 (default) of the genome will be mutated
my @chars=("a", "t", "c", "g"); #mutated bp will be marked in lowercase
my $genome;

my $k=0;
foreach (@ARGV){
	$mut_rate=$ARGV[$k+1] if /^-u$/i;
	$genome=$ARGV[$k+1] if /^-g$/i;
	$k++;
	}
open Genome, "<$genome" or die $!;

$/="\n>";
while (<Genome>){
	next if /^>\s?$/;
	chomp;
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$id=~s/\s+$//g;
	$seq=~s/\s+//g;
	my $seq_len=length $seq;
	my $mut_num=int($seq_len*$mut_rate);
	foreach (1..$mut_num){
		my $mut_pos=int(rand($seq_len));
		substr($seq,$mut_pos,1)=$chars[rand @chars]
		}
	print ">$id\n$seq\n";
	}
$/="\n";
close Genome;
