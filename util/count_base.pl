#!/usr/bin/env perl
use warnings;
use strict;

my $seperate=0;
if (defined $ARGV[1]){
	$seperate=1 if $ARGV[1] eq '-s'; #use -s to output each chromosome length
}

my $length=0; #length of all seq including gaps
my $N_length=0; #length of Ns
my $seq_count=0; #number of seqs

open File, "<$ARGV[0]" or die $!;
$/="\n>";
while (<File>){
	chomp;
	s/>//g;
	my ($id, $seq)=(split /\n/, $_, 2);
	$id = (split)[0];
	$seq=~s/\s+//g;
	my $N_count=$seq=~tr/Nn//;
	my $chr_len=length $seq;
	$length+=$chr_len;
	next unless $chr_len>0;
	$seq_count++;
	my $chr_mis=$N_count/$chr_len;
	print "$id\t$chr_len\t$N_count\t$chr_mis\n" if $seperate==1;
	$N_length+=$N_count;
	}
my $missing=$N_length/$length;
print "All\t$length\t$N_length\t$missing\t$seq_count\n";
close File;
