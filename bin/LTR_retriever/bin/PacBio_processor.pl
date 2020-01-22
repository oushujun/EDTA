#covert PacBio fastq file into fasta file with simple filtering options
#Usage: perl PacBio_processor.pl PacBio.fastq > PacBio.fasta
#06/30/2016 Shujun Ou (oushujun@msu.edu)


#!/usr/bin/env perl -w
use strict;

my $minLength=500;
my $maxLength=50000;
my $minRQ=0.8;

open FASTQ, "<$ARGV[0]" or die "ERROR: $!";

while (<FASTQ>){
	chomp;
	my $name=$_;
	next unless $name=~/^@/;#first line of the FASTA file should start with "@", otherwise next
	my $seq=<FASTQ>;
	my $strand=<FASTQ>;
	my $qual=<FASTQ>;
	next unless $seq=~/^[ATCG]+/i ;#next line should be sequence and should not contain ">", otherwise next
	next unless defined $seq;
	$seq=~s/\s+//g;
	my $length=length $seq;
	my $RQ=0;
	($name, $RQ)=($1, $2) if $name=~/^@.*\/.*\/(.*)\s+RQ=([0-9.]+)/;
	print ">read$name\n$seq\n" if ($RQ>=$minRQ and $length>=$minLength and $length<=$maxLength);
	}

