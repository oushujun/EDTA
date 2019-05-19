#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com; 05/16/2019)

my $usage = "
	Rename TIR-Learner output fasta to RepeatMasker readible format\n
	perl rename_tirlearner.pl TIR-Learner.out.fasta > TIR-Learner.out.fasta.renamed\n\n";

my $len_cutoff = 600; #TIR elements <= this value are classified as MITEs, others DNA TEs.

my $fasta = $ARGV[0];
die "\nThe input file is not found!\n$usage" unless -s $fasta;


open Fasta, "<$fasta" or die $usage;
$/ = "\n>";
while (<Fasta>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my $class = "NA";
	if (length $seq > $len_cutoff){
		$class = "DNA";
		} else {
		$class = "MITE";
		}
	my ($name, $superfam, $tsd) = ('NA', 'NA', 'NA');
	($name, $superfam, $tsd) = ($1, $2, $3) if $id =~ /^(.*)_(D.*)_TIR.*TSD:([ATCGNX_]+)_[0-9]+.*/i;
	print ">$name#$class/$superfam TSD:$tsd\n$seq\n";
	}
close Fasta;

