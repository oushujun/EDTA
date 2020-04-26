#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com; 05/16/2019)

my $usage = "
	Rename TIR-Learner output fasta to RepeatMasker readible format\n
	perl rename_tirlearner.pl TIR-Learner.out.fasta > TIR-Learner.out.fasta.renamed\n\n";

my $len_cutoff = 600; #TIR elements <= this value are classified as MITEs, others DNA TEs.

#my $fasta = $ARGV[0];
#die "\nThe input file is not found!\n$usage" unless -s $fasta;


#open Fasta, "<$fasta" or die $usage;
$/ = "\n>";
#while (<Fasta>){
while (<>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my $class = "NA";
	if (length $seq > $len_cutoff){
		$class = "DNA";
		} else {
		$class = "MITE";
		}
	my ($name, $superfam, $tsd) = ('NA', 'TIR', 'NA');
	($name, $superfam, $tsd) = ($1, $2, $3) if $id =~ /^(.*)_(D.*)_TIR.*TSD:([ATCGNX_]+)_[0-9]+.*/i; #renamed TIR-Learner output
	#e.g.: TIR-Learner_11_8901442_8902889_DTM_TIR:GGAAAAAGTA_TACTTTTTCC_100.0_TSD:TTTACTTTT_TCTACTTTT_88.89-+-1448
	($name, $superfam, $tsd) = ($1, $2, $3) if $id =~ /^(.*):[0-9mM]+:.*_(D[A-Z]+).*TSD:([ATCGNX_]+)_[0-9]+.*/i;
	#e.g.: 10:16642857:16643502:8m1M4m1M3m1M:CTCTCATG_DTA_200-+-TIR:TAGGGGTGAA_TCCACCCCTA_90.0_TSD:CTCTCATG_CTCTCATG_100.0
	if ($name ne 'NA'){
		if ($name =~ /^(.*)_([0-9]+)_([0-9]+)/){
			my ($chr, $start, $end) = ($1, $2, $3);
			$name = "$chr:$start..$end" if defined $chr;
			}
		print ">$name#$class/$superfam TSD:$tsd\n$seq\n";
		} else {
		($name, $superfam) = ($1, $2) if $id =~ /^(.*)#.*\/(D.*)/;
		$name = $id unless $name ne 'NA';
		print ">$name#$class/$superfam\n$seq\n";
		}
	}
#close Fasta;

