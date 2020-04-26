#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "\n\tperl bed2gff.pl EDTA.TE.combo.bed\n\n";
die $usage unless -s $ARGV[0];
my $date=`date -u`;
chomp ($date);
open BED, "sort -sV -k1,1 $ARGV[0] |" or die "ERROR: $!";
open GFF, ">$ARGV[0].gff" or die "ERROR: $!";
print GFF "##gff-version 2\n##date $date
##Chromosome Annotator Repeat_class/superfamily Start End Score Strand Phase Annotation\n";

my $annotator="EDTA";
while (<BED>){
	chomp;
	my ($chr, $element_start, $element_end, $element_length, $strand, $TE_ID, $TE_class, $info);
	($chr, $element_start, $element_end, $TE_class, $strand, $info) = (split)[0,1,2,4,5,6];
	$element_length = $element_end - $element_start + 1;
	print GFF "$chr\t$annotator\t$TE_class\t$element_start\t$element_end\t.\t$strand\t.\t$info\n";
	}
close BED;
close GFF;

