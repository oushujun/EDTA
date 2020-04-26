#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com) 11/15/2019

my $usage = "\n\tperl make_gff_with_intact.pl EDTA.intact.fa\n\n";
die $usage unless -s $ARGV[0];
my $date=`date -u`;
chomp ($date);
open ID, "grep \\> $ARGV[0] | sort -sV -k1,1 |" or die "ERROR: $!";
open GFF, ">$ARGV[0].gff" or die "ERROR: $!";
print GFF "##gff-version 2\n##date $date
##Chromosome Annotator Repeat_class/superfamily Start End Score Strand Phase Repeat_famliy\n";

my $annotator="EDTA";
while (<ID>){
	s/\s+//g;
	s/>//;
	my ($chr, $element_start, $element_end, $element_length, $strand, $TE_ID, $TE_class);
	($chr, $element_start, $element_end, $TE_class) = ($1, $2, $3, $4) if /^(.*):([0-9]+)\.\.([0-9]+)#(.*)$/;
	$TE_ID = "$chr:$element_start..$element_end";
	$element_length = $element_end - $element_start + 1;
	next if $element_length < 80;
	$strand = ".";
        print GFF "$chr\t$annotator\t$TE_class\t$element_start\t$element_end\t.\t$strand\t.\tID=$TE_ID#$TE_class\n";
        }
close ID;
close GFF;

