#!/usr/bin/perl -w
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
#print "$chr, $element_start, $element_end, $TE_class, $strand, $info\n";
#	$TE_ID = "$chr:$element_start..$element_end";
	$element_length = $element_end - $element_start + 1;
	print GFF "$chr\t$annotator\t$TE_class\t$element_start\t$element_end\t.\t$strand\t.\t$info\n";
	}
close BED;
close GFF;
	

#B73V4_ctg1      4992    9052    LTR     LTR/Copia       +       chr9-D-120499978;Method=homology
#B73V4_ctg1      9053    9842    TIR     DNA/CACTA       +       ZM00044_consensus;Method=homology
#B73V4_ctg1      9458    10030   TIR     DNA/DTC +       TE_00010295;Method=homology
#B73V4_ctg1      10170   10311   TIR     DNA/CACTA       +       ZM00047_consensus;Method=homology
#B73V4_ctg1      10312   10745   TIR     DNA/CACTA       -       ZM00044_consensus;Method=homology
#B73V4_ctg1      10929   11193   LTR     LTR/Gypsy       +       xilon-diguus_AC203313-7774;Method=homology
#B73V4_ctg1      11237   12375   LTR     LTR/Gypsy       +       prem1_AC206253-9147;Method=homology
#B73V4_ctg1      12816   14905   LTR     LTR/Gypsy       +       prem1_AC206253-9147;Method=homology
#B73V4_ctg1      14878   15656   LTR     LTR/Gypsy       +       TE_00005798_INT;Method=homology
#B73V4_ctg1      15630   16189   LTR     LTR/Gypsy       +       TE_00005798_INT;Method=homology
#B73V4_ctg1      16194   16315   TIR     DNA/CACTA       +       ZM00101_consensus;Method=homology
