#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

##Generate gff3 file from RepeatMasker .out file of LTR_retriever
##Usage: perl make_gff3.pl genome.fa.out
##Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University
##Version: 	3.0 07-02-2020
#		2.0 05-19-2020
#		1.0 01-27-2018


my $usage = "\n\tperl make_gff3_with_RMout.pl RM.out\n\n";
die $usage unless -s $ARGV[0];
my $annotator="LTR_retriever";

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

#read SO info and defined sequence ontology
my $SO = "$script_path/TE_Sequence_Ontology.txt";
open SO, "<$SO" or die "The sequence ontology file 'TE_Sequence_Ontology.txt' is not found in $script_path!\n";
my (%class, %SO);
while (<SO>){
	next if /#/;
	next if /^(\s+)?$/;
	my ($so_name, $so_id, $so_alias) = (split /\s+/, $_, 3);
	$so_alias =~ s/\s+//;
	$SO{$so_name} = $so_id;
	foreach my $alia ((split /,/, $so_alias)){
		$class{$alia} = $so_name;
		}
	}
close SO;

my $date=`date -u`;
chomp ($date);
open RMout, "sort -sV -k5,5 $ARGV[0]|" or die "ERROR: $!";
open GFF, ">$ARGV[0].gff3" or die "ERROR: $!";
print GFF "##gff-version 3\n##date $date
##seqid source repeat_class/superfamily start end sw_score strand phase attributes\n";

my %seq_flag;
my $i = 0; #annotation count
while (<RMout>){
	s/^\s+//;
	my ($SW_score, $div, $iden, $chr, $chr_len, $element_start, $element_end, $element_length, $left_len, $strand, $TE_ID, $TE_class);
	($SW_score, $div, $chr, $element_start, $element_end, $left_len, $strand, $TE_ID, $TE_class)=(split)[0,1,4,5,6,7,8,9,10];
	next unless defined $SW_score and $SW_score =~ /[0-9]+/;
	my $so = $TE_class;
	if (exists $class{$TE_class}){
		$so = $class{$TE_class}
		} else {
		print "$TE_class not found in the TE_SO database, will use the general term 'repeat_region\tSO:0000657' to replace it.\n";
		$so = "repeat_region";
		}

	$TE_ID=~s/_INT-int$/_INT/;
	$element_length=$element_end-$element_start+1;
	next if $SW_score<300 and $element_length<80;
	$strand="-" if $strand eq "C";
	unless (exists $seq_flag{$chr}){
		$seq_flag{$chr}=$chr;
		$left_len=~s/[\(\) ]+//g;
		$left_len=0 unless $left_len=~/^[0-9]+$/;
		$chr_len=$element_end+$left_len;
		print GFF "##sequence-region $chr 1 $chr_len\n";
		}
	$iden = 1 - $div/100;
	my $info = "ID=TE_annot_$i;Name=$TE_ID#$TE_class;Classification=$TE_class;Sequence_ontology=$SO{$so};Identity=$iden;Method=homology";
	print GFF "$chr\t$annotator\t$so\t$element_start\t$element_end\t$SW_score\t$strand\t.\t$info\n";
	$i++; #annotation count
	}
close RMout;
close GFF;
