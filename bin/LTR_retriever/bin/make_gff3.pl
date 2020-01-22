##Generate gff3 file from pass list of LTR_retriever
##Usage: perl make_gff3.pl genome.fa LTR.pass.list
##Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University
##Version: 
#	1.0 04-09-2015
#	2.0 04-11-2018 Add ltr_digest support



#!/usr/bin/env perl -w
use strict;

open FA, "<$ARGV[0]" or die "ERROR: $!";
open List, "<$ARGV[1]" or die "ERROR: $!";
open GFF, ">$ARGV[1].gff3" or die "ERROR: $!";
print GFF "##gff-version   3\n";
$/="\n>";
my $chr_info='';
my %seq_flag;
my $j=0;
while (<FA>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	my $len=length ($seq);
	print GFF "##sequence-region   $id 1 $len\n";
	$chr_info.="#$id\n";
	$seq_flag{$id}=$j;
	$j++;
	}
print GFF "$chr_info";
$/="\n";

my $annotator="LTR_retriever";
my $i=1;
while (<List>){
	next if /^#/;
	next if /^\s+/;
	my ($chr, $element_start, $element_end, $element_length, $lLTR_start, $lLTR_end, $lLTR_length, $rLTR_start, $rLTR_end, $rLTR_length, $seq_ID, $loc, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam);
	($loc, undef, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam) = (split);
#chr03:19323647..19328532        pass    motif:TGCA      TSD:GTCGC       19323642..19323646      19328533..19328537      IN:19323854..19328325   99.03
	($chr, $lLTR_start, $rLTR_end)=($1, $2, $3) if $loc=~/^(.*):([0-9]+)..([0-9]+)$/;
	($lLTR_start, $rLTR_end)=($rLTR_end, $lLTR_start) if $rLTR_end<$lLTR_start;
	($lLTR_end, $rLTR_start)=($1-1, $2+1) if $IN=~/^IN:([0-9]+)..([0-9]+)$/;
	$motif=~s/motif://gi;
	$TSD=~s/TSD://gi;
	$lTSD=~s/\.\./\t/;
	$rTSD=~s/\.\./\t/;
	$element_start=(split /\s+/, $lTSD)[0];
	$element_end=(split /\s+/, $rTSD)[1];
	my $id = "$chr:$lLTR_start..$rLTR_end#LTR/$supfam";
	if ($TSD eq "NA"){
		$element_start=$lLTR_start;
		$element_end=$rLTR_end;
		}
	my $chr_ori=$chr;
	print GFF "$chr\t$annotator\trepeat_region\t$element_start\t$element_end\t.\t$strand\t.\tID=repeat_region$i\n";
	print GFF "$chr\t$annotator\ttarget_site_duplication\t$lTSD\t.\t$strand\t.\tParent=repeat_region$i\n" unless $TSD eq "NA";
	print GFF "$chr\t$annotator\tLTR/$supfam\t$lLTR_start\t$rLTR_end\t.\t$strand\t.\tID=$id;Parent=repeat_region$i;motif=$motif;tsd=$TSD;ltr_identity=$sim;seq_number=$seq_flag{$chr_ori}\n";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$lLTR_start\t$lLTR_end\t.\t$strand\t.\tParent=$id\n";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$rLTR_start\t$rLTR_end\t.\t$strand\t.\tParent=$id\n";
	print GFF "$chr\t$annotator\ttarget_site_duplication\t$rTSD\t.\t$strand\t.\tParent=repeat_region$i\n" unless $TSD eq "NA";
	print GFF "###\n";
	$i++;
	}
