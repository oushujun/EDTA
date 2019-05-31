#Shujun Ou
#usage:  perl make_lib.pl *.scn.list *.ltrTE.lib *.ltrTE.lib.clust.info 



#!usr/bin/perl -w
use strict;

open List, "<$ARGV[0]" or die "ERROR: $!";
open Seq, "<$ARGV[1]" or die "ERROR: $!";
open Clust, "<$ARGV[2]" or die "ERROR: $!";
open Out, ">$ARGV[1].clust" or die "ERROR: $!";

my %list;
my $j=0;
while (<List>){
#Chr1:244893..246841[1]  Chr1:244893..245187     91.86
	next if /^\s+$/;
	s/^\s+//;
	my ($id, undef, $sim)=split;
	$list{"tmpseq_$j"}=[$id, $sim];
	$j++;
	}

my %bank;
my $i=0;
$/="\n>";
while (<Seq>){
	s/>//g;
	next if /^\s+$/;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	next if $seq eq '';
	$seq=~s/\s+//g;
	$bank{"tmpseq_$i"}=">$id\n$seq";
	$i++;
	}
$/="\n";

while (<Clust>){
	s/^\s+//;
	my @clust=(split /\s+/, $_);
	my $sim=0;
	my $rep='';
	foreach my $tag (@clust){
		if ($list{$tag}->[1]>$sim){
			$sim=$list{$tag}->[1];
			$rep=$tag;
			}
		}
	print Out "$bank{$rep}\n" if $rep ne '';
	}

