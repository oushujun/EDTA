#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "\nFilter sequence based on TEsorter classifications. Unclassified sequences will also be output to the clean file.
	Usage: perl cleanup_misclas.pl sequence.fa.rexdb.cls.tsv
	Author: Shujun Ou (shujun.ou.1\@gmail.com) 10/11/2019
	\n";

my $clas = $ARGV[0];
my $seq = $ARGV[1];
$seq = $1 if $clas =~ /^(.*)\.[a-z]+db.cls.tsv$/ and !defined $seq;
die $usage unless -s $clas;
die "The sequence file for $clas is not found in the same folder!\n
	You may specify the sequence file with perl cleanup_misclas.pl TEsorter.cls.tsv sequence.fa\n" unless -s $seq;

open Seq, "<$seq" or die $!;
open Clas, "<$clas" or die $!;
open Cln_clas, ">$seq.cln.list" or die $!;
open Cln_seq, ">$seq.cln" or die $!;

# translate terms
my %lib = ('DNA', 'DNA', 'MITE', 'DNA', 'Helitron', 'DNA', 'Maverick', 'DNA', 'pararetrovirus', 'DNA', 'TIR', 'DNA', 'MITE', 'DNA', 'LTR', 'LTR', 'Copia', 'Copia', 'Gypsy', 'Gypsy', 'LINE', 'LINE', 'SINE', 'SINE', 'EnSpm_CACTA', 'DTC', 'MuDR_Mutator', 'DTM', 'PIF_Harbinger', 'DTH', 'Tc1_Mariner', 'DTT', 'hAT', 'DTA', 'mixture', 'mixture', 'unknown', 'unknown');

# read sequence
$/ = "\n>";
my %seq;
my %extra;
while (<Seq>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	my $extra;
	($id, $extra) = (split /\s+/, $id, 2);
	$extra = "NA" unless defined $extra;
	$seq =~ s/\s+//g;
	$seq{$id} = $seq;
	$extra{$id} = $extra;
	}
close Seq;
$/ = "\n";

# filter out mis-classified
while (<Clas>){
	chomp;
	(print Cln_clas "$_\n" and next) if /^#/;
	my ($id, $order, $superfam) = (split)[0,1,2];
	my ($ori_cla, $ori_superfam) = ($1, $2) if $id =~ /#(.*)\/(.*)$/;
	delete $seq{$id} and next if $ori_superfam ne $lib{$superfam} and $superfam ne 'unknown';
	delete $seq{$id} and next if $lib{$ori_cla} ne $lib{$order} and $superfam eq 'unknown';
	delete $seq{$id} and next if $ori_superfam ne 'Helitron' and $order eq 'Helitron';
	delete $seq{$id} and next if $ori_superfam ne 'pararetrovirus' and $order eq 'pararetrovirus';
	delete $seq{$id} and next if $ori_superfam ne 'Maverick' and $order eq 'Maverick';
	my $extra = 'NA';
	$extra = $extra{$id} if defined $extra{$id};
	if ($extra eq 'NA'){
		$extra = '';
		} else {
		$extra = "\t$extra";
		}
	print Cln_seq ">${id}$extra\n$seq{$id}\n" and delete $seq{$id};
	print Cln_clas "$_\n";
	}

# output non-classified
foreach my $id (keys %seq){
	print Cln_seq ">$id\n$seq{$id}\n";
	my ($ori_cla, $ori_superfam) = ($1, $2) if $id =~ /#(.*)\/(.*)$/;
	$ori_cla = "TIR" if $ori_superfam =~ /^DT/;
	($ori_cla, $ori_superfam) = ("Helitron", "unknown") if $ori_superfam eq 'Helitron';
	print Cln_clas "$id\t$ori_cla\t$ori_superfam\tunknown\tnone\t?\tnone\n";
	}

close Clas;
close Cln_clas;
close Cln_seq;

