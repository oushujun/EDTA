#!/usr/bin/env perl
use warnings;
use strict;
# Shujun Ou (shujun.ou.1@gmail.com)
# v1	05/03/2023

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
open Dirt_clas, ">$seq.dirt.list" or die $!;

# translate terms
# Translaste specialized/alia -> general TE class
my %lib = ('DNA', 'DNA', 'MITE', 'DNA', 'Helitron', 'DNA', 'RC', 'DNA', 'Maverick', 'DNA', 
	'pararetrovirus', 'DNA', 'TIR', 'DNA', 'MITE', 'DNA', 'Retrovirus', 'Gypsy',
	'LTR', 'LTR', 'DIRS', 'LTR', 'Copia', 'Copia', 'Gypsy', 'Gypsy', 'LINE', 'LINE', 
	'Penelope', 'LINE', 'LINE-1', 'LINE','SINE', 'SINE', 'SINE?', 'SINE', 'SINE1', 'SINE', 
	'SINE2', 'SINE', 'SINE3', 'SINE', 'EnSpm_CACTA', 'DTC', 'MuDR_Mutator', 'DTM', 
	'PIF_Harbinger', 'DTH', 'Tc1_Mariner', 'DTT', 'hAT', 'DTA', 'mixture', 'mixture', 
	'Unassigned', 'unknown', 'unknown', 'unknown', 'Unknown', 'unknown', 'repeat', 'unknown',
	'Satellite', 'Satellite', 'snRNA', 'snRNA');

# read sequence
$/ = "\n>";
my %seq;
my %extra; # store sequence ID annotation
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

# filter out mis-classified based on TEsorter results
while (<Clas>){
	chomp;
	(print Cln_clas "$_\n" and next) if /^#/;
	my ($id, $order, $superfam) = (split)[0,1,2];
	my ($ori_order, $ori_superfam) = ('unknown', 'unknown');
	($ori_order, $ori_superfam) = ($1, $2) if $id =~ /#(.*)\/(.*)$/;
	$ori_order = $1 if $id =~ /#(.*)$/ and $ori_order eq 'unknown';
	(delete $seq{$id} and print Dirt_clas "$_\n" and next) if $lib{$ori_order} ne $lib{$order} and $lib{$ori_order} ne 'unknown';
	(delete $seq{$id} and print Dirt_clas "$_\n" and next) if $ori_superfam ne 'Helitron' and $order eq 'Helitron';
	(delete $seq{$id} and print Dirt_clas "$_\n" and next) if $ori_superfam eq 'Helitron' and $order ne 'Helitron';
	(delete $seq{$id} and print Dirt_clas "$_\n" and next) if $ori_superfam ne 'pararetrovirus' and $order eq 'pararetrovirus';
	(delete $seq{$id} and print Dirt_clas "$_\n" and next) if $ori_superfam ne 'Maverick' and $order eq 'Maverick';
	my $extra = 'NA';
	$extra = $extra{$id} if defined $extra{$id};
	if ($extra eq 'NA'){
		$extra = '';
		} else {
		$extra = "\t$extra";
		}
	if (defined $seq{$id}){
		print Cln_seq ">${id}$extra\n$seq{$id}\n" and delete $seq{$id};
		print Cln_clas "$_\n";
		}
	}

# output non-classified
foreach my $id (keys %seq){
	my $extra = 'NA';
	$extra = $extra{$id} if defined $extra{$id};
	if ($extra eq 'NA'){
		$extra = '';
		} else {
		$extra = "\t$extra";
		}
	print Cln_seq ">${id}$extra\n$seq{$id}\n";
	my ($ori_order, $ori_superfam) = ("Unknown", "Unknown");
	if ($id =~ /#(.*)\/(.*)$/){
		($ori_order, $ori_superfam) = ($1, $2);
		} elsif ($id =~ /#(.*)$/){
		($ori_order, $ori_superfam) = ($1, "Unknown");
		}
	$ori_order = "TIR" if $ori_superfam =~ /^DT/;
	($ori_order, $ori_superfam) = ("Helitron", "unknown") if $ori_superfam eq 'Helitron';
	print Cln_clas "$id\t$ori_order\t$ori_superfam\tunknown\tnone\t?\tnone\n";
	}

close Clas;
close Cln_clas;
close Cln_seq;
close Dirt_clas;
