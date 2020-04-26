#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com)
#01/04/2020

my $usage = "\n\tperl filter_gff.pl file.gff file.list > new.gff\n\n";

my $gff = $ARGV[0];
my $list = $ARGV[1];

open GFF, "<$gff" or die $usage;
open List, "<$list" or die $usage;
open RM, ">$gff.removed" or die $usage;

my %parent;
my %ID;
while (<List>){
	my ($tag, $id) = (split)[0,1];
	$parent{$id} = $id if lc $tag eq 'parent';
	$ID{$id} = $id if lc $tag eq 'id';
	}
close List;

while (<GFF>){
	chomp;
	print "$_\n" if /^#/;
	(print RM "$_\n" and next) if /^#/;
	my $info = (split)[8];
	my $parent = $1 if $info =~ /Parent=(.*?);/;
	my $ID = $1 if $info =~ /ID=(.*?);/;
	if (defined $parent and exists $parent{$parent}){
		print RM "$_\n";
		next;
		}
	if (defined $ID and exists $ID{$ID}){
		print RM "$_\n";
		next;
		}
	print "$_\n";
	}
close GFF;
close RM;

