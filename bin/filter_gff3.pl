#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com)
#01/04/2020
#05/21/2020

my $usage = "\n\tperl filter_gff3.pl file.gff3 file.list > new.gff3\n\n";

my $gff = $ARGV[0];
my $list = $ARGV[1];

open GFF, "<$gff" or die $usage;
open List, "<$list" or die $usage;
open RM, ">$gff.removed" or die $usage;

my %parent;
my %name;
my %ID;
while (<List>){
	my ($tag, $id) = (split)[0,1];
	$parent{$id} = $id if lc $tag eq 'parent';
	$name{$id} = $id if lc $tag eq 'name';
	$ID{$id} = $id if lc $tag eq 'id';
	$id =~ s/#.*//; #remove classification info and store id again
	$parent{$id} = $id if lc $tag eq 'parent';
	$name{$id} = $id if lc $tag eq 'name';
	$ID{$id} = $id if lc $tag eq 'id';
	}
close List;

while (<GFF>){
	chomp;
	print "$_\n" if /^#/;
	(print RM "$_\n" and next) if /^#/;
	my $info = (split)[8];
	my ($parent, $name, $ID) = ('NA', 'NA', 'NA');
	$parent = $1 if $info =~ /Parent=(.*?);/i;
	$name = $1 if $info =~ /Name=(.*?);/i;
	$ID = $1 if $info =~ /ID=(.*?);/i;
	if (defined $parent and exists $parent{$parent}){
		print RM "$_\n";
		next;
		}
	if (defined $name and exists $name{$name}){
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

