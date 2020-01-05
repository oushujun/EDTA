#!/usr/bin/perl -w
use strict;

# Convert gff files to enriched bed format
# Shujun Ou (shujun.ou.1@gmail.com)
# 12/19/2019

my $usage = "Usage: 
	perl gff2bed.pl file.gff [structural|homology] > file.bed
	cat *.gff | perl gff2bed.pl - [structural|homology] > file.bed \n";
die "Please indicate if the specified GFF file is generated based on structural features or homology!\n$usage\n" unless defined $ARGV[1] and $ARGV[1] =~ /^structural$|^homology$/i;
my $method = $ARGV[1];
open GFF, "<$ARGV[0]" or die $usage;

while (<GFF>){
	chmod;
	next if /^#/;
	my ($id, $type, $start, $end, $dir, $info) = (split)[0,2,3,4,6,8];
	my $class = "LTR";
	$class = "TIR" if $type =~ /DT/ or ($type =~ /DNA|MITE/ and $type !~ /Helitron|DHH/);
	$class = "Cent" if $type =~ /Cent/i;
	$class = "knob" if $type =~ /knob/i;
	$class = "LINE" if $type =~ /LINE|RIL/i;
	$class = "SINE" if $type =~ /SINE|RIS/i;
	$class = "rDNA" if $type =~ /rDNA/i;
	$class = "SAT" if $type =~ /SAT/i;
	$class = "subtelomere" if $type =~ /subtelomere/i;
	$class = "Helitron" if $type =~ /Helitron|DHH/i;
	$class = $1 if $type =~ /^(.*)\/.*/ and $1 !~ /DNA|MITE/i;
	$info = "ID=$info" if $method eq "homology";
	$info =~ s/;Method=.*//; #remove method info at the end of the annotation
	next if $type eq "repeat_region";
	next unless defined $id and defined $start;
	print "$id\t$start\t$end\t$class\t$type\t$dir\t$info;Method=$method\n";
	}

