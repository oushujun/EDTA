#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com) 07/02/2020

my $usage = "\n\tperl make_bed_with_intact.pl EDTA.intact.fa > EDTA.intact.bed \n\n";
die $usage unless -s $ARGV[0];

open ID, "grep \\> $ARGV[0] | sort -sV |" or die "ERROR: $!";

my $method = "structural";
while (<ID>){
	s/>//;
	my ($chr, $element_start, $element_end, $element_length, $TE_ID, $TE_class, $iden, $score, $strand, $phase, $extra, $id, $ext);
	($id, $ext) = (split /\s+/, $_, 2);
	($chr, $element_start, $element_end, $TE_class) = ($1, $2, $3, $4) if $id =~ /^(.*):([0-9]+)\.\.([0-9]+)#(\S+)/;
	$TE_ID = "$chr:$element_start..$element_end";
	($element_end, $element_start) = ($element_start, $element_end) if $element_start > $element_end;
	$element_length = $element_end - $element_start + 1;
	next if $element_length < 80;
	($iden, $score, $strand, $phase, $extra) = ("NA", ".", ".", ".", ".");
	if (defined $ext){
		$strand = $1 if $ext =~ /^([\-|+])\|/;
		}
        print "$chr\t$element_start\t$element_end\t$TE_ID\t$TE_class\t$method\t$iden\t$score\t$strand\t$phase\t$extra\n";
        }
close ID;

