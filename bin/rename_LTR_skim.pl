#!/usr/bin/env perl
use warnings;
use strict;

#Usage: this script is developed to output gff3 format of intact LTRs and add classification info to LTR sequences that only have bare msu locus names
#perl rename_LTR.pl genome.fa target_sequence.fa LTR_retriever.defalse
#Shujun Ou (05/21/2020)

my $usage = "\n\tperl rename_LTR_skim.pl target_sequence.fa LTR_retriever.defalse\n\n";

my $seq = $ARGV[0]; #target seq
my $anno = $ARGV[1]; #annotation of the target seq

open Seq, "<$seq" or die $usage;
open Anno, "<$anno" or die $usage;

# read annotation info
my %anno;
while (<Anno>){
	next unless /motif/;
	my $nextline = <Anno>;
	my ($loc, $lLTR_length, $rLTR_length);
	my $supfam = 'unknown';
	($loc, $supfam) = (split)[0,9];
#	chr03:19323647..19328532        pass    motif:TGCA      TSD:GTCGC       19323642..19323646      19328533..19328537      IN:19323854..19328325   99.03
#	B73V4_ctg1:5489..14866  false   motif:CGGC      TSD:CTAT        5485..5488      14867..14870    IN:6808..13534  0.8923  +       Copia   LTR     4471309

	($lLTR_length, $rLTR_length) = ($1, $2) if $nextline =~ /lLTR: ([0-9a-z]+)\s+rLTR: ([0-9a-z]+)/i;
	next unless $lLTR_length ne "NA" and $rLTR_length ne "NA" and $lLTR_length > 0 and $rLTR_length > 0;
#        Adjust: NO      lLTR: NA        rLTR: NA
#        Adjust: 3' lLTR lLTR: 810       rLTR: 804
#        Adjust: NO      lLTR: 0 rLTR: 0

	$anno{$loc} = $supfam;
	}
close Anno;


# read target seq and print formatted names
$/ = "\n>";
while (<Seq>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id =~ s/\s+.*//;
	$id =~ s/#.*//;
	$seq =~ s/\s+//g;
	my ($chr, $lLTR_start, $rLTR_end);
	($chr, $lLTR_start, $rLTR_end) = ($1, $2, $3) if $id =~ /^(.*):([0-9]+)..([0-9]+)$/;
	($lLTR_start, $rLTR_end) = ($rLTR_end, $lLTR_start) if $rLTR_end < $lLTR_start;
	my $anno = "unknown";
	next unless defined $anno{$id}; #skip the seqs that are not in %anno (duplicated, heavily nested, or lack of LTR structure)
	$anno = $anno{$id};
	print ">$id#LTR/$anno\n$seq\n";
	}
$/ = "\n";
close Seq;

