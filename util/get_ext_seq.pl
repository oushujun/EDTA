#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;

my $genome = $ARGV[0];
my $seq = $ARGV[1];
my $extlen = 30;
my $ori_name = 1; #1 will use the original coordinate as seq name.
#my $call_seq = "~/las/git_bin/EDTA/util/call_seq_by_list.pl ";
my $call_seq = "$FindBin::Bin/call_seq_by_list.pl";

open Seq, "<$seq" or die $!;
open Out, ">$seq.ext$extlen.list" or die $!;
while (<Seq>){
	next unless /^>/;
	s/>//g;
	my ($chr, $from, $to);
	($chr, $from, $to) = ($1, $2, $3) if /^(\S+)[:_\-]+([0-9]+)[_\.\:\-]+([0-9]+)/; #eg.: Chr10:3119406..3119688; Chr10_3119406_3119688; Chr10:3119406:3119688
	my ($from_ext, $to_ext) = ($from - $extlen, $to + $extlen);
	$from_ext = 1 if $from_ext < 1;
	print Out "$chr:$from..$to\t$chr:$from_ext..$to_ext\n";
	}
close Seq;
close Out;

## Get extended fasta seq
`perl $call_seq $seq.ext$extlen.list -C $genome > $seq.ext$extlen.fa`;
`perl -i -nle 's/>.*\\|/>/; print \$_' $seq.ext$extlen.fa` if $ori_name == 1;
