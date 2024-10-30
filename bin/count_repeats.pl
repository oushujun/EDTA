#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "";

my $RM_out = $ARGV[0];
my $mincov = 0.8;

open RM, "<$RM_out" or die $usage;
open Complete, ">$RM_out.complete.list";
open Fragment, ">$RM_out.fragment.list";

while (<RM>){
	s/^\s+//;
	chomp;
	s/[\(\)]+//g;
	next unless /^[0-9]+/;
	my ($id, $from, $to, $te_id, $te_class, $te_from, $te_to, $te_left) = (split)[4,5,6,9,10,11,12,13];
	my $target_len = $te_to - $te_from + 1;
	my $uncovered = $te_from + $te_left;
	my $cov = $target_len / ($target_len + $uncovered);
	if ($cov >= $mincov){
		print Complete "$id\t$from\t$to\t$te_id\t$te_class\n";
		} else {
		print Fragment "$id\t$from\t$to\t$te_id\t$te_class\n";
		}
	}

# 2352    4.5  1.0  1.9  NIP_Chr10_10008421_10018854   7444  7757  (2677) + Os0072                    MITE/Stow             3    313     (0)     1  
#  1239   17.1  0.5  0.0  NIP_Chr10_10008421_10018854   8622  8831  (1603) + Os0205                    DNAnona/MULEtir       1    211     (1)     2  
#   1019   18.8  1.9  0.5  NIP_Chr10_10008421_10018854   9143  9351  (1083) C Os0205                    DNAnona/MULEtir     (0)    212       1     3  
#   19418    1.1  0.9  0.1  NIP_Chr10_10189714_10203541      1  2194 (11634) C Os1279                    DNAauto/CACTA     (404)  13198   10989     4  
#   27716    0.9  0.1  1.2  NIP_Chr10_10189714_10203541   2240  5342  (8486) C Os1279                    DNAauto/CACTA    (2614)  10988    7920     4  
#   6
