#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;

my $usage = "\nSummarize RepeatMasker .out file results based on the genome name\n
	perl RM_summary.pl genome.fa\n";

my $genome = $ARGV[0];
my $script_path = $FindBin::Bin;

die "No fasta specified!\n$usage" unless -s $genome;
die "No RM files associated with the $genome found!\n" unless -s "$genome.out";

my $size_info = `perl $script_path/count_base.pl $genome`;
my ($total, $missing) = (split /\s+/, $size_info)[1,2];
my $actual = $total - $missing;
my ($helitron, $tir, $nonltr, $ltr) = (0, 0, 0, 0);


#count helitron
$helitron = `awk '{if (\$11~/Helitron/) print \$5"\\t"\$6"\\t"\$7}' $genome.out | sort -suV | perl $script_path/count_mask.pl`;

#count non-heli TIRs
$tir = `awk '{if (\$11~/DNA|MITE/ && \$11!~/Helitron/) print \$5"\\t"\$6"\\t"\$7}' $genome.out | sort -suV | perl $script_path/count_mask.pl`;

#count SINE and LINE
$nonltr = `awk '{if (\$11~/SINE|LINE/) print \$5"\\t"\$6"\\t"\$7}' $genome.out | sort -suV | perl $script_path/count_mask.pl`;

#count LTR
$ltr = `awk '{if (\$11~/LTR/) print \$5"\\t"\$6"\\t"\$7}' $genome.out | sort -suV | perl $script_path/count_mask.pl`;

my $hel_frac = sprintf("%.2f%", $helitron/$actual*100);
my $tir_frac = sprintf("%.2f%", $tir/$actual*100);
my $ltr_frac = sprintf("%.2f%", $ltr/$actual*100);
my $nonltr_frac = sprintf("%.2f%", $nonltr/$actual*100);


print "Total(bp): $total\n\nNonLTR(bp): $nonltr\nLTR(bp): $ltr\nTIR(bp): $tir\nHelitron(bp): $helitron\n\n";
print "nonLTR: $nonltr_frac\nLTR: $ltr_frac\nTIR: $tir_frac\nHelitron: $hel_frac\n\n";
print "Category: nonLTR LTR TIR Helitron\nPercent: $nonltr_frac $ltr_frac $tir_frac $hel_frac\n";

