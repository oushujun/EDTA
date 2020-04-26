#!/usr/bin/env perl
use strict;
use warnings;

# Developed by Yi Liao (yiliao1022@gmail.com, 05/14/2019, UCI)
# Filter out blast self-alignment and convert blast results to the RECON MSP format.


my $usage = "blastn -outfmt 6 ... | $0 - > blast_msp_recon.txt\n-OR-\n".
    "$0 blasttable.bln > blast_msp_recon.txt\n";
die $usage if !@ARGV;

open IN, "$ARGV[0]" or die "$!";

while (<IN>) {
    chomp;
    my @f = split;
    die "This blast report is not formatted correctly. Exiting.\n"
        unless @f == 12;
    # comment the next line out if compairing one sequence to itself (e.g., chrI -> chrI)
    next if (($f[0] eq $f[1]) and (($f[6] == $f[8]) | ($f[7] ==$f[9])));
    if ($f[2] == 100.00) {
    printf("%06d %03d %08d %08d %s %08d %08d %s \n",
           $f[11], $f[2], $f[6], $f[7], $f[0], $f[8], $f[9], $f[1]);
    } else {
    printf("%06d %.1f %08d %08d %s %08d %08d %s \n", 
           $f[11], $f[2], $f[6], $f[7], $f[0], $f[8], $f[9], $f[1]);
    }
}
