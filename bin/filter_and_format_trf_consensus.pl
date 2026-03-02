#!/usr/bin/env perl
use strict;
use warnings;

# keep consensus trf with more than 50 copies; dupliate the seed sequence until >= 20bp; output unique sequences
# Shujun Ou (shujun.ou.1@gmail.com) facilitated by ChatGPT
# 10/19/2025

my %seen;  # store unique consensus sequences

while (<>) {
    next if /^#/;  # skip header
    chomp;
    my @f = split /\t/;
    my ($id, $n_copy, $cons) = ($f[1], $f[9], $f[13]);
    next unless $n_copy >= 50;
    next unless defined $cons and $cons ne '';

    # ensure consensus at least 20 bp long
    my $format_cons = $cons;
    while (length($format_cons) < 20) {
        $format_cons .= $cons;
    }

    # skip duplicates
    next if $seen{$format_cons}++;
    
    print ">$id\n$format_cons\n";
}

