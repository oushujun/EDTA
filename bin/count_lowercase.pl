#!/usr/bin/perl
# facilitated by ChatGPT
# Shujun Ou (shujun.ou.1@gmail.com)
# 03/15/2024

use strict;
use warnings;

# Check if a file path is provided
my $filename = $ARGV[0];
if (not defined $filename) {
    die "Usage: $0 <filename>\n";
}

# Open the FASTA file
open(my $fh, '<', $filename) or die "Could not open file '$filename' $!";

my $count = 0; # Initialize counter for lowercase base pairs

while (my $line = <$fh>) {
    chomp $line; # Remove newline characters
    # Skip sequence identifier lines
    next if $line =~ /^>/;
    # Count lowercase characters
    $count += () = $line =~ /[a-z]/g;
}

close($fh);

print "Number of lowercase base pairs: $count\n";

