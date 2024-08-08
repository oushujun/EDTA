#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

# Add "###" labels before and after intact LTR elements
# Usage: perl format_intact_gff3.pl -f EDTA.intact.gff3 > EDTA.intact.formatted.gff3
# Shujun Ou (shujun.ou.1@gmail.com) 03/26/2024 Facilitated by ChatGPT


my %opts;
getopts('f:', \%opts);

unless ($opts{f}) {
    die "Usage: $0 -f <filename>\n";
}

my $filename = $opts{f};
open(my $fh, '<', $filename) or die "Cannot open file $filename: $!";

my $flag = 0;
my $count = 0;
my $last_was_hash = 0;  # Variable to track if the last printed line was "###"

while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /ID=repeat_region_\d+/) {
        unless ($last_was_hash) {  # Check if the last line printed was "###"
            print "###\n";  # Print ### before the matched line
            $last_was_hash = 1;  # Set this to indicate that "###" was just printed
        }
        print "$line\n";
        $flag = 1;  # Set flag to start counting the next six lines
        $count = 0;  # Reset counter
    } elsif ($flag) {
        print "$line\n";
        $last_was_hash = 0;  # Reset this since we're printing normal lines now
        $count++;
        if ($count == 5) {  # After printing six lines, print ### and reset the flag
            print "###\n";
            $flag = 0;
            $last_was_hash = 1;  # Indicate that "###" was just printed
        }
    } else {
        if ($line ne "###") {  # Only print the line if it's not "###" to avoid consecutive hashes
            print "$line\n";
            $last_was_hash = 0;  # Reset this since we're printing normal lines
        }
    }
}

close $fh;

