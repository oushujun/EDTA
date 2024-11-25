#!/usr/bin/env perl

# Script to process GFF3 files and add "###" separators for each repeat ID, keeping the header
# Usage: perl format_intact_gff3.pl -f EDTA.intact.gff3 > EDTA.intact.formatted.gff3
# Shujun Ou (shujun.ou.1@gmail.com) 07/07/2024 Facilitated by ChatGPT

use strict;
use warnings;

# Check if filename is provided
if (@ARGV != 1) {
    die "Usage: $0 <input_gff3_file>\n";
}

my $input_file = $ARGV[0];

# Open the input file or read from STDIN
my $fh;
if ($input_file eq '-') {
    open($fh, '<&', 'STDIN') or die "Could not read from STDIN: $!";
} else {
    open($fh, '<', $input_file) or die "Could not open file '$input_file' $!";
}

my %repeats;
my @order;
my @header;
my $has_printed = 0;

# Read the input file line by line
while (my $line = <$fh>) {
    chomp $line;

    if ($line =~ /^\s*$/) {
        next;  # Skip empty lines
    } elsif ($line =~ /^#/) {
        next if $line eq '###';  # Skip "###" lines in the header
        push @header, $line;  # Store header lines
        next;
    }

    # Split the line into columns
    my @fields = split "\t", $line;
    my $attributes = $fields[8];
    my %attr_hash = map { split /=/, $_, 2 } grep { /=/ } split /;/, $attributes;  # Ensure valid key-value pairs

    if (exists $attr_hash{'Parent'}) {
        # It's a child entry
        push @{$repeats{$attr_hash{'Parent'}}{'children'}}, \@fields;
    } elsif (exists $attr_hash{'ID'}) {
        # It's a repeat entry
        push @order, $attr_hash{'ID'} unless exists $repeats{$attr_hash{'ID'}};
        $repeats{$attr_hash{'ID'}}{'parent'} = \@fields;
    }
}

close($fh);

# Sort the entries based on sequence ID and coordinates
@order = sort {
    $repeats{$a}{'parent'}[0] cmp $repeats{$b}{'parent'}[0] ||
    $repeats{$a}{'parent'}[3] <=> $repeats{$b}{'parent'}[3]
} @order;

# Print the header
foreach my $header_line (@header) {
    print $header_line, "\n";
}

# Print the grouped entries with separators
foreach my $id (@order) {
    if (exists $repeats{$id}{'parent'}) {
        print "###\n" if $has_printed;
        $has_printed = 1;
        print join("\t", @{$repeats{$id}{'parent'}}), "\n";
        if (exists $repeats{$id}{'children'}) {
            foreach my $child (sort {
                $a->[0] cmp $b->[0] ||
                $a->[3] <=> $b->[3]
            } @{$repeats{$id}{'children'}}) {
                print join("\t", @$child), "\n";
            }
        }
    }
}

print "###\n" if $has_printed;  # Print the final separator only if entries were printed

__END__

