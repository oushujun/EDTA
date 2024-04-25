#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# Add ID= field from the gff3 9th column to intact fasta files by matching coordinates
# Shujun Ou (shujun.ou.1@gmail.com) 03/26/2024 Facilitated by ChatGPT

# Variables for file paths
my ($fasta_file, $gff_file);

# Parse command line arguments
GetOptions(
    'fa=s' => \$fasta_file,
    'gff=s' => \$gff_file
) or die "Usage: $0 -fa <fasta_file> -gff <gff_file>\n";

# Verify that both files are provided
die "Usage: $0 -fa <fasta_file> -gff <gff_file>\n" unless $fasta_file && $gff_file;

# Hash to store GFF3 coordinates and IDs
my %gff_data;

# Process the GFF3 file
open(my $gff_fh, '<', $gff_file) or die "Could not open GFF3 file $gff_file: $!";
while (my $line = <$gff_fh>) {
    chomp $line;
    next if $line =~ /^\s*$/ || $line =~ /^#/; # Skip empty lines and comments
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $line;
    
    # Ensure start is always less than end
    ($start, $end) = ($start < $end) ? ($start, $end) : ($end, $start);

    # Extract ID from attributes
    if ($attributes =~ /ID=([^;]+)/) {
        my $id = $1;
        # Store in hash with coordinates as key and ID as value
        $gff_data{"$chr:$start..$end"} = $id;
    }
}
close $gff_fh;

# Process the FASTA file
open(my $fasta_fh, '<', $fasta_file) or die "Could not open FASTA file $fasta_file: $!";
while (my $line = <$fasta_fh>) {
    chomp $line;
    if ($line =~ /^>/) { # Sequence ID line
        my $original_line = $line; # Keep the original line
        # Extract coordinates
        if ($line =~ /(.*)\|(\w+):(\d+)\.\.(\d+)/) {
            my ($seq_id, $chr, $start, $end) = ($1, $2, $3, $4);
            
            # Normalize start and end
            ($start, $end) = ($start < $end) ? ($start, $end) : ($end, $start);
            
            # Check if normalized coordinates are in GFF3 data
            if (exists $gff_data{"$chr:$start..$end"}) {
                # Append GFF3 ID to sequence ID
                $line = "$original_line ID=$gff_data{\"$chr:$start..$end\"}";
            }
        }
    }
    print "$line\n"; # Print modified or original line
}
close $fasta_fh;

