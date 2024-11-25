#!/usr/bin/env perl
use strict;
use warnings;

# Calculate divergence table based on superfamilies and EDTA.TEanno.bed file
# Usage: perl div_table.pl $genome.mod.EDTA.anno/$genome.mod.EDTA.TEanno.bed $genome.fa $outfile_prefix
# v0.3 11/16/2023 Shujun Ou (shujun.ou.1@gmail.com), facilitated by ChatGPT 3.5

my $usage = "
Calculate divergence table based on superfamilies and EDTA.TEanno.bed file
Usage: 
  EDTA input:
    perl div_table.pl \$genome.mod.EDTA.anno/\$genome.mod.EDTA.TEanno.bed \$genome.fa \$outfile_prefix
  RepeatMasker input:
    perl div_table.pl \$genome.out \$genome.fa RMout
    * RMout is a fixed string that indicates the input is from RepeatMasker .out file.
\n";

# Define the input file names
my $bed_file = $ARGV[0];
my $genome_file = $ARGV[1];
my $outfile_prefix = $ARGV[2] // die "No output file prefix provided\n$usage";

# Define the input file types
my $EDTAbed = 1; # The default input is EDTA.TEanno.bed. If set to 0, the input is expected RepeatMasker .out
$EDTAbed = 0 if defined $ARGV[2] and $ARGV[2] eq "RMout";

# Define the div ranges
my $max_div = 50;

# Hash to store the total lengths for each class within each div range
my %summary_table;
my %class;

# Read the genome file to calculate genome size
open my $genome_fh, '<', $genome_file or die "Cannot open genome file $genome_file\n$usage";
my $genome_size = 0;
while (my $genome_line = <$genome_fh>) {
    chomp $genome_line;

    # Skip lines starting with ">"
    next if $genome_line =~ /^>/;

    # Remove whitespace
    $genome_line =~ s/\s+//g;

    # Add the length of the nucleotide sequence to the genome size
    $genome_size += length($genome_line);
}
close $genome_fh;

# Open the input file
open my $fh, '<', $bed_file or die "Cannot open input $bed_file\n$usage";

# Read the file line by line
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/^\s+//;
    next if $line =~ /^$/;

    # Extract the fields from the line
    my ($chr, $element_start, $element_end, $TE_ID, $TE_class, $method, $iden, $score, $strand, $phase, $extra, $type, $class, $div);
    if ($EDTAbed == 1) {
        ($chr, $element_start, $element_end, $TE_ID, $TE_class, $method, $iden, $score, $strand, $phase, $extra, $type, $class) = split /\s+/, $line;
        next unless $iden =~ /^[0-9.]+$/;
        $div = 100 * (1 - $iden);
    } elsif ($EDTAbed == 0) {
        $line =~ s/^\s+//;
        next unless $line =~ /^[0-9]+/;
        ($score, $div, $chr, $element_start, $element_end, $strand, $TE_ID, $class) = (split /\s+/, $line)[0,1,4,5,6,8,9,10];
    }

    # Calculate the length of the sequence
    my $length = $element_end - $element_start + 1;

    # Define divergence range
    my $div_range = int($div);

    # Update the summary table
    $summary_table{$div_range}{$class} += $length;
    $class{$class} = $class;
}

# Open the output file for writing the summary table
open my $summary_fh, '>', "$outfile_prefix.div" or die "Cannot open output file for summary format\n";

# Print the summary table header
print $summary_fh "Div\t" . join("\t", map { "$_" } sort keys %class) . "\tgenome_size\n";

# Print the summary table rows
foreach my $div (0..$max_div) {
    print $summary_fh "$div";
    foreach my $class (sort keys %class) {
        my $length = $summary_table{$div}{$class} || 0;
        print $summary_fh "\t$length";
    }
    print $summary_fh "\t$genome_size\n";
}

# Close the input and summary output files
close $fh;
close $summary_fh;

# Process the summary file to generate the div_long output
open my $div_fh, '<', "$outfile_prefix.div" or die "Cannot open div file for reading\n";

# Open the output file for long format
open my $long_fh, '>', "$outfile_prefix.div_long" or die "Cannot open output file for long format\n";

# Print column headers
print $long_fh "div\tsupfam\tbp\tgenome_size\tpcnt\n";

# Read the header line in the div file and extract supfam names
my $header_line = <$div_fh>;
chomp $header_line;
my @supfams = split /\t/, $header_line;
shift @supfams; # remove 'Div'
pop @supfams; # remove 'genome_size'

# Read the div file and process it
while (my $line = <$div_fh>) {
    chomp $line;
    my @fields = split /\t/, $line;
    my $div = shift @fields;
    my $genome_size = pop @fields;
    for (my $i = 0; $i < scalar(@fields); $i++) {
        my $bp = $fields[$i];
        my $supfam = $supfams[$i];
        my $pcnt = sprintf("%.4f", ($bp / $genome_size) * 100);
        print $long_fh "$div\t$supfam\t$bp\t$genome_size\t$pcnt\n" if $bp > 0;
    }
}

# Close the output files
close $div_fh;
close $long_fh;
