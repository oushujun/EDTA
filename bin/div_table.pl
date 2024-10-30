#!/usr/bin/env perl
use strict;
use warnings;

# Calculate divergence table based on superfamilies and EDTA.TEanno.bed file
# Usage: perl div_table.pl $genome.mod.EDTA.anno/$genome.mod.EDTA.TEanno.bed > $genome.mod.EDTA.TEanno.div
# v0.1 07/13/2023 Shujun Ou (shujun.ou.1@gmail.com), facilitated by ChatGPT 3.5

my $usage = "
Calculate divergence table based on superfamilies and EDTA.TEanno.bed file
Usage: 
  EDTA input:
	perl div_table.pl \$genome.mod.EDTA.anno/\$genome.mod.EDTA.TEanno.bed > \$genome.mod.EDTA.TEanno.div
  RepeatMasker input:
  	perl div_table.pl \$genome.out RMout > \$genome.mod.EDTA.TEanno.div
		* RMout is a fixed string that indicates the input is from RepeatMasker .out file.
\n";

# Define the input file name
my $input_file = $ARGV[0];

# define the input file type
my $EDTAbed = 1; # The default input is EDTA.TEanno.bed. If set to 0, the input is expected RepeatMasker .out
$EDTAbed = 0 if defined $ARGV[1] and $ARGV[1] eq "RMout";

# Define the div ranges
my $max_div = 50;

# Hash to store the total lengths for each class within each div range
my %summary_table;
my %class;

# Open the input file
open my $fh, '<', $input_file or die "Cannot open input $input_file\n$usage";

# Read the file line by line
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/^\s+//;
    next if $line =~ /^$/;

    # Extract the fields from the line
    my ($chr, $element_start, $element_end, $TE_ID, $TE_class, $method, $iden, $score, $strand, $phase, $extra, $type, $class, $div);
    if ($EDTAbed == 1){
	    ($chr, $element_start, $element_end, $TE_ID, $TE_class, $method, $iden, $score, $strand, $phase, $extra, $type, $class) = split /\s+/, $line;
    	next unless $iden =~ /^[0-9.]+$/;
    	$div = 100 * (1 - $iden);
	} elsif ($EDTAbed == 0){
		$line =~ s/^\s+//;
		next unless $line =~ /^[0-9]+/;
		($score, $div, $chr, $element_start, $element_end, $strand, $TE_ID, $class) = (split /\s+/, $line)[0,1,4,5,6,8,9,10];
    	}

    # Calculate the length of the sequence
    my $length = $element_end - $element_start + 1;

    # define divergence range
    my $div_range = int($div);

    # Update the summary table
    $summary_table{$div_range}{$class} += $length;
    $class{$class} = $class;

}

# Print the summary table header
print "Div\t" . join("\t", map { "$_" } sort keys %class) . "\n";

# Print the summary table rows
foreach my $div (0..$max_div) {
    print "$div";
    foreach my $class (sort keys %class) {
        my $length = $summary_table{$div}{$class} || 0;
	#	$length = $summary_table{$div}{$class} if defined $summary_table{$div}{$class};
        print "\t$length";
    }
    print "\n";
}

# Close the input file
close $fh;

