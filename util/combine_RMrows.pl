#!/usr/bin/env perl
# facilitated by ChatGPT
# Shujun Ou

use strict;
use warnings;

my $max_gap = 50; # Set this according to your requirements
my $max_div = 5; # Set this according to your requirements

# print header
print "SW_score\tdiv\tdel\tins\tchr\tstart\tend\tchr_remain\tstrand\telement\tTE_class\telement_start\telement_end\telement_remain\trm_ID\n";

my %prev_row;
while (<>) {
    chomp;
    #    print $_ and next if /^$/;
    my ($SW_score, $div, $del, $ins, $chr, $start, $end, $chr_remain, $strand, $element, $TE_class, $element_start, $element_end, $element_remain, $rm_ID) = split;

    if (%prev_row && $prev_row{'chr'} eq $prev_row{'strand'} eq $strand && $prev_row{'element'} eq $element && $prev_row{'TE_class'} eq $TE_class
        && ($start - $prev_row{'end'}) <= $max_gap 
	&& abs($prev_row{'div'} - $div) <= $max_div) {

        # Calculate weights based on the length
        my $prev_length = $prev_row{'end'} - $prev_row{'start'};
        my $current_length = $end - $start;
        my $total_length = $prev_length + $current_length;

        # Weighted averaging
        my $combined_SW_score = sprintf("%.0f", ($prev_row{'SW_score'} * $prev_length + $SW_score * $current_length) / $total_length);
        my $combined_div = sprintf("%.1f", ($prev_row{'div'} * $prev_length + $div * $current_length) / $total_length);
        my $combined_del = sprintf("%.1f", ($prev_row{'del'} * $prev_length + $del * $current_length) / $total_length);
        my $combined_ins = sprintf("%.1f", ($prev_row{'ins'} * $prev_length + $ins * $current_length) / $total_length);

        # Combine coordinates
        my $combined_start = $prev_row{'start'};
        my $combined_end = $end;

        # Other fields
        my $combined_chr_remain = $chr_remain;
        my $combined_rm_ID = $rm_ID;
	my $combined_element_start;
	my $combined_element_end;
	my $combined_element_remain;
	if (1) {
		#if ($strand eq '+') {
		$combined_element_start = $prev_row{'element_start'};
        	$combined_element_end = $element_end;
	        $combined_element_remain = $element_remain;
	}
	#} else {
	#}
        
        print join("\t", $combined_SW_score, $combined_div, $combined_del, $combined_ins, $prev_row{'chr'}, $combined_start, $combined_end, $combined_chr_remain, $prev_row{'strand'}, $prev_row{'element'}, $prev_row{'TE_class'}, $combined_element_start, $combined_element_end, $combined_element_remain, $combined_rm_ID), "\n";
        
        %prev_row = ();
    }
    else {
	 print join("\t", @prev_row{qw(SW_score div del ins chr start end chr_remain strand element TE_class element_start element_end element_remain rm_ID)}), "\n" if %prev_row and defined $prev_row{'rm_ID'};
        %prev_row = ('SW_score' => $SW_score, 'div' => $div, 'del' => $del, 'ins' => $ins, 'chr' => $chr, 'start' => $start, 'end' => $end, 'chr_remain' => $chr_remain, 'strand' => $strand, 'element' => $element, 'TE_class' => $TE_class, 'element_start' => $element_start, 'element_end' => $element_end, 'element_remain' => $element_remain, 'rm_ID' => $rm_ID);
    }
}
print join("\t", @prev_row{qw(SW_score div del ins chr start end chr_remain strand element TE_class element_start element_end element_remain rm_ID)}), "\n" if %prev_row; # Print the last row if it's not combined

