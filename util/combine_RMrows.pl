#!/usr/bin/env perl
# facilitated by ChatGPT
# Shujun Ou (shujun.ou.1@gmail.com)
# 09/20/2023

use strict;
use warnings;

my $usage = "
Combine RepeatMasker entries that appear to be from the same repeat based on both genome and repeat coordinates.
	perl combine_RMrows.pl -rmout genome.fa.out [options] > genome.fa.cmb.out
		-rmout [file]	Repeatmasker .out file for merging. Required.
		-maxgap [int]	Maximum physical distance between the two adjacent annotation
				on the genome and on the TE. Default: 35 (bp)
		-maxdiv [float]	Maximum divergence between the two adjacent annotation. 
				Default: 3.5 (%)
		-iter [int]	Number of iterations to combine rows. Default: automatic
\n";

# parameter setting
my $max_gap = 35;
my $max_div = 3.5;
my $iter = 1; 
my $user_iter = 0;
my $rmout = ''; # Repeatmasker .out file for merging

my $k=0;
foreach (@ARGV){
        $max_gap=$ARGV[$k+1] if /^-maxgap$/i and $ARGV[$k+1] !~ /^-/;
        $max_div=$ARGV[$k+1] if /^-maxdiv$/i and $ARGV[$k+1] !~ /^-/;
        $rmout=$ARGV[$k+1] if /^-rmout$/i and $ARGV[$k+1] !~ /^-/;
	$user_iter=$ARGV[$k+1] if /^-iter$/i and $ARGV[$k+1] !~ /^-/;
        $k++;
        }

# checks
die "\nERROR: Input sequence file not exist!\n$usage" unless -s $rmout;
die "\nERROR: The -iter parameter receives non-integer input!\n$usage" unless $user_iter =~ /^[0-9]+$/;

# output combine logs
open LOG, ">$rmout.log" or die $usage;

# itreatively combine rows appear to derive from the same repeat
my $num_log = 0; # count log lines at the end of each iteration
my $next = 0;
$iter = $user_iter if $user_iter != 0;
`cp $rmout $rmout.iter0`;
for (my $i=0; $i<$iter; $i++){
	my $date=`date`;
	chomp ($date);
	print "$date\tCombine fragmented repeats. Working on iteration $i\n";

	# write temp results to file
	open RMout, "<$rmout.iter$i" or die $usage;
	$next = $i + 1;
	open Out, ">$rmout.iter$next" or die $!;

# print header
print Out "SW_score\tperc_div.\tperc_del.\tperc_ins.\tquery_sequence\tquery_begin\tquery_end\tquery_remain\tstrand\tmatching_repeat\trepeat_class/family\trepeat_begin\trepeat_end\trepeat_remain\tID\n";

my %prev_row;
while (<RMout>) {
	chomp;
	s/^\s+//;
	s/\s+/\t/g;
	next if /^$/;
	next unless /^[0-9]+/;
	s/[\(\)]+//g;

	my ($SW_score, $div, $del, $ins, $chr, $start, $end, $chr_remain, $strand, $element, $TE_class, $element_start, $element_end, $element_remain, $rm_ID) = split;
	if (%prev_row && $prev_row{'chr'} eq $chr && $prev_row{'strand'} eq $strand && $prev_row{'element'} eq $element && $prev_row{'TE_class'} eq $TE_class
		&& abs($start - $prev_row{'end'}) <= $max_gap && abs($prev_row{'div'} - $div) <= $max_div
		&& (($prev_row{'strand'} eq '+' && ($prev_row{'element_end'} - $element_start) <= $max_gap) 
			or ($prev_row{'strand'} eq 'C' && ($element_end - $prev_row{'element_remain'}) <= $max_gap))) {

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
        my $combined_rm_ID = "$prev_row{'rm_ID'}_$rm_ID";
	my $combined_element_start;
	my $combined_element_end;
	my $combined_element_remain;
	if ($strand eq '+') {
		$combined_element_start = $prev_row{'element_start'};
        	$combined_element_end = $element_end;
	        $combined_element_remain = $element_remain;
		} else {
		$combined_element_start = $prev_row{'element_start'};
		$combined_element_end = $prev_row{'element_end'};
		$combined_element_remain = $element_remain;
		}
        
	if ($prev_row{'strand'} ne 'C'){
		print Out "$combined_SW_score\t$combined_div\t$combined_del\t$combined_ins\t$prev_row{'chr'}\t$combined_start\t$combined_end\t$combined_chr_remain\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t$combined_element_start\t$combined_element_end\t($combined_element_remain)\t$combined_rm_ID\n";
		} else {
		print Out "$combined_SW_score\t$combined_div\t$combined_del\t$combined_ins\t$prev_row{'chr'}\t$combined_start\t$combined_end\t$combined_chr_remain\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t($combined_element_start)\t$combined_element_end\t$combined_element_remain\t$combined_rm_ID\n";
		}

	# print out merged lines
	print LOG "Row1:\t".join("\t", @prev_row{qw(SW_score div del ins chr start end chr_remain strand element TE_class element_start element_end element_remain rm_ID)}), "\n";
	print LOG "Row2:\t$_\n";
	print LOG "Merged:\t$combined_SW_score\t$combined_div\t$combined_del\t$combined_ins\t$prev_row{'chr'}\t$combined_start\t$combined_end\t$combined_chr_remain\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t$combined_element_start\t$combined_element_end\t$combined_element_remain\t$combined_rm_ID\n\n";
	
        %prev_row = ();
    }
    else {

	# this row cannot combine with the previous row, print out the line
	if (%prev_row and defined $prev_row{'rm_ID'}){
		if ($prev_row{'strand'} ne 'C'){
			print Out "$prev_row{'SW_score'}\t$prev_row{'div'}\t$prev_row{'del'}\t$prev_row{'ins'}\t$prev_row{'chr'}\t$prev_row{'start'}\t$prev_row{'end'}\t$prev_row{'chr_remain'}\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t$prev_row{'element_start'}\t$prev_row{'element_end'}\t($prev_row{'element_remain'})\t$prev_row{'rm_ID'}\n";
			} else {
			print Out "$prev_row{'SW_score'}\t$prev_row{'div'}\t$prev_row{'del'}\t$prev_row{'ins'}\t$prev_row{'chr'}\t$prev_row{'start'}\t$prev_row{'end'}\t$prev_row{'chr_remain'}\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t($prev_row{'element_start'})\t$prev_row{'element_end'}\t$prev_row{'element_remain'}\t$prev_row{'rm_ID'}\n";
			}
		}

# store current row info to %prev_row before moving to the next row
        %prev_row = ('SW_score' => $SW_score, 'div' => $div, 'del' => $del, 'ins' => $ins, 'chr' => $chr, 'start' => $start, 'end' => $end, 'chr_remain' => $chr_remain, 'strand' => $strand, 'element' => $element, 'TE_class' => $TE_class, 'element_start' => $element_start, 'element_end' => $element_end, 'element_remain' => $element_remain, 'rm_ID' => $rm_ID);
    }
}

# Print the last row if it's not combined
if (%prev_row and $prev_row{'strand'} ne 'C'){
	print Out "$prev_row{'SW_score'}\t$prev_row{'div'}\t$prev_row{'del'}\t$prev_row{'ins'}\t$prev_row{'chr'}\t$prev_row{'start'}\t$prev_row{'end'}\t$prev_row{'chr_remain'}\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t$prev_row{'element_start'}\t$prev_row{'element_end'}\t($prev_row{'element_remain'})\t$prev_row{'rm_ID'}\n";
	} else {
	print Out "$prev_row{'SW_score'}\t$prev_row{'div'}\t$prev_row{'del'}\t$prev_row{'ins'}\t$prev_row{'chr'}\t$prev_row{'start'}\t$prev_row{'end'}\t$prev_row{'chr_remain'}\t$prev_row{'strand'}\t$prev_row{'element'}\t$prev_row{'TE_class'}\t($prev_row{'element_start'})\t$prev_row{'element_end'}\t$prev_row{'element_remain'}\t$prev_row{'rm_ID'}\n";
	}

	# end of the iteration
	close RMout;
	close Out;

	# automatically increase iteration based on the log result
	my $curr_log = `wc -l "$rmout.log"`;
	$curr_log = (split /\s+/, $curr_log)[0];
	if ($num_log == $curr_log){
		print "Saturated at iter$i, automatically stop.\n\n";
		last;
	} else {
		$num_log = $curr_log;
		$iter++ if $user_iter == 0;
	}
}

# copy the last iteration as the final file
`cp $rmout.iter$next $rmout.cmb`;
close LOG;
