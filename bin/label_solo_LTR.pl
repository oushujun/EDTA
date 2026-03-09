#!/usr/bin/env perl
use warnings;
use strict;

## Identify solo LTR candidates from RepeatMasker output using LTR boundary info
## Solo LTRs are genomic hits that align to the terminal LTR portion of a whole-element library entry
## but do NOT extend significantly into the internal region.
##
## Usage: perl label_solo_LTR.pl -bound <LTRbound> -rmout <RM.out> [-maxint 100] > solo_LTR.list
## Output: tab-delimited list of solo LTR candidates (genomic location, library element, LTR side)
##
## Shujun Ou (03-07-2026)

my $bound_file = '';
my $rmout_file = '';
my $maxint = 100; # max bp of internal region allowed in the match to still call it solo LTR

my $k = 0;
foreach (@ARGV){
	$bound_file = $ARGV[$k+1] if /^-bound$/i;
	$rmout_file = $ARGV[$k+1] if /^-rmout$/i;
	$maxint = $ARGV[$k+1] if /^-maxint$/i;
	$k++;
}

die "Usage: perl label_solo_LTR.pl -bound <LTRbound> -rmout <RM.out> [-maxint $maxint]\n"
	unless -s $bound_file and -s $rmout_file;

# Read boundary file: TE_name => [total_len, lLTR_len, rLTR_len]
my %bound;
open BOUND, "<$bound_file" or die "Cannot open $bound_file: $!\n";
while (<BOUND>){
	chomp;
	my ($name, $total_len, $lLTR_len, $rLTR_len) = split /\t/;
	next unless defined $rLTR_len;
	# Strip #class for matching RM output (RM uses name#class format)
	$bound{$name} = [$total_len, $lLTR_len, $rLTR_len];
}
close BOUND;

# Parse RepeatMasker .out file and identify solo LTR candidates
open RM, "<$rmout_file" or die "Cannot open $rmout_file: $!\n";

# Print header
print "#solo_LTR_chr\tsolo_LTR_start\tsolo_LTR_end\tlib_element\tLTR_side\tlib_start\tlib_end\tlib_total_len\n";

while (<RM>){
	s/^\s+//;
	next if /^$/ or /^\s/ or /^SW|^score/;
	my @f = split /\s+/;
	next unless @f >= 15;

	my ($score, $div, $del, $ins, $chr, $qstart, $qend, $strand, $repeat, $class, $rstart, $rend, $rleft) =
		@f[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13];

	my $full_name = "$repeat#$class";

	# Check if this repeat is in the boundary file
	next unless defined $bound{$full_name};
	my ($total_len, $lLTR_len, $rLTR_len) = @{$bound{$full_name}};

	# Determine the match range in library coordinates
	# RM .out format: for + strand: rstart=begin, rend=end, rleft=(left)
	# For C strand: rstart=(left), rend=end, rleft=begin
	my ($lib_start, $lib_end);
	if ($strand eq '+'){
		$lib_start = $rstart;
		$lib_end = $rend;
	} else {
		# complement strand: rleft is in parentheses format
		$rleft =~ s/[()]//g;
		$rstart =~ s/[()]//g;
		$lib_start = $rleft;
		$lib_end = $rend;
	}

	# Check if the match is confined to the left LTR region
	# Left LTR: positions 1 to lLTR_len
	# Internal: positions lLTR_len+1 to total_len-rLTR_len
	# Right LTR: positions total_len-rLTR_len+1 to total_len
	my $int_start = $lLTR_len + 1;
	my $int_end = $total_len - $rLTR_len;

	# Calculate how far into internal region the match extends
	my $int_overlap = 0;
	if ($lib_end > $int_start && $lib_start < $int_end){
		# match overlaps internal region
		my $ov_start = $lib_start > $int_start ? $lib_start : $int_start;
		my $ov_end = $lib_end < $int_end ? $lib_end : $int_end;
		$int_overlap = $ov_end - $ov_start + 1;
	}

	# Solo LTR: match is mostly in LTR region with minimal internal overlap
	next if $int_overlap > $maxint;

	# Determine which LTR side
	my $side = "unknown";
	if ($lib_end <= $lLTR_len + $maxint){
		$side = "left_LTR";
	} elsif ($lib_start >= $int_end + 1 - $maxint){
		$side = "right_LTR";
	} elsif ($lib_start <= $lLTR_len && $lib_end >= $int_end + 1){
		# spans both LTRs but no internal — unlikely, skip
		next;
	}

	print "$chr\t$qstart\t$qend\t$full_name\t$side\t$lib_start\t$lib_end\t$total_len\n";
}
close RM;