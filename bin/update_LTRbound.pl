#!/usr/bin/env perl
use warnings;
use strict;

## Update the LTR boundary file with renamed TE IDs from rename_TE.pl, and reconcile
## the boundary lengths against the FINAL (post-filtering) library.
##
## Why reconciliation is needed:
##   LTR_retriever measures total_len/lLTR_len/rLTR_len on the PRE-filtering intact LTR
##   elements. EDTA advance filtering (TE_purifier / cleanup_tandem / cleanup_nested) can
##   trim those sequences afterwards, so the recorded coordinates may no longer fit the
##   sequence that actually ships. Left unreconciled, a large fraction of entries disagree
##   with the final library and many have lLTR_len + rLTR_len > total_len -- arithmetically
##   impossible, and it silently corrupts label_solo_LTR.pl (int_end = total_len - rLTR_len).
##
## Policy (conservative; correct-by-omission):
##   * A boundary row is emitted only when the shipped sequence was NOT trimmed
##     (final library length == recorded total_len) and the row is internally consistent
##     (lLTR_len + rLTR_len <= total_len). Such rows are emitted with the new TE ID.
##   * Every other row is DROPPED: its terminal-repeat coordinates no longer describe the
##     shipped sequence, and emitting a wrong boundary is worse than emitting none.
##   * The number of dropped rows is reported to STDERR (never silent).
##
## To also RECOVER boundaries for trimmed entries (higher coverage) instead of dropping
## them, re-measure lLTR/rLTR by self-aligning each final library sequence
## (see scripts/ltrbound_denovo.py in the oat workflow).
##
## Usage: perl update_LTRbound.pl <rename_map> <LTRbound> <final_lib.fa> > <updated_LTRbound>
##   rename_map format: TE_N#class \t original_name#class
##   LTRbound   format: original_name#class \t total_len \t lLTR_len \t rLTR_len

die "Usage: perl update_LTRbound.pl <rename_map> <LTRbound> <final_lib.fa>\n" unless @ARGV == 3;

open my $MAP,   '<', $ARGV[0] or die "Cannot open $ARGV[0]: $!\n";
open my $BOUND, '<', $ARGV[1] or die "Cannot open $ARGV[1]: $!\n";
open my $LIB,   '<', $ARGV[2] or die "Cannot open $ARGV[2]: $!\n";

# reverse map: original_base (no #class) => new TE ID (kept with #class for output)
my %map;
while (<$MAP>){
	chomp;
	my ($new_name, $ori_name) = split /\t/;
	next unless defined $ori_name;
	(my $ori_base = $ori_name) =~ s/#.*//;
	$map{$ori_base} = $new_name;
}
close $MAP;

# true length of each sequence in the final library, keyed by bare TE ID (class stripped)
my %len;
{
	my $id;
	while (<$LIB>){
		chomp;
		if (/^>(\S+)/){
			($id = $1) =~ s/#.*//;
			$len{$id} = 0;
		} elsif (defined $id){
			s/\s//g;
			$len{$id} += length;
		}
	}
}
close $LIB;

# rewrite IDs and reconcile lengths
my ($kept, $dropped) = (0, 0);
while (<$BOUND>){
	chomp;
	my ($name, $total_len, $lLTR_len, $rLTR_len) = split /\t/;
	next unless defined $rLTR_len;
	(my $base = $name) =~ s/#.*//;
	my $new_name = $map{$base};
	next unless defined $new_name;                 # not renamed into the final library

	(my $new_base = $new_name) =~ s/#.*//;
	my $final_len = $len{$new_base};

	# keep only rows that survived filtering unchanged and are internally consistent
	if (defined $final_len and $final_len > 0
	    and $total_len =~ /^\d+$/ and $lLTR_len =~ /^\d+$/ and $rLTR_len =~ /^\d+$/
	    and $final_len == $total_len
	    and $lLTR_len + $rLTR_len <= $total_len){
		print "$new_name\t$total_len\t$lLTR_len\t$rLTR_len\n";
		$kept++;
	} else {
		$dropped++;
	}
}
close $BOUND;

print STDERR "update_LTRbound.pl: kept $kept boundary entries consistent with the final library, "
	. "dropped $dropped (trimmed or inconsistent; re-measure by self-alignment to recover).\n";
