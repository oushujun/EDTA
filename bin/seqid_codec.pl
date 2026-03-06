#!/usr/bin/env perl
use warnings;
use strict;

# seqid_codec.pl - Encode/decode sequence IDs for the EDTA pipeline
# Shujun Ou (shujun.ou.1@gmail.com)

my $usage = "
seqid_codec.pl - Encode/decode sequence IDs for the EDTA pipeline

Modes:
  encode_fasta <input.fa> <mapfile> <output.fa>
    Clean and optionally encode FASTA seq IDs. Dynamically calculates the
    max allowed ID length based on the longest scaffold size (to fit within
    the 50-char rmblastn limit when coordinates are appended). If any cleaned
    ID exceeds this limit, all IDs are encoded with short base-62 codes.
    If encoding is performed, the mapping file is written; otherwise no
    mapping file is created.

  encode_text <input.txt> <mapfile> <output.txt>
    Replace original seq IDs with encoded codes throughout a text file using
    the mapping file. Useful for converting candidate files (e.g., .scn) to
    use encoded IDs that match an encoded genome.

  decode_fasta <input.fa> <mapfile> <output.fa>
    Restore original seq IDs in FASTA header lines using the mapping file.

  decode_text <input.txt> <mapfile> <output.txt>
    Restore original seq IDs throughout a text file (GFF3, GTF, BED, .out, .sum)
    using the mapping file.

Shujun Ou (shujun.ou.1\@gmail.com)
";

my $mode = shift @ARGV or die $usage;

if ($mode eq 'encode_fasta') {
	my ($infile, $mapfile, $outfile) = @ARGV;
	die "Usage: seqid_codec.pl encode_fasta <input.fa> <mapfile> <output.fa>\n"
		unless defined $infile && defined $mapfile && defined $outfile;
	encode_fasta($infile, $mapfile, $outfile);
} elsif ($mode eq 'encode_text') {
	my ($infile, $mapfile, $outfile) = @ARGV;
	die "Usage: seqid_codec.pl encode_text <input.txt> <mapfile> <output.txt>\n"
		unless defined $infile && defined $mapfile && defined $outfile;
	encode_text($infile, $mapfile, $outfile);
} elsif ($mode eq 'decode_fasta') {
	my ($infile, $mapfile, $outfile) = @ARGV;
	die "Usage: seqid_codec.pl decode_fasta <input.fa> <mapfile> <output.fa>\n"
		unless defined $infile && defined $mapfile && defined $outfile;
	decode_fasta($infile, $mapfile, $outfile);
} elsif ($mode eq 'decode_text') {
	my ($infile, $mapfile, $outfile) = @ARGV;
	die "Usage: seqid_codec.pl decode_text <input.txt> <mapfile> <output.txt>\n"
		unless defined $infile && defined $mapfile && defined $outfile;
	decode_text($infile, $mapfile, $outfile);
} else {
	die "Unknown mode: $mode\n$usage";
}


#################
# encode_fasta
#################
sub encode_fasta {
	my ($infile, $mapfile, $outfile) = @_;

	# Base-62 alphabet: 0-9, a-z, A-Z
	my @b62 = (0..9, 'a'..'z', 'A'..'Z');

	# Read FASTA, clean IDs, measure scaffolds
	my @entries;       # [clean_id, orig_id, seq]
	my %seen_orig;     # detect duplicate original IDs
	my %seen_clean;    # detect duplicate cleaned IDs
	my $max_seq_len = 0;
	my $dup_clean = 0; # flag: cleaning created duplicate IDs

	open my $ifh, '<', $infile or die "ERROR: Cannot open $infile: $!\n";
	local $/ = "\n>";
	while (<$ifh>) {
		chomp;
		s/^>//;
		next if /^\s*$/;
		my ($header, $seq) = (split /\n/, $_, 2);
		next unless defined $header && $header =~ /\S/;
		$seq //= '';
		$seq =~ s/\s+//g;

		# Original ID: first word of the header
		my $orig_id = (split /\s+/, $header)[0];
		next unless defined $orig_id && $orig_id ne '';

		# Check for duplicate original IDs
		die "ERROR: Duplicate sequence ID '$orig_id' in $infile! Please resolve and try again.\n"
			if exists $seen_orig{$orig_id};
		$seen_orig{$orig_id} = 1;

		# Clean: replace special chars with _, collapse, trim
		(my $clean_id = $orig_id) =~ s/[~!@#\$%\^&\*\(\)\+\-=\?\[\]\{\}:;",<\/\\|]+/_/g;
		$clean_id =~ s/_+/_/g;
		$clean_id =~ s/^_//;
		$clean_id =~ s/_$//;

		# Track duplicate cleaned IDs (encoding resolves this)
		$dup_clean = 1 if exists $seen_clean{$clean_id};
		$seen_clean{$clean_id} = 1;

		# Convert non-ATGCN bases to N in sequence
		$seq =~ s/[^ATGCNatgcn]/N/g;

		# Track longest scaffold
		my $slen = length($seq);
		$max_seq_len = $slen if $slen > $max_seq_len;

		push @entries, [$clean_id, $orig_id, $seq];
	}
	close $ifh;
	$/ = "\n";

	die "ERROR: No sequences found in $infile\n" unless @entries;

	# Dynamic id_len_max from longest scaffold size.
	# In LTR_retriever, the composite ID format is chr:s..e|s..e
	#   total_len = len(chr) + 4*D + 6  where D = max coordinate digits
	#   For total <= 50 (rmblastn limit): len(chr) <= 44 - 4*D
	my $max_coord_digits = length("$max_seq_len");
	my $id_len_max = 44 - 4 * $max_coord_digits;

	die "ERROR: Longest scaffold (${max_seq_len} bp) is too large — coordinate digits ($max_coord_digits) " .
		"leave no room for sequence IDs (id_len_max=$id_len_max). Max supported scaffold is ~10 Gb.\n"
		if $id_len_max < 1;

	# Decide: encode if any clean ID is too long or if cleaning created duplicates
	my $need_encoding = $dup_clean;
	unless ($need_encoding) {
		for my $e (@entries) {
			if (length($e->[0]) > $id_len_max) {
				$need_encoding = 1;
				last;
			}
		}
	}

	if ($need_encoding) {
		# Prefix to prevent false-positive matches during decode
		my $PREFIX = '_J';
		my $prefix_len = length($PREFIX);

		# Calculate minimum base-62 digits needed to encode all sequences
		my $n = scalar @entries;
		my $min_digits = 1;
		my $cap = 62;
		while ($cap < $n) { $min_digits++; $cap *= 62; }

		# Use maximum digit width, capped at 9 total chars for LTR_retriever compatibility
		my $max_code_len = 9;
		my $digit_width = (($id_len_max < $max_code_len) ? $id_len_max : $max_code_len) - $prefix_len;

		die "ERROR: Too many sequences ($n) for id_len_max=$id_len_max " .
			"(need $min_digits base-62 digits + ${prefix_len}-char prefix, " .
			"but only $digit_width digits available)\n"
			if $digit_width < $min_digits;

		# Write mapping file and encoded FASTA
		open my $mfh, '>', $mapfile or die "ERROR: Cannot write $mapfile: $!\n";
		open my $ofh, '>', $outfile or die "ERROR: Cannot write $outfile: $!\n";

		for my $i (0 .. $#entries) {
			my $code = $PREFIX . _int_to_base62($i, $digit_width, \@b62);
			print $mfh "$code\t$entries[$i][1]\n";
			print $ofh ">$code\n$entries[$i][2]\n";
		}

		close $mfh;
		close $ofh;

		my $max_clean_len = 0;
		for my $e (@entries) {
			my $l = length($e->[0]);
			$max_clean_len = $l if $l > $max_clean_len;
		}
		my $total_code_len = $prefix_len + $digit_width;
		print STDERR "Encoded $n sequences with ${total_code_len}-char codes " .
			"(prefix '${PREFIX}' + ${digit_width} base-62 digits, " .
			"longest clean ID: $max_clean_len, id_len_max: $id_len_max, " .
			"longest scaffold: ${max_seq_len} bp)\n";
	} else {
		# No encoding needed — write cleaned FASTA, no mapping file
		open my $ofh, '>', $outfile or die "ERROR: Cannot write $outfile: $!\n";
		for my $e (@entries) {
			print $ofh ">$e->[0]\n$e->[2]\n";
		}
		close $ofh;

		# Remove stale mapping file from a previous run if present
		unlink $mapfile if -e $mapfile;

		print STDERR "No encoding needed — all " . scalar(@entries) . " IDs fit within " .
			"id_len_max=$id_len_max (longest scaffold: ${max_seq_len} bp)\n";
	}
}


# Convert an integer to a zero-padded base-62 string
sub _int_to_base62 {
	my ($num, $width, $alphabet) = @_;
	my @digits;
	for (1 .. $width) {
		unshift @digits, $alphabet->[$num % 62];
		$num = int($num / 62);
	}
	return join '', @digits;
}


#################
# encode_text
#################
sub encode_text {
	my ($infile, $mapfile, $outfile) = @_;
	my ($map_ref, $prefix, $digit_width) = _load_map($mapfile);

	# Build reverse map: orig_id → code
	my %rmap;
	for my $code (keys %$map_ref) {
		$rmap{$map_ref->{$code}} = $code;
	}

	open my $ifh, '<', $infile or die "ERROR: Cannot open $infile: $!\n";
	open my $ofh, '>', $outfile or die "ERROR: Cannot write $outfile: $!\n";

	while (<$ifh>) {
		chomp;
		my @parts = split /(\s+)/;  # preserve delimiters
		for my $part (@parts) {
			$part = $rmap{$part} if exists $rmap{$part};
		}
		print $ofh join('', @parts), "\n";
	}

	close $ifh;
	close $ofh;
}


#################
# decode_fasta
#################
sub decode_fasta {
	my ($infile, $mapfile, $outfile) = @_;
	my ($map_ref, $prefix, $digit_width) = _load_map($mapfile);

	my $regex;
	if ($prefix ne '') {
		$regex = qr/(?<![0-9a-zA-Z_])(\Q$prefix\E[0-9a-zA-Z]{$digit_width})(?![0-9a-zA-Z_])/;
	} else {
		$regex = qr/(?<![0-9a-zA-Z])([0-9a-zA-Z]{$digit_width})(?![0-9a-zA-Z])/;
	}

	open my $ifh, '<', $infile or die "ERROR: Cannot open $infile: $!\n";
	open my $ofh, '>', $outfile or die "ERROR: Cannot write $outfile: $!\n";

	while (<$ifh>) {
		if (/^>/) {
			s/$regex/exists $map_ref->{$1} ? $map_ref->{$1} : $1/ge;
		}
		print $ofh $_;
	}

	close $ifh;
	close $ofh;
}


#################
# decode_text
#################
sub decode_text {
	my ($infile, $mapfile, $outfile) = @_;
	my ($map_ref, $prefix, $digit_width) = _load_map($mapfile);

	my $regex;
	if ($prefix ne '') {
		$regex = qr/(?<![0-9a-zA-Z_])(\Q$prefix\E[0-9a-zA-Z]{$digit_width})(?![0-9a-zA-Z_])/;
	} else {
		$regex = qr/(?<![0-9a-zA-Z])([0-9a-zA-Z]{$digit_width})(?![0-9a-zA-Z])/;
	}

	open my $ifh, '<', $infile or die "ERROR: Cannot open $infile: $!\n";
	open my $ofh, '>', $outfile or die "ERROR: Cannot write $outfile: $!\n";

	while (<$ifh>) {
		s/$regex/exists $map_ref->{$1} ? $map_ref->{$1} : $1/ge;
		print $ofh $_;
	}

	close $ifh;
	close $ofh;
}


# Load the mapping file; returns (hashref, prefix, digit_width)
sub _load_map {
	my ($mapfile) = @_;
	my %map;
	my $prefix;
	my $digit_width;

	open my $mfh, '<', $mapfile or die "ERROR: Cannot open map file $mapfile: $!\n";
	while (<$mfh>) {
		chomp;
		my ($code, $orig) = split /\t/, $_, 2;
		next unless defined $code && defined $orig;
		$map{$code} = $orig;
		unless (defined $digit_width) {
			if ($code =~ /^(_J)/) {
				$prefix = $1;
				$digit_width = length($code) - length($prefix);
			} else {
				$prefix = '';
				$digit_width = length($code);
			}
		}
	}
	close $mfh;

	die "ERROR: Empty mapping file $mapfile\n" unless defined $digit_width;
	return (\%map, $prefix, $digit_width);
}
