#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Convert EDTA GFF3 file to TEtranscripts GTF file
# GFF to GTF Mapping:
#   ID   -> transcript_id
#   Name -> gene_id
#   3rd column (type) -> family_id
#   feature -> exon (constant)
#   classification -> class_id
#
# Filtering:
#   Remove records whose 3rd column matches any in --remove list
#   Default remove list: repeat_region,long_terminal_repeat,target_site_duplication
#
# Usage examples:
#   perl gff_to_gtf.pl --gff in.gff3 --out out.gtf
#   perl gff_to_gtf.pl --gff in.gff3 > out.gtf
#   cat in.gff3 | perl gff_to_gtf.pl --gff - > out.gtf
#   perl gff_to_gtf.pl --gff - --remove knob,repeat_region < in.gff3 > out.gtf
# Shujun Ou (shujun.ou.1@gmail.com) and ChatGPT 5.2
# v1: 01/04/2025

my $gff_in  = undef;   # file path or '-' for STDIN
my $out_fn  = undef;   # file path or undef for STDOUT
my @remove  = ();      # user-specified remove list

GetOptions(
    'gff=s'    => \$gff_in,
    'out=s'    => \$out_fn,
    'remove=s' => \@remove,   # can be passed multiple times; each may be comma-separated
    'help'     => sub { usage(0) },
) or usage(1);

if (!@ARGV && !defined $gff_in && !defined $out_fn && !@remove) {
    usage(0);
}

# Default input is STDIN if not provided
$gff_in = '-' if !defined $gff_in;

# Build removal set:
# If --remove is provided, use ONLY that list.
# Otherwise, use defaults.
my @remove_list;
if (@remove) {
    for my $r (@remove) {
        push @remove_list, grep { length($_) } split /,/, $r;
    }
} else {
    @remove_list = qw(repeat_region long_terminal_repeat target_site_duplication);
}

my %REMOVE = map { $_ => 1 } @remove_list;

# Input handle
my $IN;
if ($gff_in eq '-') {
    $IN = *STDIN;
} else {
    open($IN, '<', $gff_in) or die "Cannot open --gff '$gff_in': $!\n";
}

# Output handle
my $OUT;
if (defined $out_fn && $out_fn ne '-') {
    open($OUT, '>', $out_fn) or die "Cannot open --out '$out_fn': $!\n";
} else {
    $OUT = *STDOUT;
}

while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;

    my @f = split /\t/, $line, -1;
    if (@f < 9) {
        warn "Skipping (not 9 columns): $line\n";
        next;
    }

    my ($seqname, $source, $type, $start, $end, $score, $strand, $phase, $attr) = @f;

    # Filter by 3rd column
    next if exists $REMOVE{$type};

    # Parse GFF3 attributes: key=value;key2=value2
    my %a;
    for my $kv (split /;/, $attr) {
        next if $kv eq '';
        my ($k, $v) = split /=/, $kv, 2;
        next unless defined $k && defined $v;
        $k =~ s/^\s+|\s+$//g;
        $v =~ s/^\s+|\s+$//g;
        $a{$k} = $v;
    }

    # force + strand for missing values
    #    if (!defined $strand || $strand eq '' || $strand eq '.' || $strand eq '?') {
    #	 $strand = '+';
    #	 }

    my $id   = $a{ID};
    my $name = $a{Name};
    my $class = $a{classification} || $a{Classification};
    $class =~ s/\/.*//;

    if (!defined $id || $id eq '') {
        warn "Skipping (missing ID=): $line\n";
        next;
    }
    if (!defined $name || $name eq '') {
        warn "Skipping (missing Name=): $line\n";
        next;
    }
    if (!defined $class || $class eq '') {
        warn "Skipping (missing classification=): $line\n";
        next;
    }

    # GTF: feature is always exon; keep score from GFF
    my $feature = "exon";

    my $gtf_attr = join " ",
        qq(gene_id "$name";),
        qq(transcript_id "$id";),
	qq(family_id "$name";), # $name = family
	#qq(family_id "$type";), # $type = superfamily
        qq(class_id "$class";);

    print $OUT join("\t",
        $seqname, $source, $feature, $start, $end, $score, $strand, $phase, $gtf_attr
    ), "\n";
}

close $IN  if ($gff_in ne '-');
close $OUT if (defined $out_fn && $out_fn ne '-');

exit 0;

sub usage {
    my ($exit_code) = @_;
    print STDERR <<"USAGE";
Usage:
  perl gff_to_gtf.pl --gff <in.gff3|-> [--out <out.gtf|->] [--remove a,b,c] [--remove x,y]

Notes:
  - Use --gff - to read from STDIN (pipe-friendly).
  - Omit --out (or use --out -) to write to STDOUT.
  - Keeps the original GFF score column (6th column).
  - Filters out records whose 3rd column (type) matches --remove list.
  - Default removed types (when --remove not provided):
      repeat_region,long_terminal_repeat,target_site_duplication

Examples:
  perl gff_to_gtf.pl --gff in.gff3 > out.gtf
  cat in.gff3 | perl gff_to_gtf.pl --gff - --out out.gtf
  perl gff_to_gtf.pl --gff in.gff3 --remove knob,repeat_region > out.gtf
USAGE
    exit($exit_code);
}

