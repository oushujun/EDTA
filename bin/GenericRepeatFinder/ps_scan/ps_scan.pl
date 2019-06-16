#!/usr/bin/env perl

# ps_scan - a PROSITE scanning program
#
# Revision: 1.86
#
# Copyright (C) 2001-2011 Swiss Institute of Bioinformatics
# Authors:
#   edouard.decastro@isb-sib.ch
#   Alexandre Gattiker
#   Beatrice Cuche (evaluated_by post-processing)
#   Lorenzo Cerutti (FTREP method)
#
# Contributions:
#   Lorenza Bordoli (repeat method)

# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place - Suite 330, Boston, MA  02111-1307, USA.


use 5.005_03; # uses perl >= v5.5
use IO::File;
use Carp qw(confess cluck);
use vars qw(@ISA $VERSION $errpos $errstr);
use strict;


################################################################################
# integrated subs taken from Prosite.pm module:


# scan a sequence with a perl pattern
sub scanPattern {
    my ($pattern, $sequence, $behavior, $max_x, $opt_miniprofiles) = @_;
    $behavior ||= 0;
    my $allowOverlap = !($behavior & 1);
    my $allowInclude = $behavior & 2 ;
    $max_x ||= 0;
    my @hits;
    my $pos = 0;
    my @comb = $pattern =~ /\((.*?)\)/g;
    if ($pattern) {
        my $prevstop = -1;
        my @tok;
        my $nter_anchor = $pattern=~/^\^/ ? 1 : 0;
        while (@tok = (substr($sequence, $pos) =~ /^(.*?)($pattern)/)) {
            my $prematch = shift @tok;
            my $subseq = shift @tok;
            my $length = length $subseq;
            my $number_x = 0;
            if (@tok == @comb and @tok) {
                $subseq = "";
                for (my $i = 0; $i < @tok; $i++) {
                    if ($comb[$i] =~ /\./) {
                        $tok[$i] =~ tr/A-Z/a-z/;
                    }
                    # don't count X's that match pattern elements which exclude
                    # certain AA's
                    elsif ($comb[$i] =~ /^\[\^/) {
                    }
                    # count number of X's matching a non-x position in the
                    # pattern
                    elsif (my $x_count = $tok[$i] =~ tr/Xx/Xx/) {
                        $number_x += $x_count;
                    }
                    my $tok = $comb[$i] =~ /\./ ? lc($tok[$i]) : $tok[$i];
                    $subseq .= $tok;
                    if (my @numbers = $comb[$i] =~ /(\d+)/g) {# insertion
                        my $biggest = pop @numbers;
                        $subseq .= "." x ($biggest - length $tok);
                    }
                }
            }
            elsif (@tok != @comb) {
                cluck "Internal error with regular expression $pattern\n";
            }
            my $shift = ($length || 1) - 1;
            my $stop = $pos + $length + length $prematch;
            $pos = $stop;
            $pos -= $shift if $allowOverlap;
            # may happen if empty pattern match
            last if $pos > length $sequence;
            if ($length) {
                if ($allowInclude or $stop > $prevstop) {
                    if ($number_x <= $max_x or $max_x<0) {
                        push @hits, [$subseq, $stop - $length + 1,
                            $stop, undef,
                            ($opt_miniprofiles && !@hits ? $sequence : undef) ];
                                # also push sequence for ps_scan option -b
                        $prevstop = $stop;
                    }
                }
            }
            else {
                $pos++;# empty pattern match
            }
            last if $nter_anchor and $pos;
        }
    }
    return \@hits;
}

sub scanProfiles {
    my ($file, $level_min) = @_;
    $level_min = -99 unless defined $level_min;
    my $results = {};
    local $/ = "\n";
    my $must_open_file = !UNIVERSAL::isa($file, "GLOB");
    my $pfscan_h;
    if ($must_open_file) {
        $pfscan_h = new IO::File $file or confess "Could not open $file : $!";
    }
    else {
        $pfscan_h = $file;
    }
    my $lastac = "";
    my $last_entry;
    while (defined (local $_=<$pfscan_h>)) {
        if (my ($id1, $level, $levelna, $nscore, $rawscore, $from, $to,
            $pffrom, $pfto, $repregion, $repnumber, $ac, $de) = m/
        (
        >\S*
        )? #ID of the profile, if option -x or -s
        (?<!L=\d) #These three lines are to fix a bug in the output of pfscan if >99 matches are found: the id and the level are then pasted together with an intervening "1" (e.g. Q9LGZ9 vs. preprofile PS50079). FIXME: This still doesn't work with pfsearch -xlr
        (?<!L=[-\d]\d)
        (?<!L=[-\d]\d\d)
        \s*
        (?:L=(\S+)|(NA))? #level, if option -l
        \s*
        (?:(\S+)\s+)? #normalized score, unless option -r
        (\S+) #raw score
        \s+
        pos\.
        \s*
        (\d+) #seq from
        \s*-\s*
        (\d+) #seq to
        \s*
        (?:
        \[\s*
        (\d+) #profile from
        ,\s*
        (-\d+) #profile to
        \]
        \s*
        )? #if option -z
        (REGION\d+\s)? #repeat region, if option -m
        (\d+REP\d+\s)? #repeat number, if option -m
        (\S+) #profile or sequence AC and-or ID
        (?:\s*(.+))? #sequence description when running pfsearch, may be absent
            /x) {
            # fix bug in pfsearch/pfscan which report "******" if nscore>999.999
            $nscore = "999.999" if $nscore eq "*******";
            $level = $level_min if !defined($level) and defined $levelna;
            #            0   1      2    3    4        5      6          7        8       9      10
            my $entry = ["", $from, $to, $ac, $pffrom, $pfto, $rawscore, $nscore, $level, undef, $de, []];
            if ($repnumber) { push @{$results->{$lastac}->[-1]->[11]}, $entry }
            else { push @{$results->{$ac}}, $entry }
            $lastac = $ac;
            $last_entry = $entry;
        }
        else {
            # sequence
            next unless $last_entry;
            $last_entry->[0] .= $_;
        }
    }
    if ($must_open_file) {
        close $pfscan_h or confess "Error $? closing $file";
    }
    return $results;
}


sub prositeToRegexp_old {
    local $_ = shift;
    my $notGreedy = shift;
    my $ungreed = $notGreedy ? "?" : "";
    s/\.$//;# possible end dot from parsing
    s/-//g;
    1 while s/\{([^\}^\^]*)([^\^^\}][^\}]*\})/\{^$1$2/; #insert ^
    s/\}/]/g;
    s/\{/[/g;
    s/\(/{/g;
    s/\)/}$ungreed/g;
    s/x/./ig;
    s/B/[ND]/g;
    s/Z/[QE]/g;
    # PS00539 is P-R-L-[G>] which would be converted to PRL[G$], but "$" is
    # not valid in a character range, we have to convert this to (?:[G]|$)
    s/\[([^\[\]]*)([<>])([^\[\]]*)\]/(?:[$1$3]|$2)/g;
    # do this after the previous step, so as not to mix ^ for excluding
    # character classes with ^ to indicate anchoring
    s/</^/g; # anchors
    s/>/\$/g;
    # insert () around match states and insertions
    s/ (\[[^\]]*\]|[\w.]) ( \{ \d+(,\d+)? \} )/)($1$2)(/xg;
    return "($_)";
}

# same, using tokenizing parser
sub prositeToRegexp {
    my $string = shift;
    my $notGreedy = shift;
    my $ungreed = $notGreedy ? "?" : "";
    my $preventX = shift;
    $errstr = undef;
    $errpos = undef;
    my $pushback = "";
    my $ntok = 0;
    my $regexp = "";
    my $get = sub {
        $ntok++;
        return ($pushback =~ s/(.)// ? $1 : ($string =~ s/(.)// ? $1 : undef));
    };
    while (defined (my $tok = &$get)) {
        my $state;
        my $not;
        if ($tok eq "-") {
        # ignore
        }
        elsif ($tok eq "[") {# RANGE
            while(defined (my $tok = &$get)) {
                last if $tok eq "]";
                $state .= $tok;
            }
        }
        elsif ($tok eq "{") {# negative RANGE
            $not = 1;
            while(defined(my $tok = &$get)) {
                last if $tok eq "}";
                $state .= $tok;
            }
        }
        elsif ($tok =~ /[A-Za-z]/) {# single char
            $state = $tok;
        }
        elsif ($tok eq "<") {
            $regexp .= "^";
        }
        elsif ($tok eq ">") {
            $regexp .= '$';
        }
        else {
            $errstr = "Parsing error";
            $errpos = $ntok-1;
            return undef;
        }
        if (defined $state) {
        # read range, e.g. "x(2,5)"
            my $range;
            my $range_char;
            if (defined(my $tok = &$get)) {
                if ($tok eq "(") {
                    while (defined(my $tok = &$get)) {
                        last if $tok eq ")";
                        $range .= $tok;
                    }
                }
                elsif ($tok eq "*") {# support e.g. "<{C}*>"
                    $range_char = $tok;
                }
                else {
                    $pushback .= $tok;
                    $ntok--;
                }
            }

            if ($state =~ /x/i) {
                $state = ".";
            }
            else {
            # handle B/Z unsure amino acids both in pattern and sequence
                if ($not) {
                    $state =~ s/B/NDB/g;
                    $state =~ s/Z/QEZ/g;
                }
                else {
                    $state =~ s/B/NDB/g or $state =~ s/([ND])/$1B/g;
                    $state =~ s/Z/QEZ/g or $state =~ s/([QE])/$1Z/g;
                    $state .= "X" unless $preventX;
                }
            }
            my $mod = $1 if $state =~ s/([<>])//g;
            # "$" is not valid in a character range, we have to convert this
            # to (?:[GH]|$)
            $regexp .= "(";
            $regexp .= "(?:" if $mod;
            $regexp .= "[" if length($state) > 1 or $not;
            $regexp .= "^" if $not;
            $regexp .= $state;
            $regexp .= "]" if length($state) > 1 or $not;
            $regexp .= "|" . ($mod eq "<" ? "^" : '$') . ")" if $mod;
            $regexp .= "{$range}$ungreed" if defined $range;
            $regexp .= "$range_char" if defined $range_char;
            $regexp .= ")";
        }
    }
    return $regexp;
}


# Checks that a user-entered pattern is parseable, returning an error message
# or undef TODO: <A-T-[<GE] should be an error
sub checkPatternSyntax {
    my $pattern = shift;
    my $c1 = 0;
    my $c2 = 0;
    my ($c_open_square, $c_open_curly, $c_open_paren) = (0, 0, 0);
    my ($c_close_square, $c_close_curly, $c_close_paren) = (0, 0, 0);
    if ($pattern =~ /(-){2,}/ || $pattern =~ /([,\-\(\)\{\}\[\]\<\>]){2,}]/) {
        return "duplicate character \"$1\"";
    }
    if ($pattern !~ /[a-zA-Z]/) {
        return "pattern has no characters";
    }
    if ($pattern =~ /-\(/) {
        return "dash before (";
    }
    if ($pattern =~ /([JOU])/i) {
        return "pattern contains letter \"$1\" which is not an amino acid";
    }
    if (length($pattern) > 200) {
        return "pattern is longer than the limit of 200 characters";
    }
    elsif ($pattern =~ /^\[[a-z]+\]$/i) {
        return "pattern is too degenerate";
    }
    else {
        my $ambig;
        my $ambig_complement;
        my $range;
        my %count;
        my @ambig;
        foreach (split(//,$pattern)) {
            unless ($c1 == $c2 || $c1 == $c2 +1) {
            # always close parentheses before opening a new one!
                return "nested parentheses are forbidden";
            }
            if (/[\[\{\(]/) {
                $c1++;
                if (/\[/) {
                    $ambig = " ";
                    $c_open_square++;
                }
                elsif (/\{/) {
                    $ambig_complement = " ";
                    $c_open_curly++;
                }
                elsif (/\(/) {
                    $range = " ";
                    $c_open_paren++;
                }
            }
            elsif (/[\]\}\)]/) {
                $c2++;
                if (/\]/) {
                    %count = ();
                    if (length($ambig) < 3) {
                        return "no real ambiguity inside []";
                    }
                    else {
                        @ambig = split(//, $ambig);
                        for (@ambig) {
                            $count{$_}++;
                        }
                        for (sort keys %count) {
                            if ($count{$_} ne 1) {
                            return "string inside square brackets \"
                                $ambig\" contains duplicates";
                            }
                        }
                    }
                    $ambig = "";
                    $c_close_square++;
                }
                elsif (/\}/) {
                    $ambig_complement = "";
                    $c_close_curly++;
                }
                elsif (/\)/) {
                    $range =~ s/^\s+//;
                    if ($range =~ /^(\d+),(\d+)$/) {
                        if ($1 >= $2) {
                            return "range \"$range\" is invalid
                                (second term must be greater than first)";
                        }
                    }
                    elsif ($range !~ /^\d+$/){
                        return "range \"$range\" is invalid";
                    }
                    $range = "";
                    $c_close_paren++;
                }
            }
            elsif ($range) {
                if (!/[\d,]/) {
                    return "incorrect range \"$range\"";
                }
                $range .= $_;
            }
            elsif ($ambig) {
                if (!/[A-Z<>]/) {# [G>] is allowed, e.g. PS00267, PS00539
                    return "wrong syntax for ambiguity : \"$ambig \"";
                }
                if (/([BZ])/) {
                    return "ambiguous amino acid \"
                        $1\" not allowed within ambiguity";
                }
                $ambig .= $_;
            }
            elsif ($ambig_complement) {
                if (!/[A-Z]/) {
                    return "wrong syntax for ambiguity :
                        \"$ambig_complement\"";
                }
                $ambig_complement .= $_;
            }
            else {# amino acid or anchor, or * quantifier
                if (/([^A-Zx\-<>*])/) {
                    return "invalid character : \"$1\"";
                }
            }
        }
        if ($c1 != $c2) {
            return "unbalanced (), [] or {}";
        }
        elsif ($c_open_square != $c_close_square) {
            return "unbalanced []";
        }
        elsif ($c_open_curly != $c_close_curly) {
            return "unbalanced {}";
        }
        elsif ($c_open_paren != $c_close_paren) {
            return "unbalanced ()";
        }
    }
    return undef;
}

sub parseProsite {
    local $_ = shift;
    my $ac = $1 if /^AC   (\w+)/m;
    my ($id, $type) = ($1, $2) if /^ID   (\w+); (\w+)/m;
    my $de_line = join "\n", /^(DE.*\S)/mg;
    my $de = join " ", $de_line =~ /^DE   (.*\S)/mg;
    my $pa = join "", /^PA   (.*\S)/mg;
    $pa =~ s/\.$//;
    my $rule = join "\\\n", /^RU   (.*\S)/mg;# obsolete
    my $skip = /^CC   \/SKIP-FLAG=TRUE/m;
    my $pdoc = $1 if /^DO   (\w+)/m;
    my $tax = $1 if /^CC   \/TAXO-RANGE=(.*?);/m;
    my $rep = $1 if /^CC.*\/MAX-REPEAT=(\d+)/m;
    my @sites;push (@sites, [$1, $2])
        while(m/\/SITE=(\d+),(.*?);/g);
    my @cutoffs = map {my %a=map {split /=/, $_, 2}
        split /; */, $_; \%a } /^MA   \/CUT_OFF: (.*)/mg;
    my @pps;push (@pps, {'type'=>$1,'discriminators'=>[split(';\s*',$2)]})
        while(m/^PP   \/?(\w+):\s*([^\r\n]+)/mg);
    return ($ac, $id, $type, $de, $pa, $rule,
        $pdoc, $skip, $tax, $rep, \@sites, \@cutoffs, \@pps);
};



################################################################################
# integrated subs taken from FTRep.pm module:

use constant PI           => 3.14159265;
use constant PI2          => 2*PI;

#-------------------------------------------------------------#
# sinus cardinal
sub sinc
{
    my $x = shift;
    return ($x == 0)? 1 : sin($x)/$x;
}

#-------------------------------------------------------------#
# Fourier transform
sub ft
{
    my ($x,$score,$d,$size) = @_;
    return 0 if ($size < 1);
    my $t  = PI2*($x-$d)/$size;
    my $sint = sinc($t);
    my $sint2 = $sint*$sint;
    return $score * $sint2;
}

sub FT
{
    my ($matchlist,$groups,$motif_param,$amplitude) = @_;

    $amplitude = 0; # Turn it off

    # ft and autocorrelation
    my $ft = []; # for plot
    my $fa = []; # for plot and integral

    for (my $i = 0; $i < scalar(@$groups); $i++)
    {
        for (my $x = $groups->[$i]->{start}; $x <= $groups->[$i]->{end}; $x++) # walk around the group region in the sequence
        {
            my $y  = 0; # store ft function result
            my $yk = 0; # store ft with lag function
            my $n_members = scalar(@{$groups->[$i]->{members}});
            if ($n_members > 1)
            {
                for (my $g = 0 ; $g < $n_members; $g++)
                {
                    my $k = ($g == 0)? 1 : -1; # deal with left border g-1,g becomes g,g+1

                    my $m = $groups->[$i]->{members}->[$g];
                    my $n = $groups->[$i]->{members}->[$g+$k];
                    my $m_amp = ($amplitude)? $amplitude : $matchlist->[$m]->{match_rscore}/$motif_param->{maxscore};
                    my $n_amp = ($amplitude)? $amplitude : $matchlist->[$n]->{match_rscore}/$motif_param->{maxscore};

                    $y  += ft($x,$m_amp,$matchlist->[$m]->{match_pos},$matchlist->[$m]->{match_size});
                    $yk += ft($x+$motif_param->{motif_len}*$k,$n_amp,$matchlist->[$n]->{match_pos},$matchlist->[$n]->{match_size});
                }
            }
            $ft->[$i][$x] += $y;
            $fa->[$i][$x] += $y * $yk;
        }
    }
    return ($ft,$fa);
}

#-------------------------------------------------------------#
# autocorrelation
sub AC
{
    my ($matchlist,$groups,$fa) = @_;
    my $ac = [];

    # integral of autocorrelation for each group
    for (my $i = 0; $i < scalar(@$groups); $i++)
    {
        $ac->[$i] = 0;
        for (my $x = $groups->[$i]->{start}+1; $x <= $groups->[$i]->{end}; $x++)
        {
            $ac->[$i] += ($fa->[$i][$x] + $fa->[$i][$x-1]) / 2;
        }
    }
    return $ac;
}

#-------------------------------------------------------------#
# select match
sub select_match
{
    my ($matchlist,$groups,$ac,$cutoff,$filter) = @_;
    my $max_prot_score = -999999;

    for (my $i = 0; $i < scalar(@$groups); $i++)
    {
        foreach my $g (@{$groups->[$i]->{members}})
        {
            if ( $cutoff < $ac->[$i] )
            {
                $matchlist->[$g]->{selected} = 1;
            }

            if ($matchlist->[$g]->{match_nscore} > $max_prot_score)
            {
                $max_prot_score = $matchlist->[$g]->{match_nscore}
            }
        }
    }

    if (defined($filter) && $max_prot_score < $filter)
    {
        for (my $i = 0; $i < scalar(@$groups); $i++)
        {
            foreach my $g (@{$groups->[$i]->{members}})
            {
                $matchlist->[$g]->{selected} = 0;
            }
        }
    }

    return (1);
}

#-------------------------------------------------------------#
# Sort positions
sub sortlist
{
    my ($matchlist) = shift;
    my %sortindex;
    for (my $i = 0; $i < scalar(@$matchlist); $i++)
    {
        $sortindex{$matchlist->[$i]->{match_start}} = $i;
    }

    my @sortindex;
    foreach (sort {$a<=>$b} keys %sortindex)
    {
        push(@sortindex,$sortindex{$_});
    }

    return @sortindex;
}

#-------------------------------------------------------------#
# Group matches
sub group_matches
{
    my ($matchlist,$grouping_size) = @_;

    # sort matches by their position
    my @index = sortlist($matchlist);

    my $keep_prot = 0; # used for filter
    my ($start,$end);
    my $group_index = 0;
    my $groups = [];

    # First match
    push(@{$groups->[$group_index]->{members}},$index[0]);
    $start = $matchlist->[$index[0]]->{match_start};

    for (my $i = 1; $i < scalar(@index); $i++)
    {
        if ($matchlist->[$index[$i]]->{match_pos} - $matchlist->[$index[$i-1]]->{match_pos} > $grouping_size)
        {
            # Set limits for current group
            $end = $matchlist->[$index[$i-1]]->{match_end};
            $groups->[$group_index]{start} = $start;
            $groups->[$group_index]{end}   = $end;
            $start = $matchlist->[$index[$i]]->{match_start};

            # Create new group
            $group_index++;
        }
        elsif (!$keep_prot)
        {
            $keep_prot = 1; # keep the protein if at least 2 match have distance < grouping_size
        }

        $groups->[$group_index]->{match_index} = 0;
        push(@{$groups->[$group_index]->{members}},$index[$i]);
    }

    # Deal with last match
    $end = $matchlist->[$index[scalar(@index)-1]]->{match_end};
    $groups->[$group_index]->{start} = $start;
    $groups->[$group_index]->{end}   = $end;

    if (!$keep_prot)
    {
        #$groups = []
    }
    return $groups;
}




################################################################################
# ps_scan_source.pl core:


# initializations & parameters processing

BEGIN {
   $VERSION = '1.86';
}

# Can we use the IPC::Open2 module to communicate with
# pfscan via pipes (instead of temp files) ?S
eval { require IPC::Open2 };
my $NO_DIRECT_PIPE=$? || $^O eq "MSWin32" ? 1 : 0;
# change this to the absolute path to the programs,
# unless they are located in a directory in your PATH
# or use option --pfscan and/or --psa2msa to give the full path.
my $PFSCAN  = 'pfscan';
my $PSA2MSA = 'psa2msa';
my $errcode = 0;
my $MOTIF_AC_REGEXP = '\w+\d+';#'PS\d{5}';
$|= 1;
use Getopt::Long;
use Data::Dumper;
use IO::File;
my @formats = qw(scan fasta psa msa gff pff epff sequence matchlist ipro);
my $formats_string = join " | ", @formats;

sub usage {
    my $progname =  $0;
    $progname    =~ s/.*[\\\/]//;
    print <<EOF;
$progname [options] sequence-file(s)
ps_scan version $VERSION options:
-h : this help screen

Input/Output:
  -e <string> : specify the ID or AC of an entry in sequence-file
  -o <string> : specify an output format : scan | fasta | psa | msa |
                gff | pff | epff | sequence | matchlist
  -d <file>   : specify a prosite motif file
  -p <string> : specify a pattern or the AC of a prosite motif
  -f <string> : specify a motif AC to scan against together with all its
                related post-processing motifs (but show only specified
                motif hits)

Selection:
  -r          : do not scan profiles
  -m          : only scan profiles
  -s          : skip frequently matching (unspecific) patterns and
                profiles
  -l <number> : cut-off level for profiles (default : 0)

Pattern option:
  -x <number> : specify maximum number of accepted matches of X's in
                sequence (default=0)
  -g          : Turn greediness off
  -v          : Turn overlaps off
  -i          : Allow included matches
  -b <file>   : use profiles from <f> to evaluate pattern matches
                (a LevelTag is assigned to each pattern match). If no
                file is specified, evaluator.dat will be used (searched
                in the paths \$PROSITE/ and \$PROSITE/prosite/)

profile options:
  -w pfsearch : Compares a query profile against a protein sequence
                library

Other expert options:
  --nopp             : do not post-process matches
  --reverse          : randomize the sequence database by taking the
                       reverse sequence of each individual entry
  --shuffle          : randomize the sequence database by local shuffling
                       of the residues in windows of 20 residues
  --gff <file>       : use an existing ps_scan gff result <file> as input to
                       e.g. post-process it and/or to reformat it into another
                       format (defined by -o <string>), instead of doing a
                       real scan.
  -w pfsearch -R     : use raw scores rather than normalized scores for
                       match selection
  -w pfsearch -C <x> : report only match scores higher than the specified
                       parameter <x>; an integer argument is interpreted
                       as a raw score value, a decimal argument as a
                       normalized score value. -R and -C options can be
                       combined.

Note:
  The sequence-file may be in Swiss-Prot or FASTA format.

  If no prosite motif file is specified, prosite.dat will be used
  (searched in the paths \$PROSITE/ and \$SPROT/prosite/).

  There may be several -d, -p and -e arguments.
EOF
    exit 1;
}

my $opt_noprofiles;
my $opt_onlyprofiles;
my $opt_skip;
my $opt_skiponly;
my $opt_help;
my $opt_format;
my $opt_max_x;
my $opt_nongreedy;
my $opt_nooverlaps;
my $opt_miniprofiles;
my $opt_includes;
my $opt_level = 0;
my $opt_rep_pp_4allprofiles;
my $opt_pfsearch;
my $opt_cutoff;
my $opt_raw;
my $opt_minhits;
my $opt_maxhits;
my $opt_filterheader;
my $opt_reverse;
my $opt_shuffle;
my $opt_no_postprocessing;


# list of prosite 'dat' files to be scanned
my @prosite_files;
# !can contain a list of user pattern (and/or)
# list of prosite AC (profile/pattern)!
my @motifAC_or_userpattern;
# list of ID/AC
my @entries;
# list of motif AC for which one wants to retrieve (and scan)
# all other motif linked by post-processing
my @followpp;
# list of ps_scan gff result files to be scanned
my @external_gff_files;

my $SLASH  = $^O eq "MSWin32" ? "\\" : "\/";
my $TMPDIR = ".";
for my $dir ( $ENV{TMPDIR}, $ENV{SP_TEMP}, $ENV{TMP}, $ENV{TEMP},
    "/tmp", "c:\\temp", "c:\\tmp" ) {
    if ( defined($dir) and -d $dir ) {
        $TMPDIR = $dir;
        last;
    }
}
my $TMP_COUNTER = 1;
my $TMP_BASE    = int( rand( 1000000 ) ) + int( substr( abs( $$ ), -6 ) );

my $scan_profiles;
my $scan_pattern;

my $last_profile_tmp_filename;

Getopt::Long::Configure ("bundling", "no_ignorecase");
GetOptions (
    "r"   => \$opt_noprofiles,
    "m"   => \$opt_onlyprofiles,
    "s"   => \$opt_skip,
    "h"   => \$opt_help,
    "v"   => \$opt_nooverlaps,
    "i"   => \$opt_includes,
    "g"   => \$opt_nongreedy,
    "b:s" => \$opt_miniprofiles,
    "x=i" => \$opt_max_x,
    "l=i" => \$opt_level,
    "o=s" => \$opt_format,
    "d=s" => \@prosite_files,
    "p=s" => \@motifAC_or_userpattern,
    "e=s" => \@entries,
    "f=s" => \@followpp,
    # a: undocumented, just for PROSITE team to test repeat post-processing
    # on any normal ('!' level0) profile...
    "a"   => \$opt_rep_pp_4allprofiles,
    "w=s" => \$opt_pfsearch,
    "C=f" => \$opt_cutoff,
    "R"   => \$opt_raw,
    "pfscan=s"       => \$PFSCAN,
    "psa2msa=s"      => \$PSA2MSA,
    "minhits=i"      => \$opt_minhits,
    "maxhits=i"      => \$opt_maxhits,
    "filterheader=s" => \$opt_filterheader,
    "reverse"        => \$opt_reverse,
    "shuffle=i"      => \$opt_shuffle,
    "skipflag-only"  => \$opt_skiponly,
    "nopp"           => \$opt_no_postprocessing,
    "gff=s"          => \@external_gff_files
) or &usage;

&usage if $opt_help;
&usage if !@ARGV && -t STDIN && !@external_gff_files;
if ( $opt_pfsearch ) {
    die "OPTION CONFLICT: can't use option".
    " -reverse together with -w option (pfsearch)"
        if $opt_reverse;
    die "OPTION CONFLICT: can't use option".
    " -shuffle together with -w option (pfsearch)"
        if $opt_shuffle;
    die "OPTION CONFLICT: can't use option".
    " -e together with -w option (pfsearch)"
        if @entries;
    &usage if !@prosite_files &&
        !grep {/^$MOTIF_AC_REGEXP$/} @motifAC_or_userpattern;
    # integer cutoff forces raw scores
    $opt_raw = $1 if defined($opt_cutoff) and $opt_cutoff=~ /^(\d+)$/mg;
}
my $use_pfsearchV3 = index( $opt_pfsearch, 'pfsearchV3') != -1 ? 1 : 0;


my $scan_behavior;
$scan_behavior |= 1 if $opt_nooverlaps;
$scan_behavior |= 2 if $opt_includes;

$opt_format = "scan" unless defined $opt_format;
$opt_format =~ tr/A-Z/a-z/;
die "ERROR:Output format must be one of $formats_string\n"
    unless grep {$_ eq $opt_format} @formats;
my $opt_psa_or_msa = $opt_format eq "msa" || $opt_format eq "psa";

# user patterns specified with -p option
my @userpat;
# hash of specified (-p option) prosite ac (without specified user patterns)
my $specifiedPrositeMotifByAc={};
map {
    ( /^$MOTIF_AC_REGEXP$/ ?
        $specifiedPrositeMotifByAc->{$_} = 1 : push @userpat,$_ )
} @motifAC_or_userpattern;
# if no prosite ac specified: undef struct...
keys(%$specifiedPrositeMotifByAc) or $specifiedPrositeMotifByAc = undef;

# find default prosite.dat file
if ( !@prosite_files && !@userpat ) {
    if (defined $ENV{PROSITE} and -e "$ENV{PROSITE}/prosite.dat") {
        @prosite_files = "$ENV{PROSITE}/prosite.dat";
    }
    elsif ( defined $ENV{PROSITE} and -e $ENV{PROSITE} ) {
        @prosite_files = $ENV{PROSITE};
    }
    elsif ( defined $ENV{SPROT} and -e "$ENV{SPROT}/prosite/prosite.dat" ) {
        @prosite_files = "$ENV{SPROT}/prosite/prosite.dat";
    }
    elsif ( -e "prosite.dat" ) {
        @prosite_files = "prosite.dat";
    }
    else {
        die "prosite.dat file not found, please use the -d option";
    }
}

# -b option with no pathname specified: find default evaluator.dat
if ( defined($opt_miniprofiles) && !$opt_miniprofiles ) {
    if ( defined $ENV{PROSITE} and -e "$ENV{PROSITE}/evaluator.dat" ) {
        $opt_miniprofiles = "$ENV{PROSITE}/evaluator.dat";
    }
    elsif ( defined $ENV{SPROT} and -e "$ENV{SPROT}/prosite/evaluator.dat" ) {
        $opt_miniprofiles = "$ENV{SPROT}/prosite/evaluator.dat";
    }
    elsif ( -e "evaluator.dat" ) {
        $opt_miniprofiles = "evaluator.dat";
    }
    else {
        die "evaluator.dat file not found,".
        " please specify which evaluator.dat to use after the -b option";
    }
}

my %SkipFlag;
my %KnownFalsePos;
my @MotifInfo;

my $hideMotifByPSAC          = {};
my $files_miniprofiles       = {};
my $postProcessingByPSAC     = {};# PP data by PSAC (then by PP type)
my $motifRank4PostProcessing = {};

# dispatch table for post processing of (matching) target
# dataset against effector dataset (within hits on one sequence)
my $postProcessDispatchTable = {
    # match data =
    #  0         1     2   3     4(if profile, sequence if 1st match on a pattern) 5     6         7       8      9         10
    # [matchseq, from, to, pfac, pffrom (or full entry seq if pattern + option b), pfto, rawscore, nscore, level, leveltag, de

    # promote weak match (add 1 to level)
    # if effector has at least one hit at level>=$minlevel or 0
    # $target = ref to array that contains match data
    # for target motif on 1 sequence
    # $effector = ref to array that contains match data
    # for effector motif on 1 sequence
    'PROMOTED_BY' => sub {
        my ( $target, $effector, $minlevel ) = @_;
        return unless $target and $effector;

        $minlevel ||= 0;
        ( my $ac = $effector->[0]->[3] || '' ) =~ s/\|\w+$//;
        map { $_->[8]++, $_->[9] = "PROMOTED_BY_$ac" if $_->[8] < 0 } @$target
            if grep { defined( $_->[8] ) && $_->[8] >= $minlevel } @$effector;
    },

    # gives a "level" (in fact only shown in LevelTag; displayed level will
    # still be empty like for any pattern) to pattern matches:
    # if the effector detects the match at level 0, the match is given level 0,
    # LevelTag "(0)", else LevelTag "(-1)"
    # NOTE: effector hits are calculated here!
    'EVALUATED_BY' => sub {
        my ( $target, $effector, $mp_ac ) = @_;
        # note $effector is empty as matches will be calculated here!
        return unless $target and $mp_ac;

        my $sequence = $target->[0]->[4] or return;

        my $mini_tmpfile;
        my $mini_exist;
        my $miniprofile_dat = $opt_miniprofiles;
        if ( exists $files_miniprofiles->{ $mp_ac } ) {
            $mini_tmpfile = $files_miniprofiles->{ $mp_ac };
            $mini_exist   = "true";
        }
        else {# each new minifile is put in $files_miniprofiles
            $mini_tmpfile = tmpnam();
            open ( IN, $miniprofile_dat ) or
                die "can't open evaluator.dat file [$miniprofile_dat]";
            $/ = "\n//";# entry separator
                # note full sep. is \n//\n but in DOS format is \r\n//\r\n, so
                # only \n// matches both
            while ( <IN> ) {
                ( my $entry = $_ ) =~ s/^\s+//;
                if  ( $entry =~ /$mp_ac;/ ) {
                    $mini_exist = 1;
                    open (OUT,">$mini_tmpfile");
                    print OUT $entry."\n";
                    close OUT;
                    last;
                }
            }
            close IN;

            $files_miniprofiles->{ $mp_ac } = $mini_tmpfile if -e "$mini_tmpfile";
        }
        $/ = "\n";
        # if no miniprofile correspond to pattern then no evaluation
        return unless $mini_exist;

        my $struct = do_profile_scan( $mini_tmpfile, undef, $sequence );
        # returned: hash ref where key = mini profile AC or AC|ID!
        # take first (and only) value
        $effector = ( values %$struct )[0];
        foreach my $target_hit ( @$target ) {
            my $level = -1;
            foreach my $effector_hit ( @$effector ) {
                last if $effector_hit->[1] > $target_hit->[2];
                $level = 0, last if $effector_hit->[1] <= $target_hit->[1]
                                    and $effector_hit->[2] >= $target_hit->[2];
            }
            # $target_hit->[8] = $level; do not change level, just tag!
            $target_hit->[9] = "($level)";
        }
    },

    # downgrade good match (substract 1 to level)
    # if effector has at least one hit at level>=$minlevel or 0
    'DEMOTED_BY' => sub {
        my ( $target, $effector, $minlevel ) = @_;
        return unless $target and $effector;

        $minlevel ||= 0;
        map { $_->[8]--, $_->[9] = 'DEMOTED' if $_->[8] >= 0 } @$target
            if grep { defined( $_->[8] ) && $_->[8] >= $minlevel } @$effector;
    },

    # downgrade target/effector overlapping matches if their score is
    # smaller than competitor score
    'COMPETES_HIT_WITH' => sub {
        my ( $hit_set_a,$hit_set_b,$overlap ) = @_;
        return unless $hit_set_a and $hit_set_b;

        my @set_a = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] }
                    #grep {$_->[1]||=0;$_->[2]||=0;$_->[8]>=0}
                    @$hit_set_a;
                            # create array of () hits sorted by
                            # position in set a
        my @set_b = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] }
                    #grep {$_->[1]||=0;$_->[2]||=0;$_->[8]>=0}
                        @$hit_set_b;
                            # create array of () hits sorted by
                            # position in set b
        return unless @set_a and @set_b;
        ( my $ac_a = $hit_set_a->[0]->[3] || '' ) =~ s/\|\w+$//;
        ( my $ac_b = $hit_set_b->[0]->[3] || '' ) =~ s/\|\w+$//;
        $overlap=0 if !$overlap || $overlap!~m/^\d+$/;
            # allowed overlap size (only consider overlap matches if
            # overlap size > this value)
        foreach my $hit_a ( @set_a ) {
        # for each (sorted) level>=0 hits in set a
            my $a_start = $hit_a->[1] + $overlap;
            my $a_stop  = $hit_a->[2] - $overlap;
            next if $a_start > $a_stop;
            foreach my $hit_b ( @set_b ) {
            # for each (sorted) level>=0 hits in set b
                next if $a_start > $hit_b->[2];
                    # hit b is before examined hit a: next
                last if $a_stop  < $hit_b->[1];
                    # hit b is after examined hit a: last

                # overlap: demote hit with lowest score
                if ( $hit_b->[7] < $hit_a->[7] )   {
                    $hit_b->[8]--;
                    $hit_b->[9] = "OUTCOMPETED_HIT_BY_$ac_a";
                }
                elsif ( $hit_a->[7] < $hit_b->[7]) {
                    $hit_a->[8]--;
                    $hit_a->[9] = "OUTCOMPETED_HIT_BY_$ac_b";
                }
            }
        }
    },

    # downgrade target/effector matches from set with the lowest best score
    'COMPETES_SEQ_WITH' => sub {
        my ( $hit_set_a,$hit_set_b ) = @_;
        return unless $hit_set_a and $hit_set_b;
        ( my $ac_a = $hit_set_a->[0]->[3] || '' ) =~ s/\|\w+$//;
        ( my $ac_b = $hit_set_b->[0]->[3] || '' ) =~ s/\|\w+$//;

        my $max_a = 0.0; map { $max_a=$_->[7] if $_->[7]>$max_a } @$hit_set_a;
        my $max_b = 0.0; map { $max_b=$_->[7] if $_->[7]>$max_b } @$hit_set_b;
        if ( $max_a>$max_b ) {
            map {
                $_->[8]--;
                $_->[9] = "OUTCOMPETED_SEQ_BY_$ac_a"
            } @$hit_set_b
        }
        else {
            map {
                $_->[8]--;
                $_->[9] = "OUTCOMPETED_SEQ_BY_$ac_b"
            } @$hit_set_a
        }
    },

    # promote weak match (add 1 to level)
    # if effector has at least one hit with a (n)score >=$minscore
    # and there are more than one target match
    # Note: the target is usually the effector. Method used for repeats
    'PROMOTE_REPEAT' => sub {
        my ( $target, $effector, $cutoff ) = @_;
        return unless $target and $effector;

        $cutoff ||= 0;
        #(my $ac = $effector->[0]->[3] || '') =~ s/\|\w+$//;
        map { $_->[8]++, $_->[9] = "PROMOTED_REPEAT" if $_->[8] < 0 } @$target
            if ( @$target >1 and
                    grep {
                        defined( $_->[7] ) && $_->[7] >= $cutoff
                    } @$effector
                );
    },

    # promote weak match (add 1 to level), if there is at least one hit at
    # level>=0 (RDM1)
    'PROMOTER_RDM1' => sub {
        my $target = shift or return;
        # ('R?' in profile MA/CUT_OFF line)
        map { $_->[8]++, $_->[9] = "R" if $_->[8] < 0 } @$target
            if grep {defined($_->[8]) && $_->[8] >= 0} @$target;
    },

    # promote weak match (add 1 to level), if the sum of their (n)score
    # is >= threshold  OR one hit at level>=0 ('RR' in profile MA/CUT_OFF line)
    'PROMOTER_RDM1_OR_RDM2' => sub {
        my ( $target, $effector, $threshold ) = @_;
        return unless $target and $threshold;

        # RDM1 (>=1 level 0 match -> promote level<0)
        if ( grep { defined( $_->[8] ) && $_->[8] >= 0 } @$target ) {
            map { $_->[8]++, $_->[9] = "R" if $_->[8] < 0 } @$target;
        }
        # RDM2 (when no hit with level>=0) promote level<0 if hit score
        # sum >= treshhold
        else {
            my $sum = 0.0; map { $sum+=$_->[7] || 0; } @$target;
            map { $_->[8]++, $_->[9] = "r" if $_->[8] < 0 } @$target
                if $sum >= $threshold;
        }
    },

    # same as PROMOTE_RDM1_OR_RDM2 (but on repeats from '!' profile,
    # only with -a option)
    'PROMOTER_RDM1_OR_RDM2_4!REP' => sub {
        my ( $target, $effector, $threshold ) = @_;
        return unless $target and $threshold;

        # RDM1 (>=1 level 0 match -> promote level<0)
        if ( grep { defined($_->[8]) && $_->[8] >= 0 } @$target ) {
            map { $_->[8]++, $_->[9] = "?R" if $_->[8] < 0 } @$target;
        }
        # RDM2 (when no hit with level >= 0) promote level < 0 if hit score
        # sum >= treshhold
        else {
            my $sum = 0.0; map { $sum += $_->[7] || 0; } @$target;
            map { $_->[8]++, $_->[9] = "?r" if $_->[8] < 0 } @$target
                                                    if $sum >= $threshold;
        }
    },

    'FTREP' => sub {
        my ( $target, $effector, $paramstr ) = @_;
        return unless $target;

        my ( $max_score, $ac_cutoff, $size, $gf, $filter ) =
            split( '\|', $paramstr );

        my $motif_param = {
            'motif_id'   => 'fake',
            'motif_acc'  => 'fake',
            'motif_len'  => $size,
            'maxscore'   => $max_score,
            'filter'     => undef,
            'normscore1' => undef,
        };

        my $matchlist = [ map {
            my $estart = $_->[1] - ($_->[4]||1)  +1;
            my $estop  = $_->[2] - ($_->[5]||-1) -1;
            my $size  = 1 + $estop - $estart;
            {
                'match_start'  => $_->[1],
                'match_end'    => $_->[2],
                'match_nscore' => $_->[7],
                'match_rscore' => $_->[6],
                'match_pos'    => $size/2.0 + $estart,
                'match_size'   => $size
            }
        } @$target ];

        # Group  matches
        my $groups = group_matches( $matchlist, $gf * $size );
        return unless @$groups;

        # Fourier
        my ( $ft, $fa ) = FT( $matchlist, $groups, $motif_param, 1 );

        # Autocorrelation
        my $ac = AC( $matchlist, $groups, $fa );

        # Get evalues and select matches
        select_match( $matchlist, $groups, $ac, $ac_cutoff, $filter );

        for ( my $i = 0 ; $i < @$target ; $i++ ) {
            my $match      = $target->[ $i ];
            my $ftrepmatch = $matchlist->[ $i ];
            if ( $match->[ 8 ] >= 0 ) {
                # good match: leave it as-is
            }
            elsif ( $ftrepmatch->{ selected } ) {# promote (to level 0)
                $match->[ 8 ] = 0;
                $match->[ 9 ] = 'PROMOTED_FTREP'
            }
            else {# demote (-1 to -2)
                $match->[ 8 ] = -2;
                $match->[ 9 ] = 'DEMOTED_FTREP'
            }
        }
    }

};

my $allowBidirectionalPP = {
    'PROMOTED_BY' => 1,
    # note: COMPETES_... are 'bidirectional' inside their
    # $postProcessDispatchTable sub but we don't wan't them to be called
    # bidirectionally (in pp_scan)...
};


main();
exit $errcode;


##############################################################################
# Methods

sub tmpnam {
    my $tmp;
    do {
        $tmp = $TMPDIR.$SLASH."ps${TMP_BASE}-".$TMP_COUNTER++.".tmp";
    } while (-e $tmp);
    return $tmp;
}

# -------------------------- output methods --------------------------

# format a field with a certain width
sub pf {
    return $_[0] . ( " " x ( $_[1] - length $_[0] ) )
}

my $HIT_COUNT = 0;
sub dispHits {# display hits (for 1 sequence, 1 motif)
    my ( $header, $sq, $hits, $seqid, $de, $aclist, $psac, $psid, $psde ) = @_;
    return if !$hits || !$psac || $hideMotifByPSAC->{ $psac };

    $sq ||= ''; $seqid ||= ''; $de ||= ''; $psid ||= ''; $psde ||= '';
    ( $de = $hits->[0]->[10] || '' ) =~ s/\.\s*$// if !$de && $opt_pfsearch;
    $de =~ s/[\n\r\t]/ /g;

    # Hit array:[matchseq, from, to, ac, pffrom, pfto, rawscore, nscore, level, leveltag, de,  [repstuff!?]]]
    #            0         1     2   3   4       5     6         7       8      9         10    11
    #   Note: might ~vary depending if pattern/profile pfscan, pfsearch!

    # remove unwanted matches with level < opt_level
    my $visible_hits;
    @$visible_hits =
        grep {
            # if nscore [7] is undef = match is from a pattern : show match
            # (even if has level < 0 as now pattern might (~canceled) have
            # level if used with -b miniprofile option!) or if match
            # level [8] >= user requested minimal level
            !defined( $_->[7] ) || $_->[8] >= $opt_level
        } @$hits;
    @$visible_hits or return;

    unroll_hits( $visible_hits );
    if ( defined $opt_maxhits ) {
        splice( @$visible_hits, $opt_maxhits - $HIT_COUNT );
        exit 0 if $HIT_COUNT >= $opt_maxhits ;
    }
    my $hit_count  = @$visible_hits;
    $HIT_COUNT    += $hit_count;
    return unless $hit_count;
    return if defined ( $opt_minhits ) and $opt_minhits > $hit_count;

    $visible_hits->[0]->[4] = undef unless defined( $visible_hits->[0]->[7] );
        # clear sequence put in [4] in 1st pattern hit!

    print $header if defined $header;

    if ( $opt_format eq "ipro" ) {
    # fake & hidden output: gff with additional post-processing for patterns!
    # (By pattern) if one pattern hit (on this entry) has been promoted by
    # a corresponding mini profile (EVALUATED_BY PP) to "(0)", promote all
    # other matches (of the same pattern, in the same entry)
    # (n.b. only patterns have '(<level>)' leveltags)
    # FIXME: why not make this a standard (& public) PP behavior
    # (e.g. in EVALUATED_BY PP)
        my $ok_by_ac = {};
        map{ $ok_by_ac->{ $_->[3] } = 1 if $_->[9] eq '(0)' } @$visible_hits;
        @$visible_hits =
            map { $_->[9] = '(0)' if  $ok_by_ac->{ $_->[3] };
                                                    $_ } @$visible_hits;
                # n.b. pattern leveltags are only (-1) or (0) (if not, this
                #  will give unexpected results!)
        $opt_format    = 'gff';
    }

    if ( $opt_format eq "fasta" or $opt_psa_or_msa ) {
        for my $hit ( @$visible_hits ) {
            my ( $subseq, $from, $to, $_psac, $pffrom, $pfto, $rawscore,
                 $nscore, $leveln, $levelt, $seqde ) = @$hit;
            my $print_level = defined( $levelt ) ? "L=$levelt " :
                defined($leveln) ? "L=$leveln " : "";
            if ( $opt_pfsearch ) {
                $opt_psa_or_msa = 1;
                print ">$seqid/$from-$to : $de : $psac $print_level\n";
            }
            else {
                print ">$seqid/$from-$to : $psac $psid $print_level\n";
            }
            if ( $subseq and $opt_psa_or_msa ) {# pfscan output
                $subseq =~ s/\n?$/\n/; # add \n to scanPattern output
                print $subseq;
            }
            else {
                $subseq = substr( $sq, $from-1, $to-$from+1 );
                while ( $subseq =~ /(.{1,60})/g ) {
                    print "$1\n";
                }
                print "\n" if $subseq eq "";
            }
        }
    }
    elsif ( $opt_format eq "pff" || $opt_format eq "epff" ) {
        for my $hit ( @$visible_hits ) {
            my @pff = @$hit;
            $pff[8] = $pff[9] if defined $pff[9];
                # move LevelTag into Level if defined
            pop @pff while @pff>9; # remove fields beyond numeric level
            my $subseq = shift @pff;  # remove subseq
            $pff[2] = $psac || $pff[2];
            if ($opt_format eq "epff") {
                push @pff, "" while @pff < 8;
                $subseq =~ s/\s//g;
                push @pff, $subseq;
            }
            print $seqid, "\t", join( "\t", @pff ), "\n";
        }
    }
    elsif ( $opt_format eq "gff" ) {
        for my $hit ( @$visible_hits ) {
            my ( $subseq, $from, $to, $_psac, $pffrom, $pfto, $rawscore,
                 $nscore, $leveln, $levelt, $seqde ) = @$hit;
            print join( "\t", $seqid, "ps_scan|v$VERSION", $psac,
                                $from, $to, $nscore || ".", ".", "." );
            my @attr;
            if ( $psid ) {
                $psid =~ s/.*\|//; # FIXME!: why!?
                push @attr, "Name \"$psid\"";
            }
            push @attr, "AccessionNumbers " .
                join( " ", map {"\"$_\""} @$aclist )
                if defined $aclist and @$aclist;
            push @attr, "Level $leveln" if defined $leveln;
            push @attr, "LevelTag \"$levelt\"" if defined $levelt;
            push @attr, "RawScore $rawscore" if defined $rawscore;
            push @attr, "FeatureFrom $pffrom"
                if defined $pffrom and $pffrom =~ /\d+/;
            push @attr, "FeatureTo $pfto" if defined $pfto;
            $subseq =~ s/\s//g;
            push @attr, "Sequence \"$subseq\"" if defined $subseq;
            push @attr, "SequenceDescription \"$de\"" if $de;
            push @attr, "SkipFlag 1" if $SkipFlag{$psac};
            push @attr, "KnownFalsePos $KnownFalsePos{$psac}"
            if exists $KnownFalsePos{ $psac };
            print "\t", join " ; ", @attr if @attr;
            print "\n";
        }
    }
    elsif ( $opt_format eq "sequence" ) {
        print ">$seqid : $de : $psac\n", map { "$_\n" }
            $sq =~ /(.{1,60})/g if @$visible_hits;
    }
    else {# default output format
        # fasta like header
        if ( $opt_pfsearch ) {
            print ">$seqid : $de : $psac\n";
        }
        else {
            print ">$seqid : $psac $psid $psde\n";
        }
        # hits
        for my $hit ( @$visible_hits ) {
            my ( $subseq, $from, $to, $_psac, $pffrom, $pfto, $rawscore,
                 $nscore, $leveln, $levelt, $seqde ) = @$hit;
            my $print_level = defined( $levelt ) ? " L=$levelt" :
                defined( $leveln ) ? " L=$leveln" : "";
            my $fromto = "$from - $to";
            print " " x ( 13-length $fromto ), $fromto;
            if ( $subseq ) {# pfscan output
                $subseq =~ s/\n?$/\n/;# add \n to scanPattern output
                $subseq =~ s/^(?<!\A)(.*)/ $1/mg;
                $subseq =~
                    s/^(.*)/$1 . (" " x (60-length $1)) . $print_level/e
                        if $print_level;
                print "  ", $subseq;
            }
            else {
                $subseq = substr( $sq, $from - 1, $to-$from + 1 );
                my $notfirst;
                while ( $subseq =~ /(.{1,60})/g ) {
                    print " " x 13 if $notfirst++;
                    print " $1\n";
                }
                print "\n" if $subseq eq "";
            }
        }
    }
}

sub showMatchList {# show (protein sequence) match list of (selected) motifs
    for my $pattern ( @MotifInfo ) {
        my ( $psac, $psid, $type, $psde, $pat, $skip ) = @$pattern;
        next if $type eq "MATRIX";# not available for profiles
        my $n_hits = 0;
        my $n_match = 0;
        print "$psid\n";
        for my $sprotfile ( @ARGV ) {
            open SPROTFILE, $sprotfile or die "Cannot open $sprotfile : $!";
            my $entry = "";
            # Note: only works on Swiss-Prot flat file format not on fasta ...
            while ( <SPROTFILE> ) {
                $entry .= $_;
                if ( /^\/\// ) {
                    if ( $entry =~
                            /^\s*ID\s+(\w+)/m ) {
                        my $id = $1;
                        if ( $entry =~ /^\s*AC\s+(\w+)/m ) {
                            my $ac = $1;
                            my @de = /^\s*DE\s+(.+)/mg;
                            my $de;
                            my $add_space = 0;
                            for ( @de ) {
                                $de .= " " if $add_space;
                                $de .= $_;
                                $add_space = !/-$/;
                            }
                            if ( $entry =~
                                    /^\s*SQ\s+SEQUENCE\b.*\n((.+\n)+)/m ) {
                                my $sq = $1;
                                $sq =~ tr/A-Z//cd;
                                $de = "$id $de" if $id and $de;
                                my $hits;
                                $hits = scanPatternWrapper([$pat, $sq,
                                    $scan_behavior, $opt_max_x,
                                    $opt_miniprofiles], $psid);
                                if ($hits and @$hits) {
                                    print pf($ac, 6), ", ", pf($id, 10),
                                        ", T;       ", scalar(@$hits), "\n";
                                    $n_hits += @$hits;
                                    $n_match++;
                                }
                            }
                            else {
                                warn "No sequence found in entry $id";
                            }
                        }
                    }
                    $entry = '';
                }
            }
        }
        print " ", pf( $n_match, 13 ), " ", $n_hits, "\n";
    }
}


sub process_external_gff {
# process existing external gff ps_scan results, post-process (unless --nopp)
# and/or redisplay them in any available output format
# (usefull when scanning profile by profile on a cluster = no correct pp
#  could be performed, then reuse ps_scan with --gff on collected results
#  where correct pp (as all the distinct motif matches are combined) could be
#  now performed. Caution: memory usage!!...
# Nov25 2010

    my $psinfo_bypsac = {};
    map{# store motif info by motif ac
        $_->[0] and $psinfo_bypsac->{ $_->[0] } = $_;
    } @MotifInfo;

    # parse input gff
    my $hits_by_psac_by_seq = {};
    my $info_by_seqid       = {};
    foreach my $gff_file ( @external_gff_files ) {
        next unless -e $gff_file;
        open( FH, $gff_file ) or die "can't open [$gff_file] file";
        while( <FH> ) {# read gff line by line
            s/\r?\n$//;
            # parse core elements + extra (attributes)
            my ( $id, $nil1, $psac, $from, $to, $nscore,
                 $nil2, $nil3, $extra                       ) = split /\t/;
            $extra ||= '';

            # parse gff attributes
            #my $psid = $extra =~ /Name \"(\w+)\"/ ? $1 : undef;
            my $aclist = [];
            if ( $extra =~ /AccessionNumbers (.+?) ;/g ) {
                @$aclist = map { s/['"]//g; $_ } split / /, $1
            }
            my $leveln   = $extra =~ /Level ([-\d]+)/ ? $1 : undef;
            my $levelt   = $extra =~ /LevelTag \"(.+?)\"/ ? $1 : undef;
            my $rawscore = $extra =~ /RawScore (\d+)/ ? $1 : undef;
            my $pffrom   = $extra =~ /FeatureFrom ([-\d]+)/ ? $1 : undef;
            my $pfto     = $extra =~ /FeatureTo ([-\d]+)/ ? $1 : undef;
            my $subseq   = $extra =~ /Sequence \"(.+?)\"/ ? $1 : undef;
            my $seqde    = $extra =~ /SequenceDescription \"(.+?)\"/ ? $1 : '';

            push @{ $hits_by_psac_by_seq->{ $id }->{ $psac } },
                [ $subseq, $from, $to, $psac, $pffrom, $pfto, $rawscore,
                  $nscore, $leveln, $levelt, '' ];
            $info_by_seqid->{ $id }->{ de }  = $seqde;
            $info_by_seqid->{ $id }->{ acs } = $aclist;
        }
    }
    # for each seq do pp
    foreach my $seqid ( keys %$hits_by_psac_by_seq ) {
        postProcess( $hits_by_psac_by_seq->{ $seqid } )
                                unless $opt_no_postprocessing;
                # if no pp: usefull only for format conversion
        foreach my $psac ( sort { $a cmp $b }
                            keys %{ $hits_by_psac_by_seq->{ $seqid } } ) {
        # display hits (foreach motif ac)
            my $hits   = $hits_by_psac_by_seq->{ $seqid }->{ $psac };
            my $seqde  = $info_by_seqid->{ $seqid }->{ de } || '';
            my $aclist = $info_by_seqid->{ $seqid }->{ acs } || [];
            my $psid   = $psinfo_bypsac->{ $psac }->[ 1 ] || $psac;
            my $psde   = $psinfo_bypsac->{ $psac }->[ 3 ] || '';
                # n.b. psde not in gff source, but in motif file
            dispHits( undef, '', $hits, $seqid, $seqde, $aclist, $psac,
                                                                $psid, $psde );
        }
        $hits_by_psac_by_seq->{ $seqid } = undef;
        $info_by_seqid->{ $seqid }       = undef;
    }
}


# -------------------------- scan methods --------------------------

# scan entries in specified sequence file 'collection' against motif
# collection and display results through called subs)
sub scanSeqFile {
    my $seqfile = shift or return;

    my $entry = "";
    my $opt_fasta;
    my @motifs_4_normal_scan;
    my @motifs_4_pfsearch;
    my $psinfo_bypsac = {};
    map {
        # store motif info by motif ac
        $_->[0] and $psinfo_bypsac->{$_->[0]} = $_;
        $opt_pfsearch && $_->[2] && $_->[2] eq 'MATRIX' ?
            push @motifs_4_pfsearch, $_ : push @motifs_4_normal_scan, $_
    } @MotifInfo;
    my $all_hits_bypsac_byseqid = {};
    # -W PFSEARCH SCAN (every profile against sequence collection)
    if (grep {$_->[7]} @motifs_4_pfsearch) {
        # with pfsearch (faster to scan each profile
        # against a sequence collection):
        # use whole sequence collection files, do not extract each entry!
        foreach my $prf_info (@motifs_4_pfsearch) {# for each Motifinfo
            next if $prf_info->[5]; # next if skip...
            my $psac = $prf_info->[0] or next;
            # needs an associated tmp profile file
            my $profile_tmp_file = $prf_info->[7] or next;
            my $psid = $prf_info->[1] || '';
            # scan seq collection file against 1 profile with pfsearch
            #! with pfsearch do_profile_scan
            # will return hits in a hash (ref) by seq ac|id
            # (not by psac like with pfscan)!
            my $pfsearchhits =
                do_profile_scan($profile_tmp_file, $seqfile, undef, $prf_info);
            # group all diff hits by psac and by seq :
            foreach my $seq_id (keys %$pfsearchhits) {
                # for each hit
                foreach my $hit (@{$pfsearchhits->{$seq_id}}) {
                    # set Name field (in -o gff) to psid (instead of seqid)
                    # = same as with normal scan
                    # (my $seqid=$seq_id)=~s/^\w+\|//;
                    # FIXME: should be done at the do_profile_scan level
                    # fix seqid (do_profile_scan sometimes
                    # returns ac|id instead of id)
                    $hit->[3] = $psid;
                    push @{$all_hits_bypsac_byseqid->{$seq_id}->{$psac}},
                        $hit;
                }
            }
        }
        # post process and show results (all profile hits grouped by protein)
        # note: $all_hits_bypsac_byseqid only used for -w pfsearch profile.
        # for normal scan post-processed and displayed is done in scanSeq!
        # side effect: with -w pfsearch,
        # post processing between profiles and patterns
        # won't work

        # for each matched sequence ac
        foreach my $seq_id ( keys %$all_hits_bypsac_byseqid ) {
            my $pfhits_bypsac = $all_hits_bypsac_byseqid->{ $seq_id } or next;
            # do postprocessing (on all sequence hits)
            postProcess( $pfhits_bypsac ) unless $opt_no_postprocessing;
            foreach my $psac ( keys %$pfhits_bypsac ) {# display hits
                my $psid = $psinfo_bypsac->{$psac}->[1] || $psac;
                # won't be displayed (in pfsearch mode)!
                # but in case this will be changed
                my $psde = $psinfo_bypsac->{$psac}->[3] || '';
                dispHits( undef, '',
                    $all_hits_bypsac_byseqid->{$seq_id}->{$psac},
                    $seq_id, '', undef, $psac, $psid, $psde );
            }
        }
    }
    # 'NORMAL' SCAN (every sequence against profile collection)
    if ( @motifs_4_normal_scan ) {# scan entries in sequence file collection
    	my $seqfile_h = new IO::File $seqfile
            or die "Cannot open $seqfile: $!";
        while ( <$seqfile_h> ) {# read sequence file
            my $line = $_;
            if ( /^>(.*)/ ) {# entry is in fasta
                #(my $seqid = $1) =~ s/\s+\.*$//;
                # $all_hits_bypsac_byseqid->{$seqid}
                # = scanFromFastaEntry($entry);
                scanFromFastaEntry( $entry );
                # note: do not store hits (by ac) into
                # $all_hits_bypsac_byseqid
                # might use too much memory!; side effect: with -w pfsearch,
                # post processing between profiles and patterns won't work
                $opt_fasta = 1;
                $entry     = '';
            }
            $entry .= $line;
            # entry separator (in Swiss-Prot format) => an entry has been
            # fully read
            if ( /^\/\// ) {
                $opt_fasta = 0;
                my @id = $entry =~ /^\s*ID\s+(\w+)/mg;
                my $ac_lines;
                $ac_lines .= $1 while $entry =~ /^\s*AC\s+(.*)/mg;
                my @ac;
                while ( $ac_lines =~ /(\w+)/g ) { push @ac, $1; }
                if ( @id ) {
                    if ( not (@entries) or
                        grep { my $ent = $_;
                                grep{ $_ eq $ent } @id, @ac } @entries) {
                        my $id = $id[ 0 ];
                        # get entry DE:
                        my $de = #? 'new' format
                            $entry =~ /^DE   (?:Rec|Sub)Name: Full=(.+);/m ?
                                $1 : '';
                        unless ( $de ) {# old DE format
                            my @de = $entry =~ /^\s*DE\s+(.+)/mg;
                            my $add_space = 0;
                            for ( @de ) {
                                $de .= " " if $add_space;
                                $de .= $_;
                                $add_space = !/-$/;
                            }
                        }

                        if ($entry =~ /^\s*SQ\s+SEQUENCE\b.*\n((.+\n)+)/m) {
                            my $sq = $1;
                            $sq =~ tr/A-Z//cd;
                            # $all_hits_bypsac_byseqid->{$id} =
                            # scanSeq($id, \@ac, $de || $id, $sq);
                            # note: do not store hits (by ac) into
                            # $all_hits_bypsac_byseqid would use too much
                            # memory; side effect: with -w pfsearch,
                            # post processing between a profile and and
                            # pattern won't work
                            scanSeq($id, \@ac, $de || $id, $sq);
                        }
                        else {
                            warn "No sequence found in entry $id";
                        }
                    }
                }
                # ignore entries which have "id" in lowercase
                elsif ($entry =~ /^\s*id /m) {
                }
                elsif ($entry =~ /(.*\S.*)/) {
                    warn "Bad sequence found in file, first line: $1\n";
                    $errcode = 1;
                }
                $entry = "";
            }
        }
        close $seqfile_h;
        if ($entry =~ /^>/) {# process last fasta entry
            scanFromFastaEntry($entry);
        }
        elsif ($entry =~ /(.*\S.*)/) {
            warn "Bad sequence found in file, first line : $1\n";
            $errcode = 1
        }
    }
}

# scan from 1 fasta entry (from a sequence entry collection)
sub scanFromFastaEntry {
    my $entry = shift or return;
    return unless $entry =~ s/^>((\S*).*)\n//;
    my ($fasta_header, $primary_id) = ($1, $2);
    return if defined($opt_filterheader)
        and $fasta_header !~ /$opt_filterheader/o;
    if (not (@entries) or grep {$_ eq $primary_id} @entries) {
        $entry =~ tr/A-Z//cd;
        return scanSeq($primary_id, [], $fasta_header, $entry);
    }
}

# scans one sequence against prosite motifs (only normal scan; if
# opt_pfsearch, all scans are done before in scanSeqFile)
sub scanSeq {
    my ($id, $aclist, $de, $sq) = @_;
    if ($opt_reverse) {# reverses sequence (if opt on)
        $sq = reverse $sq;
    }
    if ($opt_shuffle) {# shuffles sequence (if opt on)
        srand 0;
        my @seq = grep {$_ ne "\n"} split(//,$sq);
        $sq = "";
        for (my $start_win = 0; $start_win < @seq;
                $start_win += $opt_shuffle) {
            my $stop_win = $start_win + $opt_shuffle - 1;
            $stop_win = @seq - 1 if $stop_win >= @seq;
            my @residues = @seq[$start_win..$stop_win];
            while (@residues) {
                $sq .= splice(@residues, int(rand(scalar @residues)) ,1);
            }
        }
    }
    # scan sequence against profile first
    my $all_pfhits = {};
    # if there are some profile to be scanned...
    if ($scan_profiles && !$opt_pfsearch) {
        # in normal scan (not -w pfsearch) only 1 temp file is used (or
        # user specified prosite files)
        my $file_source = ($last_profile_tmp_filename ?
            [$last_profile_tmp_filename] : \@prosite_files);
        foreach my $profile_file (@$file_source) {#
            my $pfhits = do_profile_scan($profile_file, undef, $sq);
                # store results (do_profile_scan result is a hash ref
                # to [[hitdata],[],[]...] by profile AC)
            foreach my $ps_acid (keys %$pfhits) {
                (my $psac = $ps_acid) =~ s/\|.+//;
                $all_pfhits->{$psac} = $pfhits->{$ps_acid};
            }
        }
    }
    # loop through each motif, scan sequence against MATRIX, builds hit
    # structure
    my $hits_by_motif_ac = {};
    foreach (@MotifInfo) {# for each motif info 'data'
    	my ($psac, $psid, $type, $psde, $pat, $skip) = @$_;
        next if $skip;
        my $hits;
        # MATRIX (profile) (profiles were run just before...)
        if ($type eq "MATRIX") {
            # if opt_pfsearch: scan and results display already done
            # (in scanSeqFile sub)
            next if $opt_pfsearch;
            $hits = $hits_by_motif_ac->{$psac} = $all_pfhits->{$psac} or next;
        }
        else {# PATTERN
            warn("Empty pattern for $psac\n"), next unless $pat;
            $hits = scanPatternWrapper([$pat, $sq, $scan_behavior,
                        $opt_max_x, $opt_miniprofiles], $psid);
            $hits_by_motif_ac->{$psac} = $hits if @$hits;
        }
        dispHits(undef, $sq, $hits, $id, $de, $aclist, $psac, $psid, $psde)
            if $opt_no_postprocessing;
            # if no pp matches can be displayed right now (as they go)
            # otherwise, first have to collect all matches (this loop), then pp
    }
    unless ($opt_no_postprocessing) {
    # if pp: only can be done after having collected all matches
        # do post-processing (between all hits within one entry)
        # note: messy: pp.+disp. performed here, but also before in scanSeqFile
        # for -w pfsearch profiles!
        postProcess($hits_by_motif_ac);
        # display results (after post processing) for each motif
        # (keep order of motifs in prosite motif file or order of
        # specified (-p) motifs)
        foreach (@MotifInfo) {
            my ($psac, $psid, $type, $psde, $pat, $skip) = @$_;
            my $hits = $hits_by_motif_ac->{$psac} or next;
            dispHits(undef, $sq, $hits, $id, $de, $aclist, $psac, $psid, $psde);
        }
    }
    return $hits_by_motif_ac;
}


# run pfscan on 1 sequence against profile collection or pfsearch (if -w...):
# 1 profile against sequence collection sequence can be specified either
# as a filename or a string
sub do_profile_scan {
    my ($PROSITE, $seqfile_to_scan, $sequence, $prf_info) = @_;
    #return {} unless grep {$_->[2] eq "MATRIX"} @MotifInfo;
    my $PFSCAN_TMP = tmpnam();
    my $level_arg = defined($opt_level) ? "$opt_level" : "0";
    my $level_min = 0;
    my(@pre_command, @post_command);
    if ($opt_pfsearch) {# pfsearch (1 prf against a seq collection)
    # retrive cut off value from profile to set pfsearch -C option
    # (not needed with pfsearch 2.3 (can use -l <level>)...
    # but to be compatible with 2.2...)
        my $pfsearch_cutoff;
        # get cut-off corresponding to level -1
        # (needed for the detection of the repeats)
        map {
            $pfsearch_cutoff = ($opt_raw ? $_->{SCORE} : $_->{N_SCORE})
                if ($_->{LEVEL} eq '-1');
            $level_min = $_->{LEVEL} if ($_->{LEVEL} < $level_min);
        } @{$prf_info->[6]}
            if ($prf_info && $prf_info->[6]);
        # if option pfsearch is selected, get the user specified C=?? parameter
        # if not defined: pfsearch is run as default at L=-1
        # (pfsearch aut. detects if C is SCORE (integer) or N_SCORE(float))
        my $cutoff = defined($opt_cutoff) ? "$opt_cutoff" :
            "$pfsearch_cutoff";$cutoff ||= '0';
        # detect format
        my $fasta = $sequence ? "-f" : "";
        if ( open( DETECT, $seqfile_to_scan ) ) {
            while ( <DETECT> ) {
                next unless /\S/;
                $fasta = /^\s*>/ ? "-f" : "";
                last;
            }
            close DETECT;
        }
        @pre_command = ( $use_pfsearchV3 ? 
        	"$opt_pfsearch -o1 -c$cutoff $PROSITE" : "$opt_pfsearch $fasta -lxz $PROSITE" );
        	# p.s. if input is not fasta, pfsearchV3 will fail!
        @post_command = ( $use_pfsearchV3 ? "" : "C=$cutoff" );
    }
    else {# 'normal' pfscan (1seq against a prf collection)
        # if the user select a Level L higher than L=0 for match detection
        # the methods for the detection of repeats are not applied otherwise
        # pfscan is run as default at L<=-1
        if ($level_arg eq 0) {
            $level_arg = -1;
        }
        @pre_command = "$PFSCAN -flxz -v";
        @post_command = "$PROSITE L=$level_arg";
    }
    my $out;
    if ( $use_pfsearchV3 || $NO_DIRECT_PIPE || defined $seqfile_to_scan ) {
        my $seqfile;
        unless (defined $seqfile_to_scan) {
            $seqfile = tmpnam();
            open SEQ_TMP, ">$seqfile" or die "Cannot create $seqfile : $!";
            print SEQ_TMP ">seq for pfscan\n";
            print SEQ_TMP "$1\n" while $sequence =~ /(.{1,60})/g;
            close SEQ_TMP;
        }
        else {
            $seqfile = $seqfile_to_scan;
        }
        my $cmd = "@pre_command $seqfile @post_command > $PFSCAN_TMP";
        # launch pftool scan command
        system $cmd and die "Could not execute $cmd";
        unlink $seqfile unless defined $seqfile_to_scan;
        my $pfscan_fh = new IO::File($PFSCAN_TMP)
            or die "Cannot open $PFSCAN_TMP: $!";
        $out = scanProfiles($pfscan_fh, $level_min-1);# parse scan output
        close $pfscan_fh or die "Error $? with $PFSCAN_TMP";
    }
    else {
        # directly feed data to pfscan via pipe, do not use any temporary
        # files
        require IPC::Open2;
        my ($reader, $writer);
        my $cmd = "@pre_command - @post_command";
        my $pid =
        eval {
            IPC::Open2::open2($reader, $writer, $cmd)
            or die "Could not fork pipe to $cmd: $!";
        };
        if ($@) {
            die "$@\n" . '>' x 62 . "\nERROR: 'pfscan' execution failed.
            Check pfscan is in your PATH\n" . '>' x 62 . "\n";
        }
        local $/ = \32767;# buffer size
        print $writer ">seq for pfscan\n";
        print $writer "$1\n" while $sequence =~ /(.{1,60})/g;
        close $writer;
        $out = scanProfiles($reader, $level_min-1);
        close $reader or die "Error $? with $cmd";
        waitpid $pid, 0; #avoid defunct kid processes
    }
    unlink $PFSCAN_TMP;
    if ($opt_format eq "msa") {# if output format is MSA run psa2msa
        for my $ac (keys %$out) {
        # for each prosite profile ac write to tmp file (reuse $PFSCAN_TMP)
            open PFSCANTMP, ">$PFSCAN_TMP"
                or die "Cannot create $PFSCAN_TMP : $!";
            # for each hit (by ps ac)
            for (my $i = 0; $i < @{$out->{$ac}}; $i++) {
                my $hit = $out->{$ac}->[$i];
                print PFSCANTMP ">$i\n$hit->[0]";
            }
            close PFSCANTMP;
            my $PSA2MSA_TMP = tmpnam();
            my $cmd = "$PSA2MSA $PFSCAN_TMP > $PSA2MSA_TMP";# run psa2msa
            system $cmd and die "Cannot execute $cmd";
            open MSATMP, $PSA2MSA_TMP or die "Cannot read $PSA2MSA_TMP : $!";
            my %msa;
            my $cur_pos;
            while(defined (local $_ = <MSATMP>)) {# read psa2msa result output
                if (/^>(\d+)/) {
                    $cur_pos = $1;# hit number
                }
                elsif (defined $cur_pos) {
                    $msa{$cur_pos} .= $_;# sequence 'alignment'
                }
            }
            close MSATMP;
            unlink $PSA2MSA_TMP;
            while (my ($number, $seq) = each %msa) {
                # change match sequence (?) for hit (for prosite $ac)
                # number $number
                $out->{$ac}->[$number]->[0] = $seq;
            }
        }
    }
    return $out;
}

# replace hits on a repeat region by individual repeat elements
# (with pfsearch/pfscan 2.3)
sub unroll_hits {
    my ($hits) = @_;
    return unless(grep {$_->[11] && @{$_->[11]}} @$hits);
    for (my $i = 0; $i < @$hits; $i++) {
        my ($subseq, $from, $to, $pfid, $pffrom, $pfto, $rawscore,
            $nscore, $leveln, $levelt, $seqde, $subhits) = @{$hits->[$i]};
        next unless($subhits && @$subhits);
        map {$_->[8] = $leveln; $_->[9] = $levelt; $_->[10] = $seqde}
            @$subhits;
        splice @$hits, $i--, 1, @$subhits;
    }
}

sub scanPatternWrapper {
    my ($args, $id) = @_;
    my $out = scanPattern(@$args);
    # in PSA format, remove the '.' character from inserts.
    # these can be reintroduced with the 'psa2msa' program.
    if ($opt_format eq "psa") {
        $_->[0] =~ s/\.//g for @$out;
    }
    if ($id) {
        $_->[3] = $id for @$out;
    }
    return $out;
}

sub prositeToRegexpWrapper {
    my $out = prositeToRegexp(@_);
    unless (defined $out) {
        print STDERR "ps_scan.pl: Syntax error in pattern".
        " at position $Prosite::errpos\n";
        print STDERR "$_[0]\n";
        print STDERR " " x $Prosite::errpos, "^--- $Prosite::errstr\n";
        exit 1;
    }
    return $out;
}

# -------------------------- post processing --------------------------

sub pp_scan {
    my ($hits, $all_psac_in_seq, $intra_not_inter) = @_;
    return unless($hits && $all_psac_in_seq);
    my $seen = {};
    foreach my $target_psac (@$all_psac_in_seq) {
        # loops through all motifs in result structure
        my $hit_target_set = $hits->{$target_psac} or next;
        @$hit_target_set or next;
        foreach my $potential_pp (grep {($intra_not_inter ? $_->{effector} eq
            $target_psac : $_->{effector} ne $target_psac) }
            @{$postProcessingByPSAC->{$target_psac}}) {
            # loops through all motif associated pp
            # (order = order of PP lines...)
            # (+filter for intra/inter pp...)
            #print STDERR "##PP ($intra_not_inter) target [$target_psac]",
            #"-> pp effector[$potential_pp->{effector}] ",
            #"type[$potential_pp->{type}] value[$potential_pp->{value}]\n";
            my $effector_psac = $potential_pp->{effector} or next;
            my $pp_type = $potential_pp->{type} or next;
            next if $seen->{$pp_type}->{$target_psac}->{$effector_psac};
            #next if !$allowBidirectionalPP->{$pp_type} &&
            #            $seen->{$target_psac}->{$effector_psac};
                # next if 'reverse' (target<->effector) pp was already seen!
                # might happen with COMPETE_ PPs where PP is specified in
                # both effector & target...
            $seen->{$pp_type}->{$effector_psac}->{$target_psac} = 1
                unless $allowBidirectionalPP->{$pp_type};

            my $hit_effector_set = $hits->{$effector_psac};
            next unless $hit_effector_set or
                $potential_pp->{allow_no_effector_matches};
            next if $hit_effector_set && !@$hit_effector_set;

            my $pp_value = $potential_pp->{value};
            # perform post-processing on target:
            if ($postProcessDispatchTable->{$pp_type}) {
                &{$postProcessDispatchTable->{$pp_type}}
                    ($hit_target_set,$hit_effector_set,$pp_value);
                     # call pp subs (from dispatch table)
            } else {
                print STDERR "unknown post-processing key [$pp_type]: ignored";
            }
        }
    }
}

# post process results (inside same entry)
sub postProcess {
    my $hits = shift or return;# input: ref to hash of hits by motif ac...
    my @all_psac_in_seq = sort {
        ($motifRank4PostProcessing->{$b}||0) <=>
            ($motifRank4PostProcessing->{$a}||0)
    } keys %$hits;
    # perform intra motif pp first (repeats...)
    pp_scan($hits, \@all_psac_in_seq, 1);
    # perform inter motif pp
    pp_scan($hits, \@all_psac_in_seq, 0);
}


#############################################################################
# CORE

my $isPPLinkedTo = {};
my $isInSpecifiedPPGroup = {};

sub addPPLinkedMotifs2FetchStruct {
    foreach my $psfile (@prosite_files) {
        # loop through all specified prosite files
        open PSFILE, $psfile or die "Cannot open $psfile : $!";
        my $ac = '';
        while (<PSFILE>) {# read prosite file
            $ac = $1 if m/^AC (\w+);/;
            $ac = '' if m/^\/\//;
            map {
                s/\(.+\)$//;
                $isPPLinkedTo->{$ac}->{$_} = 1, $isPPLinkedTo->{$_}->{$ac} = 1
                if ($_ && $ac ne $_)
            } split(';\s*',$1)
                if $ac && /^PP \/?\w+:\s*([^\r\n]+)/;
        }
    }
    close PSFILE;
    # collect all profile 'linked by pp' to specified psac
    my @followpp_stack = @followpp;
    while (my $psac = pop @followpp_stack) {# for each specified psac
        $specifiedPrositeMotifByAc->{$psac} = 1;
        # add motif ac to 'motif to fetch/use' struct
        foreach my $linked_motifac (keys(%{$isPPLinkedTo->{$psac}})) {
            # add linked profile ac to profile to be hidden at display:
            $hideMotifByPSAC->{$linked_motifac} = 1
                unless grep {$linked_motifac eq $_} @followpp,
                @motifAC_or_userpattern;
            # add linked profile ac to list of profile to be 'pp-followed'
            # (unless already in specifiedPrositeMotifByAc struct):
            push @followpp_stack, $linked_motifac
                unless $specifiedPrositeMotifByAc->{$linked_motifac};
        }
    }
}

# processes PROSITE motif data/info (into global @MotifInfo)
sub processMotif {
    my $ps_entry = shift or return;# prosite entry (string)
    # parse entry
    my ($ac, $id, $type, $de, $pa, $rule, $pdoc, $skipflag,
        $tax, $rep, $sites, $cutoffs, $pps) = parseProsite($ps_entry);
    $ac or $id or die "ERROR: can't parse prosite entry";
    # store post processing data (if any)
    map {
        my $type = $_->{type};
        map {
            my ($effector,$value) = (m/^(\w+\d+)(?:\((.+)\))*/);
            push @{$postProcessingByPSAC->{$ac}},
                {'type'=>$type,'effector'=>$effector,'value'=>$value};
                # print "PP: ac[$ac], type[$type], effector[$effector],
                # value[$value]\n";
        }
        @{$_->{discriminators}}
    } @$pps
        if ($pps && @$pps);
    my $mp_ac;
    if ($type eq "PATTERN" and $opt_miniprofiles
        and $ac=~/PS(\d{5})/ and !$skipflag
        and $opt_format ne "sequence" and $opt_format ne "matchlist") {
        $mp_ac = "MP$1";
        unshift @{$postProcessingByPSAC->{$ac}},
            {'type'=>'EVALUATED_BY', 'effector'=>$mp_ac,
             'value'=>$mp_ac, 'allow_no_effector_matches'=>1 };
    }
    # store implicit 'repeat' post-processing data
    # (RR, R? in Text element of MA CUT_OFF profile lines)
    my $sum=0.0;
    foreach my $cutoff (@$cutoffs) {
        $sum += $cutoff->{N_SCORE} || 0.0 if ($cutoff->{LEVEL} >= -1);
        if ($cutoff->{TEXT} eq "'R?'") {
            # RDM1 pp: promote weak match (add 1 to level)
            # if there is at least one hit at level>=0
            unshift @{$postProcessingByPSAC->{$ac}},
                {'type'=>'PROMOTER_RDM1','effector'=>$ac,'value'=>undef};
        }
        elsif ($cutoff->{TEXT} eq "'RR'") {
            # RDM1 or RDM2 pp: RDM2=
            unshift @{$postProcessingByPSAC->{$ac}},
                {'type'=>'PROMOTER_RDM1_OR_RDM2',
                'effector'=>$ac,'value'=>$sum};
        }
        elsif ($opt_rep_pp_4allprofiles && $cutoff->{TEXT} eq "'!'") {
            # if -a option (opt_rep_pp_4allprofiles):
            # apply RDM1 RDM2 on any hits
            # from a level0 '!' profile (within a seq)!
            unshift @{$postProcessingByPSAC->{$ac}},
                {'type'=>'PROMOTER_RDM1_OR_RDM2_4!REP',
                'effector'=>$ac,'value'=>$sum};
        }
    }
    # sets skipflag, knowfalsepos flags (hashes by ac)
    $SkipFlag{$ac} = 1 if $skipflag;
    my $skip = $skipflag && $opt_skip;
    $skip = !$skipflag if $opt_skiponly;
    if ($ps_entry =~ /^DR.*, T;/m) {
        my $nbfp = 0;
        for my $line ($ps_entry =~ /^DR(.*)/mg) {
            $nbfp += $line =~ s/, F;//g;
        }
        $KnownFalsePos{$ac} = $nbfp;
    }
    # build motif info 'arrays'
    if ($type eq "PATTERN") {# for PATTERN
        $scan_pattern = 1;
        #die "ERROR: can't scan against PATTERN motifs
        #when -w (pfsearch) option is used"
        #if ($opt_pfsearch);
        $scan_pattern = 1;
        push @MotifInfo,
            [$ac, $id, $type, $de,
             prositeToRegexpWrapper($pa, $opt_nongreedy,
                defined $opt_max_x && !$opt_max_x ? 1 : 0),
             $skip];
    }
    elsif ($type eq "MATRIX") {# for profiles
        $scan_profiles = 1;
        if ($specifiedPrositeMotifByAc || $opt_pfsearch) {
            # if not all motifs are used (when some were specified) or
            # if $opt_pfsearch: save data to tmp profile
            if (!$last_profile_tmp_filename || $opt_pfsearch) {
                # open new profile temp
                close PROFILE_TMP if ($last_profile_tmp_filename);
                # close previous one
                $last_profile_tmp_filename=tmpnam();
                # get new temp file name
                open PROFILE_TMP, ">$last_profile_tmp_filename"
                    or die "Cannot open $last_profile_tmp_filename: $!";
            }
            # save profile to temp file
            print PROFILE_TMP $ps_entry
                or die "can't print to $last_profile_tmp_filename: $!\n";
        }
        push @MotifInfo,
            [$ac, $id, $type, $de, undef, $skip, $cutoffs,
            $last_profile_tmp_filename];
    }
    else {
        warn "Unknown prosite entry type $type";
    }
}

sub main {

    ######################## Motif data parsing/processing #####################
    addPPLinkedMotifs2FetchStruct() if !$opt_no_postprocessing && @followpp;
    my $rank = 0;
    for my $psfile ( @prosite_files ) {
    # loop through all specified prosite files
        open PSFILE, $psfile or die "Cannot open $psfile : $!";
        my $ps_entry = ""; my $ac = ""; my $id = ""; my $type = "";
        my $pos = 0;

        PROSITE: while ( <PSFILE> ) {# read prosite file
            $ps_entry .= $_;
            $ac=$1 if !$ac && m/^AC   (\w+);/;
            $id=$1, $type=$2
                if !$id && !$type && m/^ID   (\w+);\s+(\w+)\./;

            if (/^\/\//) {# end of an entry: PROCESS MOTIF DATA
                my $use_motif = 1;
                $use_motif = 0 if (!$id || !$ac);# skip 'bad' entries
                if ($use_motif && $opt_noprofiles && $type eq "MATRIX") {
                # skip profiles if opt_noprofiles
                    $use_motif = 0;
                    die "ERROR: profile(s) specified".
                        " with -p or -f option, but -r".
                        " (do not scan profile) option used"
                        if $specifiedPrositeMotifByAc;
                }
                # skip non profiles if opt_onlyprofiles
                if ($use_motif && $opt_onlyprofiles && $type ne "MATRIX") {
                    $use_motif = 0;
                    die "ERROR: pattern(s) specified with -p or -f option,".
                        " but -m (only scan profile) option used"
                        if $specifiedPrositeMotifByAc;
                }
                $use_motif = 0 if $use_motif && $specifiedPrositeMotifByAc &&
                                    !$specifiedPrositeMotifByAc->{$ac};
                # if AC were specified (PSxxxx format), skip entry that do
                # not match
                delete ($specifiedPrositeMotifByAc->{$ac})
                    if $specifiedPrositeMotifByAc;
                # delete specified AC so that we can check later if they
                # were all found...
                $motifRank4PostProcessing->{$ac} =++ $rank
                    if $use_motif && !$opt_no_postprocessing;
                # ... precedence order for PostProcessing = inverse of motif
                # position in file (lower will be pp first)
                processMotif( $ps_entry ) if $use_motif;
                                    # process motif data into @MotifInfo
                $ps_entry = $ac = $id = $type = "";
                $pos = tell PSFILE;
            }
        }
        close PSFILE;
    }
    close PROFILE_TMP if $last_profile_tmp_filename;
    if ( $specifiedPrositeMotifByAc ) {
    # look if there are unfound specified motif AC
        my @notfound;
        foreach my $ac_id_not_found ( keys( %$specifiedPrositeMotifByAc ) ) {
            push @notfound, $ac_id_not_found;
        }
        die "Prosite entry [@notfound] not found in".
            " specified prosite file(s)\n"
                if @notfound;
    }
    # process -p option user defined patterns
    my $user_ctr = 1;
    for ( @userpat ) {# transform userpat into correct array representation
        my $i = "0" x ( 3 - length( $user_ctr ) ) . $user_ctr++;
        push @MotifInfo, [ "USER$i", undef, "USER", undef,
                        prositeToRegexpWrapper( $_, $opt_nongreedy, 1 ), 0 ];
    }

    ###################### Perform scan (and show results) #####################

    if ( @external_gff_files ) {
        # do not scan: just (post-process and/or convert) gff results
        process_external_gff();
    }
    elsif ( $opt_format ne "matchlist" ) {# SCAN (non matchlist output format)
        # get & show motif hits
        unshift( @ARGV, '-' ) unless @ARGV;
        while ( my $seqfile = shift @ARGV ) {
            # for each specified sequence file: scan
            scanSeqFile($seqfile);
        }
    }
    else {# SCAN (matchlist output format) (fugly!?)
        # FIXME: used!? + seen in comments: only works with sp flat input!?
        showMatchList();
    }

    ############################ exit ############################

    my @tmp_mp_files = values %$files_miniprofiles;
    # delete each temporary miniprofile file
    foreach my $tmp_mp( @tmp_mp_files ) {
        unlink $tmp_mp if $tmp_mp;
    }
    # delete profile temp file(s)
    foreach my $tmp_prf ( map { $_->[7] } @MotifInfo ) {
        unlink $tmp_prf if $tmp_prf;
    }
}
