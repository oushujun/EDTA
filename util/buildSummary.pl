#!/usr/bin/env perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) buildSummary
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Summarize the annotations in one or more *.out files.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
# ChangeLog
#
#     $Log$
#     Add support for count_base.pl output to replace -genome_size and 
#     	-seq_count. Shujun Ou (shujun.ou.1@gmail.com) 03/22/2022
#
###############################################################################
#
# To Do:
#

=head1 NAME

buildSummary - Summarize the annotations in one or more *.out files

=head1 SYNOPSIS

  buildSummary [-version] [-minDiv #] [-maxDiv #] [-species <species_name>]
               [-genome <*.tsv or *.2bit>] [-useAbsoluteGenomeSize]
               [-stats <stats from count_base.pl>] <file1.out>[.gz]

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=item -minDiv #

=item -maxDiv #

=item -genome_size #

=item -seq_count #

=item -species <species_name>

=item -stats <stats from count_base.pl>

If the *.out file is manually constructed (i.e., not produced by
RepeatMasker), you need to specify -genome_size and -seq_count so 
that % can be correctly calculated. If you don't want to specify 
these two values manually, you may use this option to privide a 
stats file of the genome produced by count_base.pl in the EDTA 
package. When -stats is provided, the -genome_size and -seq_count 
values will be ignored.

=item -genome <*.tsv or *.2bit>

The default genome size used by this program is obtained using the
sequence ranges in the *.out files.  This is often an overcount, as
many assemblies contain long stretches of Ns and should not count
towards the %masked output.  This option allows the user to specify
either a tab separated file of sequence IDs and their non-ambiguous 
bp size, or a 2bit file containing the sequences themselves.

This is also a handy way to filter out sequences which shouldn't be
tabulated in the output.  Ie. If you have chrUn_* in your *.out file
but want it filtered from the output, simply don't include chrUn_* 
in your *.tsv file or *.2bit file.

NOTE: This program relies on the UCSC program twoBitInfo when you
      use 2bit files.  Please make sure this program is in your path if
      you plan to use this type of file.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2008-2012 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/..";
use Getopt::Long;
use Data::Dumper;
#use Taxonomy;
#use EMBL;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name:  $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name:  $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name:  $';
my $CVSIdTag = '$Id: buildSummary.pl,v 1.36 2017/02/01 21:01:56 rhubley Exp $';
my $Version  = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',                 # print out the version and exit
                    '-useAbsoluteGenomeSize',
                    '-species=s',
                    '-genome=s',
                    '-maxDiv=s',
                    '-minDiv=s',
		    '-stats=s',
                    '-genome_size=s',
                    '-seq_count=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

my $minDiv;
my $maxDiv;
$minDiv = $options{'minDiv'} if ( defined $options{'minDiv'} );
$maxDiv = $options{'maxDiv'} if ( defined $options{'maxDiv'} );

# Total genome length
my $rawSeqLen = 0;
my $gapLen = 0;
my $totalSeqLen = 0;
my $totalSeqNum = 0;
if ( defined $options{'stats'} ){
	my $info = `grep -P '^All' $options{'stats'}`;
	($rawSeqLen, $gapLen, $totalSeqNum) = (split /\s+/, $info)[1,2,4];
	$totalSeqLen = $rawSeqLen - $gapLen;
} else {
	$totalSeqLen = $options{'genome_size'} if ( defined $options{'genome_size'} );
	$totalSeqNum = $options{'seq_count'} if ( defined $options{'seq_count'} );
}

my $taxDB;
my $repDB;
if ( defined $options{'species'} ) {
  $taxDB =
      Taxonomy->new(
                taxonomyDataFile => "$FindBin::Bin/../Libraries/taxonomy.dat" );
  $repDB =
      EMBL->new(
                fileName => "$FindBin::Bin/../Libraries/RepeatMaskerLib.embl" );
}

my %seqUnambigSizes = ();
if ( defined $options{'genome'} ) {
  if ( $options{'genome'} =~ /.*\.2bit/ ) {
    open IN, "twoBitInfo -noNs $options{'genome'} stdout |"
        or die "Could not open $options{'genome'} for reading: $!\n";
    while ( <IN> ) {
      if ( /^(\S+)\s+(\d+)/ ) {
        $seqUnambigSizes{$1} = $2;
      }
    }
    close IN;
  }
  elsif ( $options{'genome'} =~ /.*\.tsv/ ) {
    open IN, "<$options{'genome'}"
        or die "Could not open $options{'genome'} for reading: $!\n";
    while ( <IN> ) {
      if ( /^(\S+)\s+(\d+)/ ) {
        $seqUnambigSizes{$1} = $2;
      }
    }
    close IN;
  }
  else {
    print "\n\nGenome parameter must end in either .tsv or .2bit!\n\n";
    usage();
  }
}

#
# Collect data from out file
#
my %eleStats = ();
my %typeSubTypeCount;
my %typeSubTypeMasked;
my %eleCount;
my %eleMasked;
my %classCount;
my %classCountSeq;
my %classMasked;
my %classMaskedSeq;
my %idsSeen;
my %seqs;
my $lastID        = -1;
my $lastSeqBeg    = 0;
my $lastSeqEnd    = 0;
my $lastSeq1Name  = "";
my $countOUTFiles = 0;
my $warnOnce      = 0;

if ( !@ARGV ) {
  usage();
}

my $file = $ARGV[ 0 ];
if ( $file =~ /.*\.gz/ ) {
  open IN, "gunzip -c $file|"
      or die "Could not open $file using gunzip!\n";
}
else {
  open IN, "<$file"
      or die "Could not open $file!\n";
}
while ( <IN> ) {
  if ( /^\s*\d+\s+\d+\.\d+/ ) {
    my @fields  = split;
    my $seqName = $fields[ 4 ];
    my $div     = $fields[ 1 ];
    my $id      = $fields[ 14 ];
    $id = $lastID + 0.1 if ( !defined $id );

    if ( $seqName ne $lastSeq1Name ) {
      %idsSeen    = ();
      $lastSeqBeg = 0;
      $lastSeqEnd = 0;
      $lastID     = -1;
    }

    if (    defined $minDiv && $div < $minDiv
         || defined $maxDiv && $div > $maxDiv )
    {
      $lastID = $id;
      next;
    }

    # TODO: Remove...this is only relevant to the benchmarks seqs
    #$seqName =~ s/-shuffled//g;
    #$seqName =~ s/-1bp//g;
    #$seqName =~ s/-10bp//g;
    ## TODO: Eventually add an option not to filter.
    if ( defined $options{'genome'} && !exists $seqUnambigSizes{$seqName} ) {
      $lastID = $id;
      if ( !$warnOnce ) {
        warn "INFO: Some results have been filtered because they are on\n"
            . "      sequences not contained in the "
            . $options{'genome'}
            . "\n      file. I.e $seqName\n";
        $warnOnce++;
      }
      next;
    }

    my ( $left ) = ( $fields[ 7 ] =~ /(\d+)/ );
    my $seqLen = $fields[ 6 ] + $left;
    if ( defined $seqs{$seqName} ) {
      if ( $seqs{$seqName} != $seqLen ) {
        warn "Sequence $seqName has changed from size = "
            . $seqs{$seqName}
            . " to size = "
            . $seqLen
            . "\nThis out line is the first instance of the change:\n$_";
      }
    }
    else {
      $seqs{$seqName} = $seqLen;
    }

    my $ele = $fields[ 9 ];

    #
    # Some elements have a "_" appended to the end.
    #
    $ele =~ s/_$//g;

    my $class = $fields[ 10 ];

    #
    #  Break the class down into type/subtype and confidence values
    #
    my $type = "";
    $type = $1 if ( $class =~ /(\S+)\// || $class =~ /(\S+)$/ );
    my $typeConf = 1;
    if ( $type ) {
      $typeConf = 0 if ( $type =~ /\?/ );
      $type =~ s/\?//g;
    }

    my $subType = "";
    $subType = $1 if ( $class =~ /\/(\S+)/ );
    my $subTypeConf = 1;
    if ( $subType ) {
      $subTypeConf = 0 if ( $subType =~ /\?/ );
      $subType =~ s/\?//g;
    }

#print "Class = $class Type = $type ( $typeConf ), SubType = $subType ( $subTypeConf )\n" if ( $DEBUG );

    my $strand   = $fields[ 8 ];
    my $seqBegin = $fields[ 5 ];
    my $seqEnd   = $fields[ 6 ];

    if ( !defined $eleStats{$ele} ) {
      $eleStats{$ele}->{'type'}        = $type;
      $eleStats{$ele}->{'typeConf'}    = $typeConf;
      $eleStats{$ele}->{'subType'}     = $subType;
      $eleStats{$ele}->{'subTypeConf'} = $subTypeConf;
    }

# Talk to Arian about this.  What if an overlap exists between
# two different classes.  Which one should get the larger overlap bp masked?
# also what happens when a repeat is overlapped by another repeat.  Does it get counted?
    my $bpMasked = $seqEnd - $seqBegin + 1;
    if ( $seqBegin <= $lastSeqEnd ) {
      if ( $seqEnd <= $lastSeqEnd ) {

        # Contained inside another repeat
        warn "Fully Contained Repeat: $class =  $lastSeqBeg-$lastSeqEnd "
            . "$seqBegin-$seqEnd $seqName\n"
            if ( $DEBUG );
        $bpMasked = 0;
      }
      else {
        $bpMasked = $seqEnd - $lastSeqEnd;
      }
    }

    if (
      (
        defined $idsSeen{$id} && $idsSeen{$id}->{'class'} ne $class
        && (
          !($idsSeen{$id}->{'class'} =~ /Simple/ && $class =~ /SINE|srpRNA|LINE/
          )
        )
        && (
          !($idsSeen{$id}->{'class'} =~ /SINE|srpRNA|LINE/ && $class =~ /Simple/
          )
        )
        && (
          !($idsSeen{$id}->{'class'} =~ /SINE|srpRNA/ && $class =~ /SINE|srpRNA/
          )
        )
      )
      || ( $id < $lastID - 50 && $id < 3 )
        )
    {

      # Concatenated *.out files...break up
      if ( $DEBUG ) {
        warn "WARNING: Detected concatenated out files...resetting id values\n";
        print STDERR "  --> id: $id prevClass = "
            . $idsSeen{$id}->{'class'}
            . " currClass = $class\n";
      }
      $countOUTFiles++;
      %idsSeen = ();
    }

    if ( !defined $idsSeen{$id}
         || ( $idsSeen{$id} && $idsSeen{$id}->{'ele'} ne $ele ) )
    {

      # This occurs when two fragments have the same id field but
      #   have different repeat names.  Often occurs with LTR names.
      #   The current setup will record each name with it's own count
      #   and bpmasked.
      #warn "WARNING: name/id mismatch $idsSeen{ $id }->{'ele'} != $ele\n"
      #   if ( $idsSeen{ $id } && $idsSeen{ $id }->{'ele'} ne $ele );
      #
      # Each id is counted as one instance of a repeat
      #
      $eleStats{$ele}->{'count'}++;

      if ( $subType ) {
        $typeSubTypeCount{$type}->{$subType}++;
      }
      $typeSubTypeCount{$type}->{'all'}++;
      $classCount{$class}++;
      $eleCount{$ele}++;

      #$classCountSeq{$seqName}->{$class}++ if ( $seqName =~ /chr/ );
      $classCountSeq{$seqName}->{$class}++;
      $idsSeen{$id}->{'class'} = $class;
      $idsSeen{$id}->{'ele'}   = $ele;
    }

    $eleStats{$ele}->{'bpmasked'} += $bpMasked;

    if ( $subType ) {
      $typeSubTypeMasked{$type}->{$subType} += $bpMasked;
    }
    $typeSubTypeMasked{$type}->{'all'} += $bpMasked;
    $classMasked{$class}               += $bpMasked;
    $eleMasked{$ele}                   += $bpMasked;

    #$classMaskedSeq{$seqName}->{$class} += $bpMasked if ( $seqName =~ /chr/ );
    $classMaskedSeq{$seqName}->{$class} += $bpMasked;
    $lastID       = $id;
    $lastSeq1Name = $seqName;
    $lastSeqEnd   = $seqEnd if ( $lastSeqEnd < $seqEnd );
    $lastSeqBeg   = $seqBegin;
  }
}
close IN;

# Total genome length
#my $totalSeqLen = 0;
if ($totalSeqLen == 0){
if ( defined $options{'genome'} ) {
  if ( defined $options{'useAbsoluteGenomeSize'} ) {
    foreach my $seq ( keys( %seqUnambigSizes ) ) {
      if ( defined $seqs{$seq} && $seqUnambigSizes{$seq} > $seqs{$seq} ) {
        die "Error: sequence $seq is larger in the genome than is"
            . " recorded in the results.\n";
      }
      $totalSeqLen += $seqUnambigSizes{$seq};
    }
  }
  else {
    foreach my $seq ( keys( %seqs ) ) {
      if ( $seqUnambigSizes{$seq} > $seqs{$seq} ) {
        die "Error: sequence $seq is larger in the genome than is"
            . " recorded in the results.\n";
      }
      $totalSeqLen += $seqUnambigSizes{$seq};
    }
  }
  %seqs = %seqUnambigSizes;
}
else {
  foreach my $seq ( keys( %seqs ) ) {
    $totalSeqLen += $seqs{$seq};
  }
}
}

#
# Lineage Specific Flag
#
if ( $taxDB ) {
  my %idSpecies = ();
  for ( my $i = 0 ; $i < $repDB->getRecordCount() ; $i++ ) {
    my $record = $repDB->getRecord( $i );
    $idSpecies{ $record->getId() } = [ $record->getRMSpeciesArray() ];
  }
  undef $repDB;

  foreach my $idKey ( sort keys %eleStats ) {

    #  Names like AluSg/q are made up ( ie. not in the database ).  So
    #  to correctly match these we should remove the suffix "/q" first.
    my $id = $idKey;
    $id =~ s/\/.*//;

    #
    # If both of these conditions are true then the repeat is a direct
    # match for the species queried.
    #
    my $isDescendant = 0;
    my $isAncestor   = 0;
    foreach my $name ( @{ $idSpecies{$id} } ) {
      $name =~ s/_/ /g;
      $isDescendant = $taxDB->isA( $name, $options{'species'} );
      $isAncestor = $taxDB->isA( $options{'species'}, $name );
    }
    if ( $isDescendant > 0 ) {
      $eleStats{$idKey}->{'lineageSpecific'} = 1;
    }

   #print "id = $id, isDescendant = $isDescendant, isAncestral = $isAncestor\n";
  }
  undef %idSpecies;
}

#
# Lineage specific repeat stats
#
my $lineageSpecificCount    = 0;
my $lineageSpecificBPMasked = 0;
my $ancestralCount          = 0;
my $ancestralBPMasked       = 0;
if ( $taxDB ) {
  foreach my $ele ( sort keys( %eleCount ) ) {
    if ( $eleStats{$ele}->{'lineageSpecific'} ) {
      $lineageSpecificCount    += $eleStats{$ele}->{'count'};
      $lineageSpecificBPMasked += $eleStats{$ele}->{'bpmasked'};
    }
    else {
      $ancestralCount    += $eleStats{$ele}->{'count'};
      $ancestralBPMasked += $eleStats{$ele}->{'bpmasked'};
    }
  }
}

##
## Class Listing
##
my %typeSubType = ();
foreach my $id ( keys %eleStats ) {
  my $type    = $eleStats{$id}->{'type'};
  my $subType = $eleStats{$id}->{'subType'};
  die "missing type for $id ... <$type>\n" if ( !$type );
  $subType = "none" if ( !$subType );

  if ( !defined $typeSubType{$type}->{$subType} ) {
    $typeSubType{$type}->{$subType}->{'count'}    = 0;
    $typeSubType{$type}->{$subType}->{'bpmasked'} = 0;
  }
  $typeSubType{$type}->{$subType}->{'count'}    += $eleStats{$id}->{'count'};
  $typeSubType{$type}->{$subType}->{'bpmasked'} += $eleStats{$id}->{'bpmasked'};
}
$totalSeqNum = scalar( keys( %seqs ) ) unless $totalSeqNum ne 0;
print "Repeat Classes\n";
print "==============\n";
print "Total Sequences: " . $totalSeqNum . "\n";
print "Total Length: $totalSeqLen bp\n";
if ( $taxDB ) {
  print "Ancestral Repeats: $ancestralCount ( $ancestralBPMasked bp )\n";
  print
"Lineage Specific Repeats: $lineageSpecificCount ( $lineageSpecificBPMasked bp )\n";
}
my $checkCount = 0;
my $checkLen   = 0;
my $checkPerc  = 0;
my %nonInterspersedClasses = (
                               "Low_complexity" => 1,
                               "RNA"            => 1,
                               "Satellite"      => 1,
                               "Simple_repeat"  => 1,
                               "rRNA"           => 1,
                               "scRNA"          => 1,
                               "snRNA"          => 1,
                               "srpRNA"         => 1,
                               "tRNA"           => 1
);
my @skippedClasses = ();
print "Class                  Count        bpMasked    \%masked\n";
print "=====                  =====        ========     =======\n";

foreach my $type ( sort keys( %typeSubType ) ) {
  if ( $nonInterspersedClasses{$type} ) {
    push @skippedClasses, $type;
    next;
  }
  if ( defined $typeSubType{$type}->{'none'} ) {
    print ""
        . sprintf(
                   "%-20s   %-10d   %-10d   %0.2f",
                   $type,
                   $typeSubType{$type}->{'none'}->{'count'},
                   $typeSubType{$type}->{'none'}->{'bpmasked'},
                   (
                     (
                       $typeSubType{$type}->{'none'}->{'bpmasked'} /
                           $totalSeqLen
                     ) * 100
                   )
        )
        . "% \n";
    $checkCount += $typeSubType{$type}->{'none'}->{'count'};
    $checkLen   += $typeSubType{$type}->{'none'}->{'bpmasked'};
    $checkPerc += (
         ( $typeSubType{$type}->{'none'}->{'bpmasked'} / $totalSeqLen ) * 100 );
  }
  else {
    print ""
        . sprintf( "%-20s   %-10s   %-10s   %-5s", $type, "--", "--", "--" )
        . "\n";
  }

  foreach my $subType ( sort keys( %{ $typeSubType{$type} } ) ) {
    next if ( $subType eq "none" );
    print ""
        . sprintf(
                   "    %-16s   %-10d   %-10d   %0.2f",
                   $subType,
                   $typeSubType{$type}->{$subType}->{'count'},
                   $typeSubType{$type}->{$subType}->{'bpmasked'},
                   (
                     (
                       $typeSubType{$type}->{$subType}->{'bpmasked'} /
                           $totalSeqLen
                     ) * 100
                   )
        )
        . "% \n";
    $checkCount += $typeSubType{$type}->{$subType}->{'count'};
    $checkLen   += $typeSubType{$type}->{$subType}->{'bpmasked'};
    $checkPerc += (
       ( $typeSubType{$type}->{$subType}->{'bpmasked'} / $totalSeqLen ) * 100 );
  }
}
print "                      ---------------------------------\n";
print ""
    . sprintf( "%22s %-10d   %-10d   %.2f",
               "   total interspersed",
               $checkCount, $checkLen, $checkPerc )
    . "%\n\n";
foreach my $type ( @skippedClasses ) {
  if ( defined $typeSubType{$type}->{'none'} ) {
    print ""
        . sprintf(
                   "%-20s   %-10d   %-10d   %0.2f",
                   $type,
                   $typeSubType{$type}->{'none'}->{'count'},
                   $typeSubType{$type}->{'none'}->{'bpmasked'},
                   (
                     (
                       $typeSubType{$type}->{'none'}->{'bpmasked'} /
                           $totalSeqLen
                     ) * 100
                   )
        )
        . "% \n";
    $checkCount += $typeSubType{$type}->{'none'}->{'count'};
    $checkLen   += $typeSubType{$type}->{'none'}->{'bpmasked'};
    $checkPerc += (
         ( $typeSubType{$type}->{'none'}->{'bpmasked'} / $totalSeqLen ) * 100 );
  }
  else {
    print ""
        . sprintf( "%-20s   %-10s   %-10s   %-5s", $type, "--", "--", "--" )
        . "\n";
  }

  foreach my $subType ( sort keys( %{ $typeSubType{$type} } ) ) {
    next if ( $subType eq "none" );
    print ""
        . sprintf(
                   "    %-16s   %-10d   %-10d   %0.2f",
                   $subType,
                   $typeSubType{$type}->{$subType}->{'count'},
                   $typeSubType{$type}->{$subType}->{'bpmasked'},
                   (
                     (
                       $typeSubType{$type}->{$subType}->{'bpmasked'} /
                           $totalSeqLen
                     ) * 100
                   )
        )
        . "% \n";
    $checkCount += $typeSubType{$type}->{$subType}->{'count'};
    $checkLen   += $typeSubType{$type}->{$subType}->{'bpmasked'};
    $checkPerc += (
       ( $typeSubType{$type}->{$subType}->{'bpmasked'} / $totalSeqLen ) * 100 );
  }
}
print "---------------------------------------------------------\n";
print ""
    . sprintf( "%-20s   %-10d   %-10d   %.2f",
               "Total", $checkCount, $checkLen, $checkPerc )
    . "%\n\n";

##
## ID Listing
##
print "Repeat Stats\n";
print "============\n";
print "Total Sequences: " . $totalSeqNum . "\n";
print "Total Length: $totalSeqLen bp\n";
if ( $taxDB ) {
  print "Ancestral Repeats: $ancestralCount ( $ancestralBPMasked bp )\n";
  print
"Lineage Specific Repeats: $lineageSpecificCount ( $lineageSpecificBPMasked bp )\n";
}
$checkCount = 0;
$checkLen   = 0;
$checkPerc  = 0;
print "      ID               Count      bpMasked      \%masked\n";
print "================       =====      ========       =======\n";
foreach my $ele ( sort keys( %eleStats ) ) {
  if ( $eleStats{$ele}->{'lineageSpecific'} ) {
    print ""
        . sprintf(
                   "*%16s      %-10d   %-10d   %0.2f",
                   $ele,
                   $eleStats{$ele}->{'count'},
                   $eleStats{$ele}->{'bpmasked'},
                   ( ( $eleStats{$ele}->{'bpmasked'} / $totalSeqLen ) * 100 )
        )
        . "% \n";

  }
  else {
    print ""
        . sprintf(
                   " %16s      %-10d   %-10d   %0.2f",
                   $ele,
                   $eleStats{$ele}->{'count'},
                   $eleStats{$ele}->{'bpmasked'},
                   ( ( $eleStats{$ele}->{'bpmasked'} / $totalSeqLen ) * 100 )
        )
        . "% \n";
  }
  $checkCount += $eleStats{$ele}->{'count'};
  $checkLen   += $eleStats{$ele}->{'bpmasked'};
  $checkPerc  += ( ( $eleStats{$ele}->{'bpmasked'} / $totalSeqLen ) * 100 );
}
print "--------------------------------------------------------\n";
print ""
    . sprintf( "%16s      %-10d   %-10d   %.2f",
               "         ", $checkCount, $checkLen, $checkPerc )
    . "%\n";

print "By Sequence\n";
print "===========\n";

print "Seq            Count      bpMasked\n";
print "=====          =====      ========\n";
foreach my $seq ( keys( %classCountSeq ) ) {
  next if ( $seq =~ /.*_random/ );
  $checkCount = 0;
  $checkLen   = 0;
  foreach my $class ( sort keys( %{ $classCountSeq{$seq} } ) ) {
    $checkCount += $classCountSeq{$seq}->{$class};
    $checkLen   += $classMaskedSeq{$seq}->{$class};
  }
  print ""
      . sprintf( "%16s   %-8d   %-8d", $seq, $checkCount, $checkLen ) . "\n";
}

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my _privateMethod( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##      Something useful.
##
##-------------------------------------------------------------------------##
sub _privateMethod {
  my %parameters = @_;

  print ""
      . ( &caller( 0 ) )[ 0 ] . "::"
      . ( &caller( 0 ) )[ 3 ] . "( "
      . @{ [ %parameters ] }
      . "): Called\n"
      if ( $DEBUG );

}

