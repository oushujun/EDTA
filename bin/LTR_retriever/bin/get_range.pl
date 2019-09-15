#!usr/bin/perl -w
use strict;

my $usage="to out put some information from the ltrharvest_screen.scn/ltrdigest_tabout.csv/RepeatMasker.out file
	Usage: perl get_range.l [0/1] LTRharvest.scn/ltrdigest_tabout.csv/RepeatMasker.out your_LTRharvest.ltrTE.fa(optional when RepeatMasker.out)
	eg. perl get_range.pl 1 Rice_TIGR7.nmtf.scn.adj Rice_TIGR7.nmtf.ltrTE.pass.list -x
	Author: Shujun Ou 01-01-2014
";

my $version="
	Version: 
		v1.9 Shujun Ou	03-17-2016 Add -fl parameter to output LTR regions with 50 bp two-end-extended coordinates, also change -x to 50bp extension.
		v1.8 Shujun Ou	04-16-2015 Add -g parameter to take whole genome file as input, instead of the *ltrTE.fa file from LTRharvest
		v1.7 Shujun Ou	04-08-2015 Add -int parameter to output the internal region coordinates to the fasta name
		v1.6 Shujun Ou	10-27-2014 Add -L parameter to output the longer LTR region, update several seq ID regex patterns
		v1.5 Shujun Ou	08-13-2014 Add -i parameter to output internal region list
		v1.4 Shujun Ou	08-04-2014 Add -x parameter to output range-extended LTR list
		v1.3 Shujun Ou	06-03-2014 Add on internal length control (default), use -N to disable the control
		v1.2 Shujun Ou	02-19-2014 Support RepeatMasker.out
		v1.1 Shujun Ou  02-12-2014 Support LRTharvest and LTRdigest output, get chr name from LTRharvest.fa file
";

if ($#ARGV<1){die "ERROR: $usage"}
my $LTR=$ARGV[0]; #0 for RepeatMasker analysis, 1 for LTRharvest and LTRdigest analysis
die "ERROR: $usage" unless $LTR=~/^[01]$/;

my $IN=1; #1 for internal length/ratio control (default), 0 for no control
my $min_inlen=100; #minimum internal length
my $max_inlen=15000;
my $min_iLratio=0.05; #dft=0.05, minimum internal/LTR region length ratio
my $max_iLratio=50; #dft=50, maximum internal/LTR region length ratio
my $boundary=0; #0 for no boundary cutting (default). 1 will walk 10bp inside the original range on both ends for boundary alignment.
my $flanking=0; #0 for no flanking output. 1 will output extra 50 bp flanking at both ends for boundary correction, TSD searching, and LTR judging.
my $extend=0; #0 for no region extension (default). 1 will extend 50bp outside the original range for boundary adjustment.
my $longer=0; #0 for no right and left LTR length comparison, 1 will output the longer LTR region only.
my $full=0; #0 will not output full length loci. 1 will output the full range loci.
my $in=0; #0 will not output internal region to the list. 1 will output lLTR, rLTR and internal region in separate lines.
my $int=0; #0 will not output internal region coordinates to the sequence ID. 1 will output start and end position of internal region in the fasta name
my $genome=0; #0 means $ARGV[2] is not genome file; 1 means $ARGV[2] is the genome file with the sequence order same as the LTRharvest sequence number
#my $select=0; #0 will output all regions, 1 will output specific regions ([1], [IN] or [2]) that indicated at the first column

my $k=0;
foreach my $para (@ARGV){
	$boundary=1 if ($para=~/^-b$/i);
	$flanking=1 if ($para=~/^-fl$/i);
	$extend=1 if ($para=~/^-x$/i);
	$longer=1 if ($para=~/^-L$/i);
	$full=1 if ($para=~/^-f$/i);
	$in=1 if ($para=~/^-i$/i);
	$IN=0 if ($para=~/^-N$/i);
	$int=1 if ($para=~/^-int$/i);
	$genome=1 if $para=~/^-g$/i;
	$max_iLratio=$ARGV[$k+1] if $para=~/^-max_ratio$/i;
	$k++;
	}

open TBL, "<$ARGV[1]" or die "ERROR: $!";
open LTRlist, ">$ARGV[1].list" unless $extend==1;
open Extend, ">$ARGV[1].extend" if $extend==1;
open Full, ">$ARGV[1].full" if $full==1;

my %chr;
if ($LTR==1 && $genome==0){
	open FA, "<$ARGV[2]" or die "ERROR: $!";
	while (<FA>){
		s/>//;
		$chr{"$2..$3"}=$1 if (/^(\S+).*\[([0-9]+),([0-9]+)\]/); #eg: >9311_chr01 (dbseq-nr 0) [101308,114181]
		$chr{"$2..$3"}=$1 if (/^(\S+)_[0-9]+.*\[([0-9]+),([0-9]+)\]/); #eg: >gi.478805111.gb.AQOG01030080.1_1 (dbseq-nr 173) [4033,8637]
		$chr{"$2..$3"}=$1 if (/^(\S+)\:([0-9]+)\.\.([0-9]+)\|([0-9]+)\.\.([0-9]+)/); #eg: >Chr1:106522..118080|106502..118100
		$chr{"$2..$3"}=$1 if (/^(\S+)\|([0-9]+)\.\.([0-9]+)/);#eg: >Chr1|106522..118080
		$chr{"$2..$3"}=$1 if (/^(\S+)\:([0-9]+)\.\.([0-9]+)\|(\S+)/);#eg: >Chr1:106522..118080|Chr1
		$chr{"$2..$3"}=$1 if (/^(\S+):([0-9]+)..([0-9]+)\[[12]\]/); #eg: >gi.478805265.gb.AQOG01029926.1:10426..15413[1]
#print "$id\n";
		$chr{"$2..$3"}=$1 if (/(\S+)\:([0-9]+)\.\.([0-9]+)/); #eg: gi.478789307.gb.AQOG01045884.1:56716..59758     pass    motif:AAAG      TSD:TGAAG
		$chr{"$2..$3"}=$1 if (/^(\S+)\:([0-9]+)\.\.([0-9]+)\|/);#eg: >Chr1:106522..118080|
#print "$1\n";
		}
	}

if ($LTR==1 && $genome==1){
	my @id=`grep \\> $ARGV[2]`;
	my @rev_id = reverse @id;
	my $i=0;
	foreach (@id, @rev_id){
		chomp;
		s/>//g;
		s/\s+//g;
		$chr{$i}=$_;
		$i++;
		}
	}

while (<TBL>){
	if (/^#/){next}
	if (/perc/ or /score/){next;}
	if (/^\s+$/){next}
	my ($element_start, $element_end, $element_length, $sequence, $lLTR_start, $lLTR_end, $lLTR_length, $rLTR_start, $rLTR_end, $rLTR_length, $lTSD_start, $lTSD_end, $lTSD_motif, $rTSD_start, $rTSD_end, $rTSD_motif, $PPT_start, $PPT_end, $PPT_motif, $PPT_strand, $PPT_offset, $PBS_start, $PBS_end, $PBS_strand, $tRNA, $tRNA_motif, $PBS_offset, $tRNA_offset, $PBS_tRNA_edist, $Protein_domain_hits, $similarity, $seq_ID, $chr);

##This is for RepeatMasker out file analysis
##    SW   perc perc perc  query     position in query              matching                    repeat          position in repeat
## score   div. del. ins.  sequence  begin    end          (left)   repeat                      class/family begin   end    (left)    ID
##121374    0.3  0.0  0.0  Chr1       2991444  3005061 (38467215) C Os2091in-Osr28_AP002539     LTR/Gypsy       (0)  13617       1  304
if ($LTR==0){
	s/^\s+//;
	(undef, $similarity, undef, undef, $seq_ID, $element_start, $element_end, undef)=split;
	$seq_ID=(split /_/, $seq_ID)[0];
	print LTRlist "$seq_ID:$element_start..$element_end\t$seq_ID:$element_start..$element_end\n";
	}

if ($LTR==1){

##This is for LTRharvest result analysis
##s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr
	s/^\s+//;
	($element_start,  $element_end,  $element_length, $lLTR_start,  $lLTR_end,  $lLTR_length, $rLTR_start,  $rLTR_end,  $rLTR_length, $similarity, $seq_ID, $chr)=(split /\s+/, $_);
#34 4594 4561 34 291 258 4335 4594 260 0.962     + CTCAC 29..33, 4595..4599 TG,TG,CA,CA
#start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len similarity seqid chr direction TSD lTSD rTSD motif superfamily family age(ya)

	$chr=$chr{"$element_start..$element_end"} if $genome==0;
#	if ($genome==1 and !$chr) {
	if ($genome==1) {
		if (exists $chr{$seq_ID}){
			$chr=$chr{$seq_ID};
			} else {
			$chr=$seq_ID;
			}
		}

	next unless defined $lLTR_length and defined $rLTR_length;
#	print "ERROR: could not recognize the following line:\n\t$_" unless defined $lLTR_length and defined $rLTR_length;
	my $long="NA";
	if ($lLTR_length>=$rLTR_length){
		$long="left";
		} else {
		$long="right"
		}

##to obtain the boundary sequences
	my ($int_str, $int_end)=($lLTR_end+1,$rLTR_start-1);
	if ($boundary==1){
		$lLTR_start+=10;
		$lLTR_end-=10;
		$rLTR_start+=10;
		$rLTR_end-=10;
		}

##for internal region length control
	my $mark=1;
	if ($IN==1){
		my $in_len=$rLTR_start-$lLTR_end;
		my $inLTR_ratio=$in_len/(($lLTR_length+$rLTR_length)/2);
		if ($in_len<$min_inlen or $in_len>$max_inlen or $inLTR_ratio<$min_iLratio or $inLTR_ratio>$max_iLratio){
			$mark=0;
			}
		}

	if ($extend==1 or $flanking==1){
		$lLTR_start-=50; #directly extend the coordinate for 50 bp
		$lLTR_end+=50;
		$rLTR_start-=50;
		$rLTR_end+=50;
		print Extend "$chr:$element_start..$element_end\t$chr:$lLTR_start..$rLTR_end\n" if defined $chr;
		}

	if ($full==1){
		print Full "$chr:$element_start..$element_end\t$chr:$element_start..$element_end\n" if defined $chr;
		}

##This is for LTRdigest result analysis 
#	($element_start,  $element_end,  $element_length,  $sequence,  $lLTR_start,  $lLTR_end,  $lLTR_length,  $rLTR_start,  $rLTR_end,  $rLTR_length,  $lTSD_start,  $lTSD_end,  $lTSD_motif,  $rTSD_start,  $rTSD_end,  $rTSD_motif,  $PPT_start,  $PPT_end,  $PPT_motif,  $PPT_strand,  $PPT_offset,  $PBS_start,  $PBS_end,  $PBS_strand,  $tRNA,  $tRNA_motif,  $PBS_offset,  $tRNA_offset,  $PBS_tRNA_edist,  $Protein_domain_hits)=split (/\t/,$_);

	if ($mark==1 and defined $chr and $extend!=1 and $int==0){
		print LTRlist "$chr:$element_start..$element_end\[1]\t$chr:$lLTR_start..$lLTR_end\n" if ($longer==0 or ($longer==1 && $long eq "left"));
		print LTRlist "$chr:$element_start..$element_end\[2]\t$chr:$rLTR_start..$rLTR_end\n" if ($longer==0 or ($longer==1 && $long eq "right"));
		print LTRlist "$chr:$element_start..$element_end\[IN]\t$chr:$int_str..$int_end\n" if $in==1;
		}
	
	if ($mark==1 and defined $chr and $extend!=1 and $int==1){
		print LTRlist "$chr:$element_start..$element_end\[1]($int_str..$int_end-$similarity)\t$chr:$lLTR_start..$lLTR_end\n";
		print LTRlist "$chr:$element_start..$element_end\[2]($int_str..$int_end-$similarity)\t$chr:$rLTR_start..$rLTR_end\n";
		}	
	}
}

close FA if $LTR==1;
close Extend if $extend==1;
close Full if $full==1;
close LTRlist unless $extend==1;
close TBL;

