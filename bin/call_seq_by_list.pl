#!/usr/bin/env perl
use warnings;
use strict;

my $usage="\tUsage: perl call_seq_by_list.pl MSU_format_list database(optional, C[custom]) range(optional, itself[default]/up_[int]/down_[int])
eg1. perl call_seq_by_LOC.pl array_list itself >result	##call LOC sequence within the MSU_r7 rice genome
eg2. perl call_seq_by_LOC.pl array_list -C your_database up_2000 >result	##call sequence of upper 2000 bp region in the list, from the provided database

Update history:
	v2.6	read the list line by line and load the genome one chromosome at a time to cut memory use; output unchanged (2026/09/05)
	v2.5	do not output sequence with all Ns (2023/01/30)
	v2.4	output sequences without headers
	v2.3	output a list of entirely excluded sequences
	v2.2	support -ex. Exclude the sequence information provided by list and output the rest. Can use with -cov
	v2.1	support - direction. If the locus position is upside down, eg. Chr1:2..1, it would be treated as negative strand request	Shujun Ou 2014/8/19
	v2.0	Increase speed and reduce memory consumption by introducing the substr function, remove the - direction support	Shujun Ou 2014/8/17
	v1.4	Fix the out of range coordiante	Shujun Ou 2014/07/25
	v1.3	addin custom database support	Shujun Ou 2013/11/27
	v1.2	get rid of the LOC file, use MSU locus format (eg. Chr01:10000..11000) to call sequence in batch	Shujun Ou 2013/11/16
	v1.1	ignore not exist LOC in the list	Shujun Ou 2013/06/21
	v1.0	fix the up/down stream and -/+ strand confusion	Shujun Ou 2013/06/19
	v0.5	birthday of the program skeleton	Shujun Ou 2013/05/17\n";

#array list format: name MSU_position(Chr:from..to, no ":" allowed in Chr)
#example:
#LTR-101308..114181[1]   Chr1:101308..102387
#What_ever_name_u_like   Chr1:113101..114181
#LTR-683241..694400[1]   9311-Chr1:683241..683682

my $position='';
my $range='itself'; ##defalut
my $header=1; #1, output sequence headers (default); 0, no headers
my $length=0; ### default get the LOC seq itself
my $rmvoid=0; #0 for output empty sequences anyways; 1 for output only non-empty sequences. Sequenece with all Ns or Xs are also considered empty.
my $exclude=0; #0 for output sequence specified by list (default); 1 for exclude sequence specified by list
my $coverage=1; #work with $exclude, if the excluded portion is too long (default 1, [0-1]), discard the entire sequence
my $purge=0; #work with $exclude, switch on=1/off=0(default) to clean up aligned region and joint unaligned sequences
my $genome;

my $i=0;
foreach my $para (@ARGV){
	if ($para=~/up_([0-9]+)/i){
		$range=$para;
		$length=$1;
		$position='up';
		}
	elsif ($para=~/down_([0-9]+)/i){
		$range=$para;
		$length=$1;
		$position='down';
		} 
	elsif ($para=~/itself/i){
		$range='itself';
		} ##default range
	$genome=$ARGV[$i+1] if $para=~/^-C$/i;
	$rmvoid=1 if $para=~/^-rmvoid$/i;
	$header=$ARGV[$i+1] if $para=~/^-header$/i;
	$coverage=$ARGV[$i+1] if $para=~/^-cov$/i;
	$purge=$ARGV[$i+1] if $para=~/^-purge$/i;
	$exclude=1 if $para=~/^-ex$/i;
	$i++;
	}

open Exclude, ">$genome.exclude" or die "\n\t\tERROR: Can not create the file $genome.exclude to output excluded seq IDs!\n" if $exclude == 1;
open Genome, "<$genome" or die "\n\t\tERROR: Genome sequence not found, or wrong parameters!\n$usage";

#index the genome FASTA by chromosome instead of loading it whole: %chr_offset stores the byte
#offset of each chromosome's record (the last record of a name wins, as with %genome before),
#@chr_order keeps the records in file order for the -ex whole-genome output
my (%chr_offset, @chr_order);
my $offset = 0;
$/="\n>";
while (<Genome>){
	my $rec_offset = $offset;
	$offset += length $_;
	next if /^>\s?$/;
	chomp;
	tr/>//d;
	s/^\s+//;
	my $pos = index ($_, "\n");
	next if $pos < 0; #no sequence under the header
	my $chr = substr ($_, 0, $pos);
	$chr=~s/\s+$//;
	push @chr_order, [$chr, $rec_offset];
	$chr_offset{$chr}=$rec_offset;
	}
$/="\n";

#load the sequence of one chromosome on demand; only the most recently loaded chromosome is
#kept in memory (list entries within a chromosome are not coordinate-sorted, so random access
#within the current chromosome is still required)
my ($cur_chr, $cur_seq);
sub genome_seq {
	my ($chr) = @_;
	return $cur_seq if defined $cur_chr and $cur_chr eq $chr;
	return unless exists $chr_offset{$chr};
	seek Genome, $chr_offset{$chr}, 0 or return;
	my $rec;
	{
		local $/="\n>";
		defined ($rec = <Genome>) or return;
		chomp $rec;
	}
	$rec=~tr/>//d;
	$rec=~s/^\s+//;
	my $pos = index ($rec, "\n");
	return unless $pos >= 0; #no sequence under the header
	my $chr_rec = substr ($rec, 0, $pos);
	substr ($rec, 0, $pos+1, ''); #drop the header line in place
	$chr_rec=~s/\s+$//;
	$rec=~s/\s+//g;
	($cur_chr, $cur_seq) = ($chr_rec, $rec);
	return $cur_seq;
	}

open List, "sort -k2,2 -suV $ARGV[0] |" or die "\n\tERROR: no LOC list!\n$usage";
my $first = <List>; #read the first line for initialization, then stream the rest

die "Warning: LOC list $ARGV[0] is empty.\n" unless defined $first;
$first = <List> if $first =~ /^\s?$/; #remove the first empty line
my $chr='';
my %chr; #store chr names being worked
$first=~s/^\s+//;
my $chr_pre=$1 if (split /\s+/, $first)[1]=~/(.*):[0-9]+\.\.[0-9]+$/;
my $str=1; #the coordinate of the first bp of a sequence
my $stp; #length of the chr being worked; always assigned before first use
my $seq='';
my $base_num = 0; # count non-N bases in $seq

my $line = $first;
while (defined $line){
	chomp $line;
	next if $line=~/^\s?$/;
	my ($loc, $pos, $strand, $start, $stop);
	$strand="+";
	$line=~s/^\s+//;
	($loc, $pos)=split /\s+/, $line;
	if ($pos=~/^(.*)\:(-?[0-9]+)\.\.(-?[0-9]+)$/){
		($chr, $start, $stop)=($1, $2, $3);
		} else {
	die "ERROR: Can not recognize this MSU position in the list (line: \"$line\")!\n";
		}
	$chr{$chr}=$chr;
	if ($start>$stop){
		($start, $stop)=($stop, $start);
		$strand="-";
		$position="down" if $position eq "up";
		$position="up" if $position eq "down";
		}
	if ($position eq 'up'){
		$stop=$start-1;
		$start=$start-$length;
		}

	if ($position eq 'down'){
		$start=$stop+1;
		$stop=$start+$length-1;
		}

	$start=1 if $start<=0;

	next unless exists $chr_offset{$chr};
#	next if $genome{$chr}=~/^\s+$/;
	next if length (genome_seq($chr)) < 10;
	$stop=length (genome_seq($chr)) if $stop>=length (genome_seq($chr));

	if ($exclude==0){
		$seq=substr (genome_seq($chr), $start-1, $stop-$start+1) if exists $chr_offset{$chr};
		$seq=" " if $seq eq '';

		if ($strand eq "-"){ #if the locus is in neigative strand (-), then get a complementary and reversed strand
			$seq=~tr/tgcaTGCA/acgtACGT/;
			$seq=reverse $seq; ### get a reverse sequence
			($start, $stop)=($stop, $start);
			}
		$base_num = $seq =~ tr/[acgtACGT]/[acgtACGT]/;
		unless ($base_num == 0 and $rmvoid==1){ ###print out target sequence
			print ">$chr:$start..$stop|$loc\n" if $header == 1;
			print "$seq\n";
			}
		}
	my $cov;
	if ($exclude==1){
		if ($chr_pre ne $chr and $chr ne ''){
			if (exists $chr_offset{$chr_pre}){
				$stp=length (genome_seq($chr_pre));
				$seq.=substr (genome_seq($chr_pre), $str-1, $stp-$str+1) if ($str<=$stp and $str!=1);
				if (($stp-length $seq)/$stp >= $coverage){
					print Exclude "$chr_pre\n";
					} elsif ($purge==1){
					print ">$chr_pre\n" if $header==1;
					print "$seq\n";
					} else {
					print ">$chr_pre\n" if $header==1;
					print genome_seq($chr_pre), "\n";
					}
				} else {
				warn "WARNING: $chr_pre is not found in the genome, skipped!\n";
				}

#			if ($purge==1){
#				print ">$chr_pre\n$seq\n" unless ($stp-length $seq)/$stp >= $coverage;
##				print ">$chr_pre|cleanup\n$seq\n" unless ($stp-length $seq)/$stp >= $coverage;
#				} else {
#				print ">$chr_pre\n$genome{$chr_pre}\n" unless ($stp-length $seq)/$stp >= $coverage;
#				}
			$chr_pre=$chr;
			$str=1;
			$stp=length (genome_seq($chr));
			$seq='';
			if ($start>1){
				$stp=$start-1;
				$seq.=substr (genome_seq($chr), $str-1, $stp-$str+1) if exists $chr_offset{$chr};
				}
			$str=$stop+1;
			next;
			}
			if ($start>1){
				$stp=$start-1;
				$seq.=substr (genome_seq($chr), $str-1, $stp-$str+1) if (exists $chr_offset{$chr} and $str<=$stp); #and $str!=1);
				}
			$str=$stop+1;
		}
	}
	continue {
	$line = <List>;
	}
close List;

	if ($exclude==1){
		if (exists $chr_offset{$chr_pre}){
			$stp=length (genome_seq($chr_pre));
			$seq.=substr (genome_seq($chr_pre), $str-1, $stp-$str+1) if ($str<=$stp and $str!=1);
			if (($stp-length $seq)/$stp >= $coverage){
				print Exclude "$chr_pre\n";
				} elsif ($purge==1){
				print ">$chr_pre\n" if $header==1;
				print "$seq\n";
				} else {
				print ">$chr_pre\n" if $header==1;
				print genome_seq($chr_pre), "\n";
				}
			} else {
			warn "WARNING: $chr_pre is not found in the genome, skipped!\n";
			}
#		if ($purge==1){
#			print ">$chr\n$seq\n" unless ($stp-length $seq)/$stp >= $coverage;
##			print ">$chr|cleanup\n$seq\n" unless ($stp-length $seq)/$stp >= $coverage;
#			} else {
#			print ">$chr\n$genome{$chr}\n" unless ($stp-length $seq)/$stp >= $coverage;
#			}
		#output chromosomes not specified in the list; these were iterated in %genome hash
		#order (random per process) before, now they follow the genome file order
		foreach my $chr_rec (@chr_order){
			my ($chr_g, $rec_offset) = @$chr_rec;
			next if $chr_offset{$chr_g} != $rec_offset; #only the last record of a chr counts
			unless (exists $chr{$chr_g}) {
				print ">$chr_g\n" if $header==1;
				print genome_seq($chr_g), "\n";
				}
			}
		}

