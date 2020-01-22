##Shujun Ou
##Usage: perl Six-frame_translate.pl sequence.nt.fa > sequence.aa.fa

#!/usr/bin/env perl -w
use strict;

my %codon2aa = qw(
	TAC  Y  TAT  Y  TAA  *  TAG  *  TGC  C  TGT  C  TGA  *  TGG  W
	CTA  L  CTC  L  CTG  L  CTT  L  CCA  P  CCC  P  CCG  P  CCT  P
	TCA  S  TCC  S  TCG  S  TCT  S  TTC  F  TTT  F  TTA  L  TTG  L
	CAC  H  CAT  H  CAA  Q  CAG  Q  CGA  R  CGC  R  CGG  R  CGT  R
	ATA  I  ATC  I  ATT  I  ATG  M  ACA  T  ACC  T  ACG  T  ACT  T
	AAC  N  AAT  N  AAA  K  AAG  K  AGC  S  AGT  S  AGA  R  AGG  R
	GTA  V  GTC  V  GTG  V  GTT  V  GCA  A  GCC  A  GCG  A  GCT  A
	GAC  D  GAT  D  GAA  E  GAG  E  GGA  G  GGC  G  GGG  G  GGT  G
	);

$/="\n>";
while (<>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$id=~s/\s+//g;
	$seq=~s/\s+//g;
	$seq=uc $seq;

#reverse complement the ori sequence
	my $seq_rev=reverse $seq;
	$seq_rev=~tr/atcgATCG/tagcTAGC/;

#six-frame translate the sequence
	my $pro1=seq2aa($seq);
	$seq=~s/.//;
	my $pro2=seq2aa($seq);
	$seq=~s/.//;
	my $pro3=seq2aa($seq);

	my $rev_pro1=seq2aa($seq_rev);
	$seq_rev=~s/.//;
	my $rev_pro2=seq2aa($seq_rev);
	$seq_rev=~s/.//;
	my $rev_pro3=seq2aa($seq_rev);
	print ">$id|aa1\n$pro1\n>$id|aa2\n$pro2\n>$id|aa3\n$pro3\n>$id|rev_aa1\n$rev_pro1\n>$id|rev_aa2\n$rev_pro2\n>$id|rev_aa3\n$rev_pro3\n";
	}
$/="\n";

sub seq2aa {
	my $seq=$_[0];
	my $aa='';
	while ($seq =~ s/(...)//){
		my $codon='';
		if (exists $codon2aa{$1}){
			$codon=$codon2aa{$1};
			} else {
			$codon='X';
			}
		$aa.=$codon;
		}
	$aa=~s/[^A-Z*]+//g;
	return $aa;
	}

