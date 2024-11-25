#!/usr/bin/env perl
use warnings;
use strict;

#modified from http://blog.sina.com.cn/s/blog_4ba035220100tkpl.html

my ($dmr, $out, $max_gap); #$max_gap, gaps between this length (bp) will be joined
($dmr, $out, $max_gap) = @ARGV;
$max_gap = 0 unless defined $max_gap;
$out = '' unless defined $out;

die usage() unless @ARGV >= 1;
open DMR,"sort -suV $dmr |" or die "$!";
open OUT, ">$out" or die "$!" if $out ne '';
while (my $line = <DMR>){
	chomp $line;
	next if $line =~ /^\s+$/;
	my ($chr1, $stt1, $end1, $info1) = (split /\s+/, $line, 4); #the first line
	$info1 = '' unless defined $info1;
	PATH:{ #loop the entire file except the first line
		$line = <DMR>;
		if (!$line){
			if ($out ne ''){
				print OUT "$chr1\t$stt1\t$end1\t$info1\n"; #print out the last line
				} else {
				print "$chr1\t$stt1\t$end1\t$info1\n"; #print out the last line
				}
			} 
		else {
			my ($chr2, $stt2, $end2, $info2) = (split /\s+/, $line, 4);
			chomp $info2;
			$info2 = '' unless defined $info2;
			if (($chr1 eq $chr2) && $stt2 >= $stt1 && $stt2 - $end1 <= $max_gap + 1){
				$end1 = $end2 if $end2 > $end1;
				$info1 = $info2;
				redo PATH;
				}
			else {
				if ($out ne ''){
					print OUT "$chr1\t$stt1\t$end1\t$info1\n";
					} else {
					print "$chr1\t$stt1\t$end1\t$info1\n";
					}
				($chr1, $stt1, $end1, $info1) = ($chr2, $stt2, $end2, $info2);
				redo PATH;
				}
			}
		}
	}

sub usage{
    my $die=<<DIE;
    perl *.pl <DMR candidate> <OUT>
DIE
}
