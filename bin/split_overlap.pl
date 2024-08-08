#!/usr/bin/env perl
use warnings;
use strict;
# Shujun Ou (07/29/2020 shujun.ou.1@gmail.com)

#modified from http://blog.sina.com.cn/s/blog_4ba035220100tkpl.html

my ($input, $output, $min_len); #after split, homology-based entries < $min_len will be removed
($input, $output, $min_len) = @ARGV;
$min_len = 50 unless defined $min_len; #default 50 bp
die usage() unless @ARGV >= 2;

# iteratively split the input gff3 file
my $iter = 1;
my ($in, $out) = ($input, "$input.iter$iter");
for (my $i=0; $i<$iter; $i++){
	open IN, "sort -suV $in |" or die "$!";
	open OUT, ">$out" or die "$!";

while (my $line = <IN>){
	chomp $line;
	next if $line =~ /^\s+$/;
	my ($chr1, $stt1, $end1, $info1) = (split /\s+/, $line, 4); #the first line
	my $method1 = (split /\s+/, $info1)[2];
	$method1 = 'NA' unless defined $method1;
	($stt1, $end1) = ($end1, $stt1) if $stt1 > $end1;
	my $len1 = $end1 - $stt1 + 1;
	next if ($len1 < $min_len and $method1 eq 'homology') or $len1 <= 2;
	PATH:{ #loop the entire file except the first line
		$line = <IN>;
		if (!$line){
			print OUT "$chr1\t$stt1\t$end1\t$info1\n" unless ($len1 < $min_len and $method1 eq 'homology') or $len1 <= 2; #print out the last line
			} 
		else {
			chomp $line;
			my ($chr2, $stt2, $end2, $info2) = (split /\s+/, $line, 4); #the next line
			my $method2 = (split /\s+/, $info2)[2];
			$method2 = 'NA' unless defined $method2;
			($stt2, $end2) = ($end2, $stt2) if $stt2 > $end2;
			my $len2 = $end2 - $stt2 + 1;
			next if ($len2 < $min_len and $method1 eq 'homology') or $len2 <= 2;
			if (($chr1 eq $chr2) && $stt2 <= $end1 && $end2 > $end1){
				my $keep = &compare($len1, $method1, $len2, $method2);
				if ($keep eq 'keep1'){
					print OUT "$chr1\t$stt1\t$end1\t$info1\n" unless ($len1 < $min_len and $method1 eq 'homology') or $len1 <= 2;
					$stt2 = $end1 + 1;
					($stt2, $end2) = ($end2, $stt2) if $stt2 > $end2;
					$len2 = $end2 - $stt2 + 1;
					($stt1, $end1, $info1, $len1, $method1) = ($stt2, $end2, $info2, $len2, $method2);
					redo PATH;
					} else {
					$end1 = $stt2 - 1;
					($stt1, $end1) = ($end1, $stt1) if $stt1 > $end1;
					$len1 = $end1 - $stt1 + 1;
					print OUT "$chr1\t$stt1\t$end1\t$info1\n" unless ($len1 < $min_len and $method1 eq 'homology') or $len1 <= 2;
					($stt1, $end1, $info1, $len1, $method1) = ($stt2, $end2, $info2, $len2, $method2);
					redo PATH;
					}
				}
			elsif (($chr1 eq $chr2) && $stt2 <= $end1 && $end2 <= $end1){ #the next line is nested within
				my $end_temp = $stt2 - 1;
				$len1 = $end_temp - $stt1 + 1;
				print OUT "$chr1\t$stt1\t$end_temp\t$info1\n" unless ($len1 < $min_len and $method1 eq 'homology') or $len1 <= 2; #first section
				print OUT "$chr2\t$stt2\t$end2\t$info2\n" unless ($len2 < $min_len and $method2 eq 'homology') or $len2 <= 2; #nested section
				($stt1, $len1) = ($end2 + 1, $end1 - $end2 + 1); #$chr1, $end1, $info1, $method1 unchanged
				redo PATH;
				}
			else {
				print OUT "$chr1\t$stt1\t$end1\t$info1\n" unless ($len1 < $min_len and $method1 eq 'homology') or $len1 <= 2;
				($chr1, $stt1, $end1, $info1, $len1, $method1) = ($chr2, $stt2, $end2, $info2, $len2, $method2);
				redo PATH;
				}
			}
		}
	}
	close OUT;
	close IN;
	my $old_stat = `wc -l "$in"`;
	my $new_stat = `wc -l "$out"`;
	$old_stat = (split /\s+/, $old_stat)[0];
	$new_stat = (split /\s+/, $new_stat)[0];
	if ($old_stat == $new_stat){
		last;
		} else {
		$iter++;
		$in = $out;
		$out = "$input.iter$iter";
		}
	}
`mv $out $output`;


# determine which annotation to keep
sub compare (){
	my ($len1, $method1, $len2, $method2) = ($_[0], $_[1], $_[2], $_[3]);
	my $result = 'NA';
	$method1 = 'homology' if $method1 !~ /homology|structural/i;
	$method2 = 'homology' if $method2 !~ /homology|structural/i;
	if ($method1 eq $method2){ #keep the longest if annotated by the same method
		if ($len1 >= $len2){ #keep the first element if has the same length
			$result = 'keep1';
			}
		else {
			$result = 'keep2';
			}
	} elsif ($method1 =~ /structural/i){ #keep structural-based method disregard length
		$result = 'keep1';
		}
	else {
		$result = 'keep2';
		}
	return $result;
	}

# usage
sub usage{
    my $die=<<DIE;
Split overlapping annotations. Prefer to keep structural-based entries and/or longer entries.
    perl split_overlap.pl <BED IN> <BED OUT>
DIE
}

