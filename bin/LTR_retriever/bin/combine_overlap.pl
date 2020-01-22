#!/usr/bin/env perl -w
use strict;

#modified from http://blog.sina.com.cn/s/blog_4ba035220100tkpl.html
#Shujun Ou (oushujun@msu.edu)

my $max_gap=0; #gaps between this length (bp) will be joined
my ($dmr,$out)=@ARGV;
die usage() unless @ARGV>=2;
my $age=0; #0 means not averaging the age of overlapping entries
$age=1 if defined $ARGV[2];
open DMR,"sort -suV -k1,3 $dmr |" or die "ERROR: $!";
open OUT,">$out" or die "ERROR: $!";
while(my $line=<DMR>){
	next if $line=~/^\s+$/;
	$line=~s/^\s+//;
	my ($chr1,$stt1,$end1,$time1,$time_avg,$i);
	$i=1;
	($chr1,$stt1,$end1,$time1)=(split(/\s+/,$line))[0,1,2,3] if $age==1;
	($chr1,$stt1,$end1)=(split(/\s+/,$line))[0,1,2] if $age==0;
    PATH:{
        $line=<DMR>;
        if(!$line){ #print out the last line of the file
            print OUT "$chr1\t$stt1\t$end1\t$time1\n" if $age==1;
            print OUT "$chr1\t$stt1\t$end1\n" if $age==0;
        }else{
	$line=~s/^\s+//;
        my ($chr2,$stt2,$end2,$time2);
	($chr2,$stt2,$end2,$time2)=(split(/\s+/,$line))[0,1,2,3] if $age==1;
        ($chr2,$stt2,$end2)=(split(/\s+/,$line))[0,1,2] if $age==0;
 	if ($stt2>=$stt1 && $stt2-$end1<=$max_gap+1 && ($chr1 eq $chr2)){ #combine overlapping lines
		$end1=$end2 if $end2>$end1;
		$time1+=$time2 if $age==1;
		$i++;
                redo PATH;
            }else{
		$time_avg=$time1/$i if $age==1;
                print OUT "$chr1\t$stt1\t$end1\t$time_avg\n" if $age==1;
                print OUT "$chr1\t$stt1\t$end1\n" if $age==0;
                ($chr1,$stt1,$end1)=($chr2,$stt2,$end2);
		$time1=$time2 if $age==1;
		$i=1;
                redo PATH;
           }
        }
    }
}

sub usage{
    my $die=<<DIE;
    ERROR: Usage: perl *.pl <DMR candidate> <OUT> [*-age]
DIE
}
