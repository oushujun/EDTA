#!usr/bin/perl -w
use strict;

#modified from http://blog.sina.com.cn/s/blog_4ba035220100tkpl.html

my $max_gap=0; #gaps between this length (bp) will be joined
my ($dmr,$out)=@ARGV;
die usage() unless @ARGV==2;
open DMR,"sort -t \$'\t' -k1,1 -k2n,2 -k3n,3 $dmr |" or die "$!";
open OUT,">$out" or die "$!";
while(my $line=<DMR>){
	next if $line=~/^\s+$/;
	my ($chr1,$stt1,$end1)=(split(/\s+/,$line))[0,1,2];
    PATH:{
        $line=<DMR>;
        if(!$line){
            print OUT "$chr1\t$stt1\t$end1\n";
        }else{
        my ($chr2,$stt2,$end2)=(split(/\s+/,$line))[0,1,2];
 	if ($stt2>=$stt1 && $stt2-$end1<=$max_gap+1 && ($chr1 eq $chr2)){
		$end1=$end2 if $end2>$end1;
                redo PATH;
            }else{
                print OUT "$chr1\t$stt1\t$end1\n";
                ($chr1,$stt1,$end1)=($chr2,$stt2,$end2);
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
