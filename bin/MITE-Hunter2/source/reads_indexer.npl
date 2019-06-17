#!/usr/local/bin/perl -w
# build the reads index file 
# Written by prj03
# 29/10/102
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:o:h:");

$Input_Reads = defined $opt_i ? $opt_i : "";
$Output      = defined $opt_o ? $opt_o : $Input_Reads.".index";
$Help        = defined $opt_h ? $opt_h : "";

usuage() if((!$Input_Reads)||($Help));
#-----------------------------------------------------
open(IF, "$Input_Reads")||die "$!\n";
$Pre_Len  = 0;
while(<IF>) {
      if(/^>(\S+)/) {
		$Name = $1;
		if(defined($Name_Loc{$Name})) {
			die("duplicated name: $Name\n");
		}else{
			$Name_Loc{$Name} = $Pre_Len + length($_);
		}
		$Pre_Len += length($_);
	  }else{
         $Pre_Len += length($_);
       }
}
close(IF);

open(ID, ">$Output\n")||die "$!\n";
foreach(keys(%Name_Loc)) {
        print (ID "$_ $Name_Loc{$_}\n");
}
close(ID);

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :reads_indexer <options> <specification> <default>

    \$Input_Reads = defined \$opt_i ? \$opt_i : "";
    \$Output      = defined \$opt_o ? \$opt_o : \$Input_Reads.".index";
    \$Help        = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}
