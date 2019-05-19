#!/usr/local/bin/perl -w
# Written by vehell
# 5/7/106
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:n:o:h:");

$Input  = defined $opt_i ? $opt_i : "";
$Number = defined $opt_n ? $opt_n : 5;
$Output = defined $opt_o ? $opt_o : "group";
$Help   = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------
$Num = 0;
open(IF, "$Input")||die "$!\n";
while(<IF>) {
	if(/^>/) {
		$Num ++;
	}
}
close(IF);

$Group_Num = int($Num / $Number) + 1;

print "$Group_Num/$Num\n";

$Num = 0;
$i = 0;
open(IF, "$Input")||die "$!\n";
while(<IF>) {
	if(/^>/) {
		if($Num % $Group_Num == 0) {
			print "$Num $Group_Num\n";
			$i ++;
			if($Num > 0) {
				close(OF);
			}
			open(OF, ">$Output.$i")||die "$!\n";	
		}
		
		$Num ++;
		print(OF "$_"); 
	}else{
		print (OF "$_");
	}
}
close(IF);
close(OF);

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :head_mover.pl <options> <specification> <default>

	\$Input  = defined \$opt_i ? \$opt_i : "";
	\$Number = defined \$opt_n ? \$opt_n : 5;
	\$Output = defined \$opt_o ? \$opt_o : "group";
	\$Help   = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}
