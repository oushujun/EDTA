#!/usr/bin/perl -w
# 8/17/2010
# add the following parameter because some candidate for MITE_Hunter contains too long low complexitys that makes the programs very slow
# $MaxLen = defined $opt_l ? $opt_l : 35;
# -----------------------------------------------------
use Getopt::Std;
# -----------------------------------------------------
getopts("i:d:p:l:o:h:");

$Input  = defined $opt_i ? $opt_i : "";
$Dusted = defined $opt_d ? $opt_d : "none";
$Per    = defined $opt_p ? $opt_p : 0.2;
$MaxLen = defined $opt_l ? $opt_l : 90;
$Output = defined $opt_o ? $opt_o : "dusted.fa";
$Help   = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------

if($Dusted eq "none") {
	system("_mdust_ $Input > $Input.dusted\n");
	$Dusted = $Input.".dusted";
}

open(DF, "$Dusted")||die"$!\n";
%Undusted_Names = ();
$First = 1;

while(<DF>) {
	chomp;
	if(/^>(\S+)/) {
		$Temp_Name = $1;
		if($First == 1) {
			$Name  = $Temp_Name;
			$First = 0;
			$Seq   = "";
		}else{
			$Len   = length($Seq);
			@BPs   = split(//, $Seq);
			$N_Num = 0;
			$Conti_Len = 0;
			$Max_Conti_Len = 0;
			foreach(@BPs) {
				if($_ eq "N") {
					$N_Num ++;
					$Conti_Len ++;
				}else{
					if($Conti_Len > $Max_Conti_Len) {
						$Max_Conti_Len = $Conti_Len;
					}
					$Conti_Len = 0;
				}
			}
			if($Max_Conti_Len >= $MaxLen) {
#				print "$Name\t$Max_Conti_Len\n";
			}

			if(($N_Num / $Len < $Per)&&($Max_Conti_Len < $MaxLen)) {
				$Undusted_Names{$Name} = 1;
			}
			$Name  = $Temp_Name;
			$Seq   = "";
		}
	}else{
		$Seq .= $_;
	}
}

$Len   = length($Seq);
@BPs   = split(//, $Seq);

$N_Num = 0;
$Conti_Len = 0;
$Max_Conti_Len = 0;
foreach(@BPs) {
	if($_ eq "N") {
		$N_Num ++;
		$Conti_Len ++;
	}else{
		if($Conti_Len > $Max_Conti_Len) {
			$Max_Conti_Len = $Conti_Len;
		}
		$Conti_Len = 0;
	}
}
if(($N_Num / $Len < $Per)&&($Max_Conti_Len < $MaxLen)) {
	$Undusted_Names{$Name} = 1;
}
close(DF);

open(UDF, "$Input")||die"$!\n";
open(FF, ">$Output")||die"$!\n";
while(<UDF>) {
	chomp;
	if(/^>(\S+)/) {
		$Name = $1;
		if(defined($Undusted_Names{$Name})) {
			$Dusted = 0;
			print (FF "\>$Name\n");
		}else{
			$Dusted = 1;
		}	
	}else{
		print (FF "$_\n") if($Dusted == 0);
	}
}
close(FF);
close(UDF);

#----------------------------------------------------
# -----------------------------------------------------
sub usuage {
    print "\n","Hi,", ' need some help? @_@',"\n";   
    print STDERR <<"    _EOT_";

    Usage: program <options> <specification file> <default>
    
	\$Input  = defined \$opt_i ? \$opt_i : "";
	\$Dusted = defined \$opt_d ? \$opt_d : "none";
	\$Per    = defined \$opt_p ? \$opt_p : 0.2;
	\$MaxLen = defined \$opt_l ? \$opt_l : 90;
	\$Output = defined \$opt_o ? \$opt_o : "dusted.fa";
	\$Help   = defined \$opt_h ? \$opt_h : "";
           
    _EOT_
    exit(1);
}