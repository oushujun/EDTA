#!/usr/bin/perl -w
# -----------------------------------------------------
use Getopt::Std;
# -----------------------------------------------------
getopts("q:b:f:m:c:o:");

$Query_Seq = defined $opt_q ? $opt_q : "";
$Blast_Lis = defined $opt_b ? $opt_b : "";
$Flank_Len = defined $opt_f ? $opt_f : 60;
$Max_Edge  = defined $opt_m ? $opt_m : 15;
$Minimal_C = defined $opt_c ? $opt_c : 2;
$Output    = defined $opt_o ? $opt_o : "confirmed_MITEs";

usuage() if((!$Query_Seq)||(!$Blast_Lis));

# -----------------------------------------------------
print "Reading data...\n";
open(QS, "$Query_Seq")||die "$!\n";
while(<QS>) {
	chomp;
	if(/^>(\S+)/) {
		$Name = $1;
	}else{
		$Name_Len{$Name} = length($_);
		$Name_Seq{$Name} = $_;
	}
}
close(QS);

open(BL, "$Blast_Lis")||die "$!\n";
while(<BL>) {
	chomp;
	next if(/^\#/);
	$Line = $_;
	@Content   = split(/\s+/, $Line);
	$Query_Name   = $Content[0];
	$Sbjct_Name   = $Content[1];
	$Query_Start  = $Content[6];
	$Query_Stop   = $Content[7];
	$Sbjct_Start  = $Content[8];
	$Sbjct_Stop   = $Content[9];
	
	$Query_Len = $Name_Len{$Query_Name};
	$Sbjct_Len = $Name_Len{$Sbjct_Name};

	if($Sbjct_Start > $Sbjct_Stop) {
		$Temp = $Sbjct_Start;
		$Sbjct_Start = $Sbjct_Stop;
		$Sbjct_Stop  = $Temp;
	}
	if((abs($Query_Start - $Flank_Len) <= $Max_Edge)&&(abs($Query_Len - $Query_Stop - $Flank_Len) <= $Max_Edge)) {
		if((abs($Sbjct_Start - $Flank_Len) <= $Max_Edge)&&(abs($Sbjct_Len - $Sbjct_Stop - $Flank_Len) <= $Max_Edge)) {
			if(defined($MITE_Found{$Query_Name})) {
				$MITE_Found{$Query_Name} ++;
				$MITE_Found_Lis{$Query_Name} .= $Line."\t"."$Query_Len $Sbjct_Len"."\n";

			}else{
				$MITE_Found{$Query_Name} = 2;
				$Seq_with_Flank = $Name_Seq{$Query_Name};
				$Seq_Without_Flank = substr($Seq_with_Flank, $Flank_Len, $Query_Len - $Flank_Len * 2);

				$MITE_Found_Seq{$Query_Name} = $Seq_Without_Flank;
				$MITE_Found_Lis{$Query_Name} = $Line."\t"."$Query_Len $Sbjct_Len"."\n";
			}
		}
	}
}
close(BL);

open(OF1, ">$Output.fa")||die "$!\n";
open(OF2, ">$Output.lis")||die "$!\n";
foreach(keys(%MITE_Found)) {
	if($MITE_Found{$_} >= $Minimal_C) {
		print (OF1 ">$_\n$MITE_Found_Seq{$_}\n");
		print (OF2 "$MITE_Found_Lis{$_}");
	}
}
close(OF1);
close(OF2);

# -----------------------------------------------------
sub usuage {
    print "\n","Hi,", ' need some help? @_@',"\n";   
    print STDERR <<"    _EOT_";

    Usage: program <options> <specification file> <default>
    
	\$Query_Seq = defined \$opt_q ? \$opt_q : "";
	\$Blast_Lis = defined \$opt_b ? \$opt_b : "";
	\$Flank_Len = defined \$opt_f ? \$opt_f : 60;
	\$Max_Edge  = defined \$opt_m ? \$opt_m : 15;
	\$Output    = defined \$opt_o ? \$opt_o : "confirmed_MITEs";
	\$Minimal_C = defined \$opt_c ? \$opt_c : 2;
	\$Help      = defined \$opt_h ? \$opt_h : "";
           
    _EOT_
    exit(1);
}

