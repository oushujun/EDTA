#!/usr/local/bin/perl -w
# 8/13/2010
# add: 
# split the given fasta file into small segments for blast, etc.
# Written by vehell
# 30/6/106
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:w:s:P:o:h:");

$Input  = defined $opt_i ? $opt_i : "";
$Window = defined $opt_w ? $opt_w : 12000;
$Step   = defined $opt_s ? $opt_s : 4000;
$Per    = defined $opt_P ? $opt_P : 1;	# increase this number (at most 1) if you are pacient enough
$Output = defined $opt_o ? $opt_o : "$Input.split";
$Help   = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------
open(IF, "$Input")||die"$!\n";
while(<IF>) {
	chomp;
	if(/^>(\S+)/) {
		$Name = $1;
	}else{
		$Name_Seq{$Name} .= $_;
	}
}
close(IF);

open(OF, ">$Output")||die"$!\n";
foreach(keys(%Name_Seq)) {
		$Name = $_;
		$Seq  = $Name_Seq{$Name};
		$Len  = length($Seq);
		if($Len > $Window) {
			$First = 1;
			for($i = 0; $i < $Len - $Window; $i += $Step) {
				if($First == 1) {
					if($i + $Window < $Len) {
						$Win = substr($Seq, $i, $Window);
					}else{
						$Win = substr($Seq, $i);
					}
					print(OF ">$Name.$i\n$Win\n");
					$First = 0;
					$Sub_Per = 0;
				}else{
					$Sub_Per += $Per;

					if($Sub_Per < 1) {
						next;
					}else{
						if($i + $Window < $Len) {
							$Win = substr($Seq, $i, $Window);
						}else{
							$Win = substr($Seq, $i);
						}
						print(OF ">$Name.$i\n$Win\n");
						$Sub_Per = 0;
					}
				}
			}
		}else{
			print(OF ">$Name.0\n$Seq\n");
		}
}
close(OF);

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :fasta_windows_maker.pl <options> <specification> <default>

	\$Input  = defined \$opt_i ? \$opt_i : "";
	\$Window = defined \$opt_w ? \$opt_w : 12000;
	\$Step   = defined \$opt_s ? \$opt_s : 4000;
	\$Per    = defined \$opt_P ? \$opt_P : 1;	# increase this number (at most 1) if you are pacient enough
	\$Output = defined \$opt_o ? \$opt_o : "\$Input.split";
	\$Help   = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}
