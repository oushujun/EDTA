#!/usr/local/bin/perl -w
# Written by vehell
# 11/15/2009
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("t:o:h:");

$Tag    = defined $opt_t ? $opt_t : "TAG";
$Output = defined $opt_o ? $opt_o : $Tag.".lib";
$Help   = defined $opt_h ? $opt_h : "";

usuage() if($Help);
#-----------------------------------------------------
@Family_Files = glob("*.fa");
open(OF, ">$Output")||die"$!\n";
$Tourist_Num = 0;
$Stowaway_Num = 0;
$hAT_Num = 0;
$Mu_Num = 0;
$CACTA_Num = 0;
$Unknown_Num = 0;
$Strange_Num = 0;
foreach(@Family_Files) {
	$File_Name = $_;
	open(IF, "$File_Name")||die"$!\n";
	%Name_Seq = ();
	while(<IF>) {
		chomp;
		if(/^>(\S+)/) {
			$Name = $1;
		}else{
			$Name_Seq{$Name} .= $_;
		}
	}
	close(IF);
	
	if($File_Name =~ /Tourist/) {
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_Tourist_".$Tourist_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$Tourist_Num ++;
		}
	}elsif($File_Name =~ /Stowaway/) {
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_Stowaway_".$Stowaway_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$Stowaway_Num ++;
		}
	}elsif($File_Name =~ /hAT/) {
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_mhAT_".$hAT_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$hAT_Num ++;
		}
	}elsif($File_Name =~ /Mu/) {
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_mMutator_".$Mu_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$Mu_Num ++;
		}
	}elsif($File_Name =~ /CACTA/) {
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_mCACTA_".$CACTA_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$CACTA_Num ++;
		}
	}elsif($File_Name =~ /Unknow/) {
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_mUnknow_".$Unknown_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$Unknown_Num ++;
		}
	}elsif($File_Name =~ /Strange/) {
		print "Strange: $File_Name\n";
		foreach(keys(%Name_Seq)) {
			$Name = $Tag."_mStrange_".$Strange_Num;
			print (OF ">$Name\n$Name_Seq{$_}\n");
			$Strange_Num ++;
		}
	}else{
		print "$File_Name can not be recognized\n";
	}
}
close(OF);

print "Tourist Num: $Tourist_Num\n";
print "Stowaway Num: $Stowaway_Num\n";
print "hAT Num: $hAT_Num\n";
print "Mu Num: $Mu_Num\n";
print "CACTA Num: $CACTA_Num\n";
print "Strange Num: $Strange_Num\n";
print "Unknown Num: $Unknown_Num\n";

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :head_mover.pl <options> <specification> <default>

	\$Tag    = defined \$opt_t ? \$opt_t : "TAG";
	\$Output = defined \$opt_o ? \$opt_o : \$Tag".lib";
	\$Help   = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}
