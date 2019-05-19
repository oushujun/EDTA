#!/usr/local/bin/perl -w
# Written by Yujun Han (yhan@plantbio.uga.edu)
# To group the protein or DNA sequnces and select the exemplars
# 18/Jan/2009
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:b:p:m:s:n:l:o:h:");

$Input    = defined $opt_i ? $opt_i : "";
$BLAST    = defined $opt_b ? $opt_b : "none";
$Protein  = defined $opt_p ? $opt_p : 0;
$Min_P    = defined $opt_m ? $opt_m : 90;
$Min_S    = defined $opt_s ? $opt_s : 80;
$Min_N    = defined $opt_n ? $opt_n : 2;
$No_Left  = defined $opt_l ? $opt_l : "0";
$Output   = defined $opt_o ? $opt_o : $Input.".$Min_P";
$Help     = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------
if($BLAST eq "none") {
	print "self blast...\n";
	if($Protein == 1) {
		system ("_formatdb_ -i $Input -o F -p T\n");
		system ("_blastall_ -i $Input -o $Input.self -e 1e-10 -F F -p blastp -d $Input -v 5000 -b 5000 -m 9\n");
	}else{
		system ("_formatdb_ -i $Input -o F -p F\n");
		system ("_blastall_ -i $Input -o $Input.self -e 1e-10 -F F -p blastn -d $Input -v 5000 -b 5000 -m 9 -G -4 -E -2 -q -3 -r 2\n");		
	}
}

print "Loading blast output ...\n";

%Seq_Group = ();

open(FF, "$Input")||die "Query file not exist\n";
while(<FF>) {
	chomp;
	if(/^>(\S+)/) {
		$Name = $1;
		$Name_Infor{$Name} = $_;
	}else{
		$Name_Seq{$Name} = $_;
		$Name_Len{$Name} = length($_);
	}
}
close(FF);

open(BF, "$Input.self")||die "BLAST result is not exist\n" if($BLAST eq "none");
open(BF, "$BLAST")||die "BLAST result is not exist\n" if($BLAST ne "none");
while(<BF>) {
	chomp;
	next if(/^\#/);

	@Content   = split(/\s+/, $_);
    $Query_Name   = $Content[0];
    $Sbject_Name  = $Content[1];
    $Ident        = $Content[2];
	$Aligned_Len  = $Content[3];

	next if($Query_Name eq $Sbject_Name);
	
    $Query_Len = $Name_Len{$Query_Name};
    $Sbjct_Len = $Name_Len{$Sbject_Name};
	
	if($Query_Len == 0) {
		print "Q0 $Query_Name\n";
	}elsif($Sbjct_Len == 0) {
		print "S0 $Sbject_Name\n";
	}


	if(($Aligned_Len/$Query_Len * 100 >= $Min_P)&&($Aligned_Len/$Sbjct_Len * 100 >= $Min_P)) {
		$Legth_Quality = 1;
	}else{
		$Legth_Quality = 0;
	}

	next if($Legth_Quality == 0);
	
	if($Ident >= $Min_S) {
		$Identy_Quality = 1;
	}else{
		$Identy_Quality = 0;
	}

	next if($Identy_Quality == 0);

	if(defined($Seq_Group{$Query_Name})) {
		if($Seq_Group{$Query_Name} =~ /$Sbject_Name /) {
			next;
		}else{
			$Seq_Group{$Query_Name} .= $Sbject_Name." ";
		}
	}else{
		$Seq_Group{$Query_Name} = $Sbject_Name." ";
	}
}
close(BF);

open(OF, ">$Output")||die "$!\n";
$All_Done = 0;
while($All_Done == 0) {
	$All_Done = 1;
	$Max_Num = 1;
	foreach(keys(%Seq_Group)) {
		$Name = $_;
		next if(defined($Name_Used{$Name}));
		@Members = split(/ /, $Seq_Group{$Name});
		$Num = 1;
		foreach(@Members) {
			$Num ++;
		}
		if($Num > $Max_Num) {
			$Max_Num = $Num;
			$Max_Mem_Name = $Name;
		}
		$All_Done = 0;
	}
	next if($All_Done == 1);

#	print "$Max_Mem_Name\t$Seq_Group{$Max_Mem_Name}\n\n";

	$Max_Mem_Seq = $Name_Seq{$Max_Mem_Name};
	$Name_Used{$Max_Mem_Name} = 1;
	@Max_Mem = split(/ /, $Seq_Group{$Max_Mem_Name});
	foreach(@Max_Mem) {
		$Name_Used{$_} = 1;
	}
	next if($Max_Num < $Min_N);
	print (OF "$Name_Infor{$Max_Mem_Name} $Max_Num\n$Max_Mem_Seq\n");
}

if($No_Left != 0) {
	foreach(keys(%Name_Seq)) {	
		if(!(defined($Name_Used{$_}))) {
			print (OF "$Name_Infor{$_}\n$Name_Seq{$_}\n");
#			print "$_\t$Name_Seq{$_}\n\n";
		}
	}
}
close(OF);
#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :big_element_finder.pl <options> <specification> <default>

	\$Input    = defined \$opt_i ? \$opt_i : "";
	\$BLAST    = defined \$opt_b ? \$opt_b : "none";
	\$Protein  = defined \$opt_p ? \$opt_p : 0;
	\$Min_P    = defined \$opt_m ? \$opt_m : 90;
	\$Min_S    = defined \$opt_s ? \$opt_s : 80;
	\$Min_N    = defined \$opt_n ? \$opt_n : 2;
	\$No_Left  = defined \$opt_l ? \$opt_l : "0";
	\$Output   = defined \$opt_o ? \$opt_o : \$Input.".$Min_P";
	\$Help     = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}