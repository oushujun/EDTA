#!/usr/local/bin/perl -w
# Written by Yujun Han (yhan@plantbio.uga.edu)
# 5/10/2010
# add the function to identify pair TEs
# use 80-80-80 (A unified classification system for eukaryotic transposable elements)
# 10/31/2009 To classify sequences into families, based on BLAST output
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:b:p:e:o:n:a:t:B:S:m:s:h:");

$Input      = defined $opt_i ? $opt_i : "";
$BLAST      = defined $opt_b ? $opt_b : "none";
$Protein    = defined $opt_p ? $opt_p : 0;
$Max_Evalue = defined $opt_e ? $opt_e : "1e-5";
$Output     = defined $opt_o ? $opt_o : "output";
$CPU_Num    = defined $opt_n ? $opt_n : 5;
$Min_Algn   = defined $opt_a ? $opt_a : 80;
$Min_Ident  = defined $opt_t ? $opt_t : 80;
$Big_Per	= defined $opt_B ? $opt_B : 0.7;
$Small_Per	= defined $opt_S ? $opt_S : 0.9;
$Min_Camira = defined $opt_m ? $opt_m : 250;
$Step       = defined $opt_s ? $opt_s : 1;	#	1: just group (after manually check) 2: group into known families
$Help       = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------
if($BLAST eq "none") {
	print "self blast...\n";
	if($Protein == 1) {
		system ("_formatdb_ -i $Input -o F -p T\n");
		system ("_blastall_ -i $Input -o $Input.self -e $Max_Evalue -p blastp -d $Input -v 5000 -b 5000 -m 9 -a $CPU_Num\n");
	}else{
		system ("_formatdb_ -i $Input -o F -p F\n");
		system ("_blastall_ -i $Input -o $Input.self -e $Max_Evalue -p blastn -d $Input -v 5000 -b 5000 -m 9 -G -4 -E -2 -q -3 -r 2\n");		
	}
}

print "Loading blast output ...\n";

open(FF, "$Input")||die "Query file not exist\n";
while(<FF>) {
	chomp;
	if(/^>(\S+)/) {
		$Name = $1;
		$Name_Infor{$1} = $_;
	}else{
		$Name_Seq{$Name} .= $_;
	}
}
close(FF);

#---------------- calculate TE length ------------------
foreach(keys(%Name_Seq)) {
	$Name = $_;
	$Name_Len{$Name} = length($Name_Seq{$Name});
}

#---------------- looking for paired TEs -----------------------
open(BF, "$Input.self")||die "BLAST result is not exist\n" if($BLAST eq "none");
open(BF, "$BLAST")||die "BLAST result is not exist\n" if($BLAST ne "none");
while(<BF>) {
	chomp;
	next if(/^\#/);

	@Content    = split(/\s+/, $_);
    $Query_Name = $Content[0];
    $Sbjct_Name = $Content[1];
	$Identity   = $Content[2];
	$Align_Len  = $Content[3];

	$Query_Len = $Name_Len{$Query_Name};
	$Sbjct_Len = $Name_Len{$Sbjct_Name};

	next if($Query_Name eq $Sbjct_Name);
	next if($Query_Len < $Min_Camira);
	if(($Align_Len / $Query_Len < $Big_Per)&&($Align_Len / $Sbjct_Len > $Small_Per)) {
		$Paired_TEs{$Query_Name} .= $Sbjct_Name." ";
	}
}
close(BF);

open(OF, ">$Output.paired")||die"$!\n";
foreach(keys(%Paired_TEs)) {
	print (OF "\>$_\t$Paired_TEs{$_}\n$Name_Seq{$_}\n");
}
close(OF);
#---------------------------------------------------------------

%Group_TEs = ();
%TE_Group  = (); 
$Group_Num = 1;

open(BF, "$Input.self")||die "BLAST result is not exist\n" if($BLAST eq "none");
open(BF, "$BLAST")||die "BLAST result is not exist\n" if($BLAST ne "none");
while(<BF>) {
	chomp;
	next if(/^\#/);

	@Content    = split(/\s+/, $_);
    $Query_Name = $Content[0];
    $Sbjct_Name = $Content[1];
	$Identity   = $Content[2];
	$Align_Len  = $Content[3];

	next if($Align_Len < $Min_Algn);

	next if($Identity < $Min_Ident);

	next if($Query_Name eq $Sbjct_Name);

	next if(defined($Paired_TEs{$Query_Name}));

	if(defined($TE_Group{$Query_Name})) {
		if(defined($TE_Group{$Sbjct_Name})) {
			if($TE_Group{$Query_Name} == $TE_Group{$Sbjct_Name}) {
				next;
			}else{
				@Sbjct_Group_TEs = split(/ /, $Group_TEs{$TE_Group{$Sbjct_Name}});
				$Group_TEs{$TE_Group{$Sbjct_Name}} = "__none__";
				foreach(@Sbjct_Group_TEs) {
					if($Group_TEs{$TE_Group{$Query_Name}} =~ /$_ /) {
						print "Wrong1\n";
					}else{
						$Group_TEs{$TE_Group{$Query_Name}} .= $_." ";
						$TE_Group{$_} = $TE_Group{$Query_Name};
					}
				}
				$TE_Group{$Sbjct_Name} = $TE_Group{$Query_Name};
			}
		}else{
			if($Group_TEs{$TE_Group{$Query_Name}} =~ /$Sbjct_Name /) {
				print "Wrong2\n";
			}else{
				$Group_TEs{$TE_Group{$Query_Name}} .= $Sbjct_Name." ";
				$TE_Group{$Sbjct_Name} = $TE_Group{$Query_Name};
			}
		}
		
	}elsif(defined($TE_Group{$Sbjct_Name})) {
		if($Group_TEs{$TE_Group{$Sbjct_Name}} =~ /$Query_Name /) {
			print "Wrong3\n";
		}else{
			$Group_TEs{$TE_Group{$Sbjct_Name}} .= $Query_Name." ";
			$TE_Group{$Query_Name} = $TE_Group{$Sbjct_Name};
		}
	}else{
		$Group_TEs{$Group_Num} = $Query_Name." ".$Sbjct_Name." ";
		$TE_Group{$Query_Name} = $Group_Num;
		$TE_Group{$Sbjct_Name} = $Group_Num;
		$Group_Num ++;
	}
}
close(BF);

%TE_Grouped   = ();
$TE_Total_Num = 0;

$Camira_Group_Num = 1;
$mPIFHarbinger_Group_Num = 1;
$mMariner_Group_Num = 1;
$mhAT_Group_Num = 1;
$mCACTA_Group_Num = 1;
$mMutator_Group_Num = 1;
$Unknown_Group_Num = 1;

$Group_Num    = 1;

foreach(keys(%Group_TEs)) {
	next if($Group_TEs{$_} eq "__none__");
	@TE_Names = split(/ /, $Group_TEs{$_});
	if($Step == 1) {
		$File_Name = $Output."_".$Group_Num.".fa";
		open(OF, ">$File_Name")||die"$!\n";
	}else{
		$mPIFHarbinger  = 0;
		$mMariner = 0;
		$mhAT     = 0;
		$mMutator = 0;
		$mCACTA   = 0;
		
		foreach(@TE_Names) {
			if(/mPIFHarbinger/) {
				$mPIFHarbinger = 1;
			}
			if(/mMariner/){
				$mMariner = 1;
			}
			if(/mhAT/) {
				$mhAT = 1;	
			}
			if(/Mutator/){
				$mMutator = 1;
			}
			if(/CACTA/){
				$mCACTA = 1;
			}
		}
		if($mPIFHarbinger + $mMariner + $mhAT + $mMutator + $mCACTA > 1) {
			$File_Name = $Output."_"."camira"."_".$Camira_Group_Num.".fa";
			$Camira_Group_Num ++;
		}elsif($mPIFHarbinger + $mMariner + $mhAT + $mMutator + $mCACTA == 1) {
			if($mPIFHarbinger == 1) {
				$File_Name = $Output."_mPIFHarbinger_".$mPIFHarbinger_Group_Num.".fa";
				$mPIFHarbinger_Group_Num ++;
			}
			if($mMariner == 1){
				$File_Name = $Output."_mMariner_".$mMariner_Group_Num.".fa";
				$mMariner_Group_Num ++;
			}
			if($mhAT == 1) {
				$File_Name = $Output."_mhAT_".$mhAT_Group_Num.".fa";
				$mhAT_Group_Num ++;
			}
			if($mMutator == 1){
				$File_Name = $Output."_mMutator_".$mMutator_Group_Num.".fa";
				$mMutator_Group_Num ++;
			}
			if($mCACTA == 1){
				$File_Name = $Output."_mCACTA_".$mCACTA_Group_Num.".fa";
				$mCACTA_Group_Num ++;
			}
		}else{
			$File_Name = $Output."_Unknown_".$Unknown_Group_Num.".fa";
			$Unknown_Group_Num ++;
		}

		open(OF, ">$File_Name")||die"$!\n";
	}

	foreach(@TE_Names) {
		$TE_Grouped{$_} = 1;
		if((defined($Name_Infor{$_}))&&(defined($Name_Seq{$_}))) {
			if($Step == 1) {
				print (OF "$Name_Infor{$_}\n$Name_Seq{$_}\n");
			}else{
				print (OF "$Name_Infor{$_}\n$Name_Seq{$_}\n");
			}
		}else{
			print "here $_\n";
		}
		$TE_Total_Num ++;
	}
	close(OF);
	$Group_Num ++;
}

#-------------------------- output singlet ------------------
if($Step == 1) {
	$File_Name = $Output."_singlet.fa";
	open(OF, ">$File_Name")||die"$!\n";
	foreach(keys(%Name_Infor)) {
		if(defined($TE_Grouped{$_})) {
			next;
		}else{
			print (OF "$Name_Infor{$_}\n$Name_Seq{$_}\n");
			$TE_Total_Num ++;
		}
	}
}else{
	foreach(keys(%Name_Infor)) {
		$Key = $_;
		if(defined($TE_Grouped{$_})) {
			next;
		}else{
			if($Key =~ /mPIFHarbinger/) {
				$File_Name = $Output."_mPIFHarbinger_".$mPIFHarbinger_Group_Num.".fa";
				$mPIFHarbinger_Group_Num ++;
			}elsif($Key =~ /mMariner/){
				$File_Name = $Output."_mMariner_".$mMariner_Group_Num.".fa";
				$mMariner_Group_Num ++;
			}elsif($Key =~ /mhAT/) {
				$File_Name = $Output."_mhAT_".$mhAT_Group_Num.".fa";
				$mhAT_Group_Num ++;
			}elsif($Key =~ /mMutator/){
				$File_Name = $Output."_mMutator_".$mMutator_Group_Num.".fa";
				$mMutator_Group_Num ++;
			}elsif($Key =~ /mCACTA/){
				$File_Name = $Output."_mCACTA_".$mCACTA_Group_Num.".fa";
				$mCACTA_Group_Num ++;
			}else{
				$File_Name = $Output."_Unknown_".$Unknown_Group_Num.".fa";
				$Unknown_Group_Num ++;
			}

			open(OF, ">$File_Name")||die"$!\n";
			print (OF "$Name_Infor{$_}\n$Name_Seq{$_}\n");
			$TE_Total_Num ++;
		}
	}
}
close(OF);

print "Total $TE_Total_Num Sequences were classified into $Group_Num groups\n";
#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :big_element_finder.pl <options> <specification> <default>

	\$Input      = defined \$opt_i ? \$opt_i : "";
	\$BLAST      = defined \$opt_b ? \$opt_b : "none";
	\$Protein    = defined \$opt_p ? \$opt_p : 0;
	\$Max_Evalue = defined \$opt_e ? \$opt_e : "1e-5";
	\$Output     = defined \$opt_o ? \$opt_o : "output";
	\$CPU_Num    = defined \$opt_n ? \$opt_n : 5;
	\$Min_Algn   = defined \$opt_a ? \$opt_a : 80;
	\$Min_Ident  = defined \$opt_t ? \$opt_t : 80;
	\$Big_Per	 = defined \$opt_B ? \$opt_B : 0.7;
	\$Small_Per	 = defined \$opt_S ? \$opt_S : 0.9;
	\$Min_Camira = defined \$opt_m ? \$opt_m : 250;
	\$Step       = defined \$opt_s ? \$opt_s : 1;	#	1: just group (after manually check) 2: group into known families
	\$Help       = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}