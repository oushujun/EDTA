#!/usr/local/bin/perl -w
# 9/24/2009 pairwise compare the edge region to detect false boundries
# get consensus seq from .aln file and output exact same seq if existed
# Written by prj0401
# 27/5/104
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:t:v:e:m:h:");

$Input   = defined $opt_i ? $opt_i : "";
$Tool    = defined $opt_t ? $opt_t : "none";
$Min_D   = defined $opt_v ? $opt_v : 0.2;	# percentage for consensus
$Edge    = defined $opt_e ? $opt_e : 60;
$Min_Num = defined $opt_m ? $opt_m : 2;
$Help    = defined $opt_h ? $opt_h : "";

usuage() if(((!$Input)&&($Tool eq "none"))||($Help));
#-----------------------------------------------------
if($Tool ne "none") {
	open(MF, ">mass_identify_TE_from_muti.sh")||die($!);
	@Alns = glob("*.aln");
	foreach(@Alns) {
		print (MF "perl $Tool -i $_ -v $Min_D -e $Edge -m $Min_Num\n");
	}
	close(MF);
	system("chmod 770 mass_identify_TE_from_muti.sh\n");
	print "mass_identify_TE_from_muti.sh is ready.\n";
	exit(0);
}

#----------------------------------- read the mutiple alignment file ------------------
open(IF, "$Input")||die "$!\n";
while(<IF>) {
      chomp;
      if(/^>(\S+)/){
         $ID  = $1;
      }else{
         $ID_Seq{$ID} .= $_;
	  }
}
close(IF);

#---------------------------------- get the length of the alignment --------------------
$Total_Len  = length($ID_Seq{$ID});
$Total_Num  = 0;

if($Total_Num == 2) {
	$Min_D -= 0.05;
}

#---------------------------------- get the total number of the alignment --------------
foreach(keys(%ID_Seq)) {
        $Seq              = $ID_Seq{$_};
        $Seq              =~ s/\-//g;
        $Total_Num ++;
}

exit(0) if($Total_Num < $Min_Num);

$Min_Pair_Score    = int($Edge * 2/3);
$Min_Same_Edge_Num = int($Total_Num / 2);

#---------------------------------- count the gap length in head and tail --------------
@Head_Gaps = ();
@Tail_Gaps = ();
foreach(keys(%ID_Seq)) {
	$Seq = $ID_Seq{$_};
	@BPs = split(//, $Seq);

	$Head_Gap = 0;
	$True_Edge = $Edge;
	for($i = 0; $i < $True_Edge; $i ++) {
		if($BPs[$i] eq "-") {
			$Head_Gap ++;
			$True_Edge ++;
		}
	}
	push(@Head_Gaps, $Head_Gap);

	$Tail_Gap = 0;
	$True_Edge = $Edge;
	for($i = $Total_Len - 1; $i > $Total_Len - $True_Edge; $i --) {
		if($BPs[$i] eq "-") {
			$Tail_Gap ++;
			$True_Edge ++;
		}
	}
	push(@Tail_Gaps, $Tail_Gap);
}

$Sum_Head_Gap = 0;
foreach(@Head_Gaps) {
	$Sum_Head_Gap += $_;
}

$Sum_Tail_Gap = 0;
foreach(@Tail_Gaps) {
	$Sum_Tail_Gap += $_;
}

$Head_Gap = int($Sum_Head_Gap / $Total_Num);
$Tail_Gap = int($Sum_Tail_Gap / $Total_Num);

$Head_Edge = $Edge + $Head_Gap;
$Tail_Edge = $Edge + $Tail_Gap;

#---------------------- check whether there are subgroups at the edges ------------
$Head_Sub_Grouped_Num = 0;
$Tail_Sub_Grouped_Num = 0;
foreach(keys(%ID_Seq)) {
	$Key1 = $_;
	$Seq1 = $ID_Seq{$Key1};
	@BPs1 = split(//, $Seq1);
	
	$Head_Edge_Same_Num = 0;
	$Tail_Edge_Same_Num = 0;

	foreach(keys(%ID_Seq)) {
		$Key2 = $_;
		next if($Key1 eq $Key2);
		$Seq2 = $ID_Seq{$Key2};
		@BPs2 = split(//, $Seq2);
		
		#--------------------- Head edge ---------------------
		if($Head_Edge_Same_Num == 0) {
			$Pair_Score = 0;
#			print "\n$Key1\t$Key2\n";
			for($i = 0; $i < $Head_Edge; $i ++) {
				if(($BPs1[$i] ne "-")&&($BPs1[$i] eq $BPs2[$i])) {
					$Pair_Score ++; 
#					print "$BPs1[$i]";
				}
			}
			if($Pair_Score > $Min_Pair_Score) {
				$Head_Edge_Same_Num = 1;
#				print "Head $Pair_Score $Head_Edge\n";
			}
		}
		#--------------------- Tail edge ---------------------
		if($Tail_Edge_Same_Num == 0) {
			$Pair_Score = 0;
			for($i = $Total_Len - 1; $i > $Total_Len - $Tail_Edge - 1; $i --) {
				$Pair_Score ++ if(($BPs1[$i] ne "-")&&($BPs1[$i] eq $BPs2[$i]));
			}

			if($Pair_Score > $Min_Pair_Score) {
				$Tail_Edge_Same_Num = 1;
#				print "Tail $Pair_Score $Tail_Edge \t$Key1\t$Key2\n\n";
			}
		}
	}
	$Head_Sub_Grouped_Num += $Head_Edge_Same_Num;
	$Tail_Sub_Grouped_Num += $Tail_Edge_Same_Num;
}

if(($Head_Sub_Grouped_Num > $Min_Same_Edge_Num)||($Tail_Sub_Grouped_Num > $Min_Same_Edge_Num)) {
	die "$Head_Sub_Grouped_Num/$Total_Num\t$Tail_Sub_Grouped_Num/$Total_Num\n";
}else{
	print "$Head_Sub_Grouped_Num/$Total_Num\t$Tail_Sub_Grouped_Num/$Total_Num\n";
}

#---------------------------------- analyzing the data, comlumn by comlumn -------------
for($i = 0; $i < $Total_Len; $i ++) {          # get the colume Bps
    #------------------------------ find the most abundant bps in each column ----------
	$Bps_Num{"A"} = 0;
    $Bps_Num{"T"} = 0;
    $Bps_Num{"G"} = 0;
    $Bps_Num{"C"} = 0;
    $Bps_Num{"-"} = 0;
    $Bps_Num{"N"} = 0;
    $Bps_Num{"X"} = 0;
	
    foreach(keys(%ID_Seq)) {
            $Bp = substr($ID_Seq{$_}, $i, 1);
            $Bps_Num{$Bp} ++;
    }

    $Max    = 0;
    foreach(keys(%Bps_Num)) {
            if($Bps_Num{$_} > $Max) {
               $Max    = $Bps_Num{$_};
               $Max_Bp = $_;
            }
    }

	#------------------ get the consensus bps --------------------
	$Highest_Per = int($Max / $Total_Num * 1000)/1000;
	$Highest_Pers[$i] = $Highest_Per;
    $Con_Bps[$i] = $Max_Bp;
}

#------------------------------- count how many BPs are conserved at both ends ----
$Head_Con_Per = 0;
$Tail_Con_Per = 0;
$Mid_Con_Per  = 0;

#------------ Head ------------
for($i = 0; $i < $Head_Edge; $i ++) {
	$Head_Con_Per += $Highest_Pers[$i] if($Con_Bps[$i] ne "-");
}

#------------ Tail ------------
for($i = $Total_Len - 1; $i > $Total_Len - $Tail_Edge; $i --) {
	$Tail_Con_Per += $Highest_Pers[$i] if($Con_Bps[$i] ne "-");
}

#------------ Middle ---------------------
$Mid_Len = 0;
$Con_BPs = "";
for($i = $Head_Edge; $i < $Total_Len - $Tail_Edge + 1; $i ++) {
	if($Con_Bps[$i] ne "-") {
		$Mid_Con_Per += $Highest_Pers[$i] ;
		$Con_BPs .= $Con_Bps[$i];
		$Mid_Len ++;
	}
}

#--------------------------
$Aver_Head_Con_Per = int($Head_Con_Per/$Head_Edge * 1000) / 1000;
$Aver_Mid_Con_Per  = int($Mid_Con_Per/$Mid_Len * 1000) / 1000;
$Aver_Tail_Con_Per = int($Tail_Con_Per/$Tail_Edge * 1000) / 1000;

print "$Aver_Head_Con_Per\t$Aver_Mid_Con_Per\t$Aver_Tail_Con_Per\n";

if(($Aver_Mid_Con_Per - $Aver_Head_Con_Per > $Min_D)&&($Aver_Mid_Con_Per - $Aver_Tail_Con_Per > $Min_D)) {
   system "cp $Input $Input.elite\n";
}

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :mutialignment_consensus_getter.pl <options> <specification> <default>

    \$Input   = defined \$opt_i ? \$opt_i : "";
    \$Tool    = defined \$opt_t ? \$opt_t : "none";
    \$Min_D   = defined \$opt_v ? \$opt_v : 0.2;
    \$Edge    = defined \$opt_e ? \$opt_e : 60;
	\$Min_Num = defined \$opt_m ? \$opt_m : 2;
    \$Help    = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}
