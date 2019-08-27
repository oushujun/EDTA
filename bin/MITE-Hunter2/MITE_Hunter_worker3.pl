#!/usr/local/bin/perl -w
# 8/19/2010
# blast -F T
# 6/16/2010
# add trucated function for checking long sequnces. A short sequence using only 200 bp from the head and the tail will be used.
# 12/11/2009
# add the function of filtering out TE flanks with additional sequences; two parameters were added: $Min_Addi and $Min_Addi_Per
# 12/10/2009
# 1) Add two parameters: $ReBlast_Num and $Long_Muscle to encrease sensitivity of BLAST and reduce running time for muscle
# 2) Set max intron = 10. or else it will introduce big gaps into the mutiple alignment file
# 3) Set default $Max_Output = 35
# 12/10/2009
# set blast to high resolution mode
# 11/26/2009
# get consensus even didn't find the TIR/TSD structure, making it easier for manual check
# 11/10/2009 change E value from 0.000001 to 1e-10
# 10/28/2009
# get the consensus, find TSD and TIR and do classification
# for Altix only, otherwise the path needs to be changed
#-----------------------------------------------------
use Getopt::Std;
use FindBin;
#-----------------------------------------------------
my $script_path = $FindBin::Bin;

getopts("i:D:e:M:m:I:n:E:t:s:T:p:P:b:r:l:f:v:a:A:c:G:o:h:");

$Input        = defined $opt_i ? $opt_i : "";
$Database     = defined $opt_D ? $opt_D : "/iob_home/srwlab/vehell/scratch/Data_Room/Maize/maize_pseudo.seq";
$Max_Evalue   = defined $opt_e ? $opt_e : "1e-10"; 
$Min_Pro      = defined $opt_M ? $opt_M : 0.7;		# for homologs
$Min_Copy     = defined $opt_m ? $opt_m : 3;
$ID_Tag       = defined $opt_I ? $opt_I : "ZM";
$Max_Output   = defined $opt_n ? $opt_n : 35;
$Edge_Len     = defined $opt_E ? $opt_E : 60;
$TIR_Len      = defined $opt_t ? $opt_t : 8;
$TIR_Mismatch = defined $opt_s ? $opt_s : 2; 
$Max_TSD      = defined $opt_T ? $opt_T : 10;
$Min_TSD_Per  = defined $opt_p ? $opt_p : 0.5;		# for TSD
$Min_Mul_Per  = defined $opt_P ? $opt_P : 0.7;		# for multiple alignment
$Max_Boundry  = defined $opt_b ? $opt_b : 6;
$ReBlast_Num  = defined $opt_r ? $opt_r : 10;
$Long_Muscle  = defined $opt_l ? $opt_l : 900;
$Fixed_TSDs   = defined $opt_f ? $opt_f : "TA_";
$Min_Diverge  = defined $opt_v ? $opt_v : 0.25;
$Min_Addi     = defined $opt_a ? $opt_a : 25;		# minimal additional sequence length
$Min_Addi_Per = defined $opt_A ? $opt_A : 0.8;		# minimal percentage for identifying additional sequence
$Truncated    = defined $opt_c ? $opt_c : 0;
$Max_Gap      = defined $opt_G ? $opt_G : 10;
$Output       = defined $opt_o ? $opt_o : "output.fa";
$Help         = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------
if($Fixed_TSDs ne "none") {
	@Fixed_Seqs = split(/_/, $Fixed_TSDs);
	foreach(@Fixed_Seqs) {
		push(@Fixed_Lens, length($_));
	}
}

#--------------------------------------------------------
open(OC, ">$Output")||die"$!\n";

#------------------------------------------------
open(IF, $Input)||die"$!\n";
while(<IF>) {
	$Line = $_;
	if(/^>(\S+)/) {
		$Name = $1;
	}else{
		$Line =~ s/\s+$//;
		$Name_Seq{$Name} .= $Line;
	}
}
close(IF);

foreach(keys(%Name_Seq)) {
	$Name = $_;
	$Raw_Seq  = $Name_Seq{$Name};

	if($Truncated == 1) {
		$Raw_Seq  = $Name_Seq{$Name};
		if(length($Raw_Seq) > 400) {
			$Head_Seq = substr($Raw_Seq, 0, 200);
			$Tail_Seq = substr($Raw_Seq, length($Raw_Seq) - 200, 200);
			$Max_Gap = length($Raw_Seq)- 350;
			$Raw_Seq = $Head_Seq.$Tail_Seq;
		}
	}

	print "$Name\n";
	
	open(QF, ">temp_query.fa")||die"$!";
	print(QF ">$Name\n$Raw_Seq\n");
	close(QF);
	
	print "BLAST...\n";

	#------------------------------ sometimes using too sensitive search will result in low quality alignment that is difficult to detect TIRs and TSDs
	system("blastall -i temp_query.fa -d $Database -e $Max_Evalue -p blastn -o temp_blast.out\n");
	
	print "Searching for homologs by PHI ...\n";
	open(TF, "temp_blast.out")||die"$!\n";
	$No_Hits_Found = 0;
	while(<TF>) {
		if(/No hits found/) {
			$No_Hits_Found = 1;
			print "$Name has no homologs!\n";
			last;
		}
	}
	close(TF);

	if($No_Hits_Found == 0) {
		print "Find homologs...\n";
		system("perl $script_path/IDTE_PHI.pl -i temp_blast.out -q temp_query.fa -D $Database -o $Name -M $Min_Pro -d $Max_Gap -I $ID_Tag -n $Max_Output -e $Max_Evalue -l $Edge_Len -r $Edge_Len\n"); #shujun
		if(-e "$Name.flank") {
			$Seq_Num = 0;
			open(NF, "$Name.flank")||die"2 $!\n";
			while(<NF>) {
				if(/^\>/) {
					$Seq_Num ++;
				}
			}
			close(NF);
			
			#----------------------------------- use high resolution blast parameters if homologs number is smaller than ReBlast_Num ---------
			if($Seq_Num < $ReBlast_Num) {
				print "$Seq_Num < $ReBlast_Num blast again using -G -4 -E -2 -q -3 -r 2 ...\n";
				system("blastall -i temp_query.fa -d $Database -e $Max_Evalue -p blastn -o temp_blast.out -G -4 -E -2 -q -3 -r 2\n");
				print "Find homologs...\n";
				system("perl $script_path/IDTE_PHI.pl -i temp_blast.out -q temp_query.fa -D $Database -o $Name -M $Min_Pro -d $Max_Gap -I $ID_Tag -n $Max_Output -e $Max_Evalue -l $Edge_Len -r $Edge_Len\n");
				$Seq_Num = 0;
				open(NF, "$Name.flank")||die"2 $!\n";
				while(<NF>) {
					if(/^\>/) {
						$Seq_Num ++;
					}
				}
				close(NF);
			}

			#------------------------------ mutiple alignment, identify TEs and find consensus ---------------------
			if($Seq_Num >= $Min_Copy) {
				print "multiple alignment...\n";
				if(length($Raw_Seq) < $Long_Muscle) {
					system("muscle -in $Name.flank -out $Name.aln 2>/dev/null"); #shujun slienced muscle output
				}else{
					print "(length < $Long_Muscle run muscle with -maxiters 1 -diags)\n";
					system("muscle -in $Name.flank -out $Name.aln -maxiters 1 -diags 2>/dev/null"); #shujun slienced muscle output
				}

				
				#-------------------------------- filter out sequences with additional sequences ---------------------
				%Temp_ID_Seq  = ();
				%Temp_Bps_Num = ();
				@Temp_Con_Bps = ();
				@Temp_Highest_Pers = ();
				%Temp_Addi_Names = ();
				%Temp_Name_Seq = ();

				open(IF, "$Name.aln")||die "$Name.aln not exits\n";	#------------ read the mutiple alignment file ------------------
				while(<IF>) {
					chomp;
				    if(/^>(\S+)/){
						$ID  = $1;
					}else{
						$Temp_ID_Seq{$ID} .= $_;
					}
				}
				close(IF);
				
				$Total_Len  = length($Temp_ID_Seq{$ID});	#---------- get the length of the alignment --------------------
				$Total_Num  = 0;
				
				foreach(keys(%Temp_ID_Seq)) {	#---------- get the total number of the alignment --------------
			        $Seq = $Temp_ID_Seq{$_};
			        $Seq =~ s/\-//g;
			        $Total_Num ++;
				}

				for($i = 0; $i < $Total_Len; $i ++) {          # get the colume Bps
					$Temp_Bps_Num{"A"} = 0;
					$Temp_Bps_Num{"T"} = 0;
					$Temp_Bps_Num{"G"} = 0;
					$Temp_Bps_Num{"C"} = 0;
					$Temp_Bps_Num{"-"} = 0;
					$Temp_Bps_Num{"N"} = 0;
					$Temp_Bps_Num{"X"} = 0;
		
					foreach(keys(%Temp_ID_Seq)) {
						$Bp = substr($Temp_ID_Seq{$_}, $i, 1);
						$Temp_Bps_Num{$Bp} ++;
					}

					$Max = 0;
					foreach(keys(%Temp_Bps_Num)) {
						if($Temp_Bps_Num{$_} > $Max) {
							$Max    = $Temp_Bps_Num{$_};
							$Max_Bp = $_;
						}
					}

					$Highest_Per = int($Max / $Total_Num * 1000)/1000;	#------------------ get the consensus bps --------------------
					$Temp_Highest_Pers[$i] = $Highest_Per;
					$Temp_Con_Bps[$i] = $Max_Bp;
				}

				$Addi_Flag = 0;
				foreach(keys(%Temp_ID_Seq)) {	#-------------------- find the sequence that have long additional sequences -------------
					$Addi = 0;
					$Temp_Name = $_;
					$Len  = 0;
					for($i = 0; $i < $Total_Len; $i ++) {
				        $Bp = substr($Temp_ID_Seq{$Temp_Name}, $i, 1);
						if(($Bp ne $Temp_Con_Bps[$i])&&($Temp_Highest_Pers[$i] > $Min_Addi_Per)) {
							if($Addi == 0) {
								$Len = 1;
								$Addi = 1;
							}else{
								$Len ++;
							}
						}else{
							if($Len > $Min_Addi) {
								$Temp_Addi_Names{$Temp_Name} = 1;
								$Addi_Flag = 1;
								last;
							}else{
								$Len = 0;
								$Addi = 0;
							}
						}
				    }
				}
				
				if($Addi_Flag == 1) {		# ------------- filter out bad sequecnes then do multiple alignment again -------------
					open(TF, "$Name.flank")||die"$Name.flank not exits\n";
					while(<TF>) {
						chomp;
						if(/^>(\S+)/) {
							$Temp_Name = $1;
						}else{
							$Temp_Name_Seq{$Temp_Name} .= $_;
						}
					}
					close(TF);
					
					open(TOF, ">$Name.flank")||die"$!\n";
					$Add_Filted_Num = 0;
					foreach(keys(%Temp_Name_Seq)) {
						if(defined($Temp_Addi_Names{$_})) {
#							print "$_ was filtered out\n";
							$Add_Filted_Num ++;
							next;
						}else{
							print (TOF ">$_\n$Temp_Name_Seq{$_}\n");
						}
					}
					close(TOF);
					print "$Add_Filted_Num with additional sequences were filtered out\n";

					if(length($Raw_Seq) < $Long_Muscle) {
						system("muscle -in $Name.flank -out $Name.aln 2>/dev/null"); #shujun slienced muscle output
					}else{
						print "(length < $Long_Muscle run muscle with -maxiters 1 -diags)\n";
						system("muscle -in $Name.flank -out $Name.aln -maxiters 1 -diags 2>/dev/null"); #shujun slienced muscle output
					}
				}

				#----------------------------------------- identify whether it is a TE or not ----------------------------
				system("perl $script_path/identify_TE_from_mutialign.pl -i $Name.aln -e $Edge_Len -v $Min_Diverge -m $Min_Copy\n"); #shujun
				
				system("rm $Name.aln\n");

				if(-e "$Name.aln.elite") {
					$CC = consensus_classifier();
					if($CC == 1) {
						print(OC "\>$Name TSD Len:$Max_TSD_Len\n$Con_Seq\n");
					}else{
						print(OC "\>$Name Unknow\n$Con_Seq\n");
					}
				}
			}
			system("rm $Name.flank\n");
		}
	}
	system("rm temp*\n");
}
close(OC);


#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :tblastn_copy_finder.pl <options> <specification> <default>

	\$Input        = defined \$opt_i ? \$opt_i : "";
	\$Database     = defined \$opt_D ? \$opt_D : "/scratch/vehell/Data_Room/Maize/maize_pseudo.seq";
	\$Max_Evalue   = defined \$opt_e ? \$opt_e : "1e-10"; 
	\$Min_Pro      = defined \$opt_M ? \$opt_M : 0.7;
	\$Min_Copy     = defined \$opt_m ? \$opt_m : 3;
	\$ID_Tag       = defined \$opt_I ? \$opt_I : "ZM";
	\$Max_Output   = defined \$opt_n ? \$opt_n : 35;
	\$Edge_Len     = defined \$opt_E ? \$opt_E : 60;
	\$TIR_Len      = defined \$opt_t ? \$opt_t : 8;
	\$TIR_Mismatch = defined \$opt_s ? \$opt_s : 2; 
	\$Max_TSD      = defined \$opt_T ? \$opt_T : 10;
	\$Min_TSD_Per  = defined \$opt_p ? \$opt_p : 0.5;
	\$Min_Mul_Per  = defined \$opt_P ? \$opt_P : 0.7;		# for multiple alignment
	\$Max_Boundry  = defined \$opt_b ? \$opt_b : 6;
	\$ReBlast_Num  = defined \$opt_r ? \$opt_r : 10;
	\$Long_Muscle  = defined \$opt_l ? \$opt_l : 900;
	\$Fixed_TSDs   = defined \$opt_f ? \$opt_f : "none";	# TA_
	\$Min_Diverge  = defined \$opt_v ? \$opt_v : 0.25;
	\$Min_Addi     = defined \$opt_a ? \$opt_a : 25;		# minimal additional sequence length
	\$Min_Addi_Per = defined \$opt_A ? \$opt_A : 0.8;		# minimal percentage for identifying additional sequence
	\$Truncated    = defined \$opt_c ? \$opt_c : 0;
	\$Output       = defined \$opt_o ? \$opt_o : "output.fa";

	\$Help         = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}

#-----------------------------------------------------
sub consensus_classifier {
	%ID_Seq = ();
	#----------------------------------- read the mutiple alignment file ------------------
	open(IF, "$Name.aln.elite")||die "$Name.aln.elite $!\n";
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

	#---------------------------------- get the total number of the alignment --------------
	$Total_Num  = 0;

	foreach(keys(%ID_Seq)) {
	    $Seq = $ID_Seq{$_};
        $Seq =~ s/\-//g;
        $Total_Num ++;
	}

	#---------------------------------- count the gap length in head and tail --------------
	@Head_Gaps = ();
	@Tail_Gaps = ();
	foreach(keys(%ID_Seq)) {
		$Seq = $ID_Seq{$_};
		@BPs = split(//, $Seq);

		$Head_Gap = 0;
		$True_Edge = $Edge_Len;
		for($i = 0; $i < $True_Edge; $i ++) {
			if($BPs[$i] eq "-") {
				$Head_Gap ++;
				$True_Edge ++;
			}
		}
		push(@Head_Gaps, $Head_Gap);

		$Tail_Gap = 0;
		$True_Edge = $Edge_Len;
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

	$Head_Edge = $Edge_Len + $Head_Gap;
	$Tail_Edge = $Edge_Len + $Tail_Gap;
	
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

	    $Max     = 0;
		$Max_Num = 0; 
		foreach(keys(%Bps_Num)) {
            if($Bps_Num{$_} > $Max) {
               $Max     = $Bps_Num{$_};
               $Max_Bp  = $_;
			   $Max_Num = $Bps_Num{$_}; 
            }
		}

		#------------------ get the consensus bps --------------------
		$Con_Bps[$i]    = $Max_Bp;
		$Con_BP_Num[$i] = $Max_Num;
	}

	#---------------------------------- find the putative start sites -------------
	@Head_Start_Locs = ();
	for($i = 20; $i < $Head_Edge + $Edge_Len/2; $i ++) {          # get the colume Bps
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
			   $Max_BP = $_;
            }
		}

		next if($Max_BP eq "-");

		if($Max / $Total_Num >= $Min_Mul_Per) {
			push(@Head_Start_Locs, $i);
		}
	}
	
	#---------------------- find the most likely 5' TIR start site -------------
	$Head_Start_Nums = 0;
	@Head_Start_Locs_New = ();
	for($i = 0; $i < @Head_Start_Locs - 1; $i ++) {
		$Gap_Num = 0;
		$Con_Num = 0;
		for($j = $Head_Start_Locs[$i]; $j < $Head_Start_Locs[$i] + 10 + $Gap_Num; $j ++) {
			if($Con_Bps[$j] eq "-") {
				$Gap_Num ++;
			}elsif($Con_BP_Num[$j] / $Total_Num > $Min_Mul_Per) {
				$Con_Num ++;
			}else{
				next;
			}
		}
		if($Con_Num > 5) {
			$Head_Start_Nums ++;
			push(@Head_Start_Locs_New, $Head_Start_Locs[$i]);
		}
		last if($Head_Start_Nums >= $Max_Boundry);
	}
	@Head_Start_Locs = @Head_Start_Locs_New;

#	print "Head: @Head_Start_Locs\n";

	#---------------------------------- find the putative tail sites -------------
	@Tail_Stop_Locs = ();
	for($i = $Total_Len - 20; $i > $Total_Len - $Tail_Edge - $Edge_Len/2; $i --) {          # get the colume Bps
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
			   $Max_BP = $_;
            }
		}
#		print "$i\t$Max_BP\t$Max / $Total_Num > $Min_Mul_Per\n";
		next if($Max_BP eq "-");
#		print "$i\t$Max / $Total_Num > $Min_Mul_Per?\n";

		if($Max / $Total_Num >= $Min_Mul_Per) {
			push(@Tail_Stop_Locs, $i);
		}
	}

#	print "Tail: @Tail_Stop_Locs\n";
	
	#---------------------- find the most likely 3' TIR end site ----------
	$Tail_Stop_Nums = 0;
	@Tail_Stop_Locs_New = ();
	for($i = 0; $i < @Tail_Stop_Locs - 1; $i ++) {
		$Gap_Num = 0;
		$Con_Num = 0;
		for($j = $Tail_Stop_Locs[$i]; $j > $Tail_Stop_Locs[$i] - 10 - $Gap_Num; $j --) {
			if($Con_Bps[$j] eq "-") {
				$Gap_Num ++;
			}elsif($Con_BP_Num[$j] / $Total_Num > $Min_Mul_Per) {
				$Con_Num ++;
			}else{
				next;
			}
		}
		if($Con_Num > 5) {
			$Tail_Stop_Nums ++;
			push(@Tail_Stop_Locs_New, $Tail_Stop_Locs[$i]);
		}
		last if($Tail_Stop_Nums >= $Max_Boundry);
	}
	@Tail_Stop_Locs = @Tail_Stop_Locs_New;

#	print "Tail: @Tail_Stop_Locs\n";

	#------------------------------- find TSDs for each pair of TIRs -------
	$Max_TSD_Num  = 0;
	$Max_TSD_Len  = 0;
	$Pair_Num     = 0;
	$Multi_Result = 0;
	
	foreach(@Head_Start_Locs) {
		$Start = $_;
		foreach(@Tail_Stop_Locs) {
			$Stop = $_;

			#-------------------- find TIRs ----------------------
			$Gap = 0;
			$TIR_L = "";
			for($j = 0; $j < $TIR_Len + $Gap; $j ++) {
				$BP = $Con_Bps[$Start + $j];
				if($BP eq "-") {
					$Gap ++;
				}else{
					if($j == 0) {
						$TIR_L = $BP;
					}else{
						$TIR_L = $TIR_L.$BP;
					}
				}
			}

			$Gap = 0;
			$TIR_R = "";
			for($j = 0; $j < $TIR_Len + $Gap; $j ++) {
				$BP = $Con_Bps[$Stop - $j];
				if($BP eq "-") {
					$Gap ++;
				}else{
					if($j == 0) {
						$TIR_R = $BP;
					}else{
#						print "BP $BP TIR_R $TIR_R\n";
						$TIR_R = $BP.$TIR_R;
					}
				}
			}

			$TIR_R_Rev = DNA_reverser($TIR_R);

#			print "$Start\t$Stop\t$TIR_L\t$TIR_R\n";
			next if(seq_matcher($TIR_L, $TIR_R_Rev, $TIR_Mismatch) == 0);		# doesn't considerate in/del
#			print "matched\n";

			
			#------------------------ find TSDs --------------
			for($i = 0; $i <= $Max_TSD; $i ++) {
				$TSD_Num[$i] = 0;
			}
			foreach(keys(%ID_Seq)) {
#				$Seq_Name = $_;
				$Seq  = $ID_Seq{$_};
		
				$TSD_L_Add = 0;
				$TSD_R_Add = 0;
				$TSD_L_Len = 0;
				$TSD_R_Len = 0;
				for($i = 2; $i <= $Max_TSD; $i ++) {
					while($TSD_L_Len < $i) {
						$TSD_L = substr($Seq, $Start - $i - $TSD_L_Add, $i + $TSD_L_Add);
						$TSD_L =~ s/-//g;
						$TSD_L_Len = length($TSD_L);
						if($TSD_L_Len < $i) {
							$TSD_L_Add ++;
						}
						last if ($TSD_L_Add > $Head_Edge);
					}

					if($Fixed_TSDs ne "none") {
						$Qualified = 1;
						foreach(@Fixed_Lens) {
							if($i eq $_) {
								$Qualified = 0;
								foreach(@Fixed_Seqs) {
									$Qualified = 1 if($TSD_L eq $_);
								}
							}
						}
						next if($Qualified == 0);						
					}
					
					while($TSD_R_Len < $i) {
						$TSD_R = substr($Seq, $Stop + 1, $i + $TSD_R_Add);
						$TSD_R =~ s/-//g;
						$TSD_R_Len = length($TSD_R);
						if($TSD_R_Len < $i) {
							$TSD_R_Add ++;
						}
						last if ($TSD_R_Add > $Head_Edge);
					}						

					if($TSD_L eq $TSD_R) {
#						print "$TSD_L\n";
						$TSD_Num[$i] ++;
#						print "$Seq_Name\t$Start\t$Stop\t$TIR_L\t$TIR_R\t$TSD_L\n";
					}
				}
			}

			for($i = 2; $i <= $Max_TSD; $i ++) {
				if(($i == 2)||($i == 3)) {
					$Long_Bonus = 0;
				}elsif(($i == 4)||($i == 5)) {
					$Long_Bonus = 0.1;
				}elsif(($i == 6)||($i == 7)) {
					$Long_Bonus = 0.2;
				}else{	
					$Long_Bonus = 0.3;
				}
				if($TSD_Num[$i] / $Total_Num >= $Min_TSD_Per - $Long_Bonus) {
					print "TSD length: $i\tStart: $Start\tStop: $Stop\tScore: $TSD_Num[$i]/$Total_Num > $Min_TSD_Per - $Long_Bonus\n";
					if($Max_TSD_Num < $TSD_Num[$i]) {
						if($Max_TSD_Num == 0) {
							$Max_TSD_Num = $TSD_Num[$i];
							$Max_TSD_Len = $i;
							$TE_Start = $Start;
							$TE_Stop  = $Stop;
							$Pair_Num = 1;
						}elsif($Max_TSD_Len < $i) {
							$Max_TSD_Num = $TSD_Num[$i];
							$Max_TSD_Len = $i;
							$TE_Start = $Start;
							$TE_Stop  = $Stop;
							$Pair_Num = 1;
						}elsif(($Max_TSD_Len == $i)&&($Max_TSD_Num < $TSD_Num[$i])) {	
							$Max_TSD_Num = $TSD_Num[$i];
							$Max_TSD_Len = $i;
							$TE_Start = $Start;
							$TE_Stop  = $Stop;
							$Pair_Num = 1;
						}elsif(($Max_TSD_Len == $i)&&($Max_TSD_Num == $TSD_Num[$i])){
							$Multi_Result = 1;
						}else{
							next;
						}
					}
				}
			}
		}
	}

#	print "return\n";

	$Con_Seq = join("", @Con_Bps);
	if(($Multi_Result == 1)||($Pair_Num == 0)) {
		$Con_Seq = substr($Con_Seq, $Head_Edge, $Total_Len - $Head_Edge - $Tail_Edge + 1);
		$Con_Seq =~ s/-//g;
		return(0);
	}else{
		$Con_Seq = substr($Con_Seq, $TE_Start, $TE_Stop - $TE_Start + 1);
		print "Con seq : $TE_Start - $TE_Stop\n";
		$Con_Seq =~ s/-//g;
		return(1);
	}
}

#-----------------------------------------------------
sub DNA_reverser {
    my($Seq) = @_;
	$Seq = reverse $Seq;
	$Seq =~ tr/ACGTacgt/TGCAtgca/;
    return($Seq);
}

#------------------------------------------------------------------
sub seq_matcher {
	my($Seq1, $Seq2, $Mis_Match) = @_;
	my(@BP1, @BP2, $Len, $Score, $i);
	@BP1 = split(//, $Seq1);
	@BP2 = split(//, $Seq2);
	$Len = length($Seq1);
	$Score = 0;
	for($i = 0; $i < $Len; $i ++) {
		if($BP1[$i] eq $BP2[$i]) {
			$Score ++;
		}
	}
	if($Len - $Score <= $Mis_Match) {
		return(1);
	}else{
		return(0);
	}
}
