#!/usr/bin/perl -w
# 8/19/2010
# connect step 2 and 3
# -F T
# 8/17/2010
# add max low complex seq length
# 
# 8/13/2010
# add sample function: $Sample_Per  = defined $opt_P ? $opt_P : 0.25;
# add "$Truncated   = defined $opt_C ? $opt_C : 0;"
# 10/28/2009 blast: -G 4 -E 2 -q -3 -r 2
# get putative MITE-like sequences from given genomic sequences
# -----------------------------------------------------
use Getopt::Std;
use FindBin; #shujun
# -----------------------------------------------------
my $script_path = $FindBin::Bin; #shujun

getopts("i:F:w:s:g:n:c:d:f:t:M:l:p:L:I:m:T:C:P:A:S:h:");

$Genomic_Seq = defined $opt_i ? $opt_i : "";
$Formated    = defined $opt_F ? $opt_F : 0;

$Window_Len  = defined $opt_w ? $opt_w : 2000;
$Step_Len    = defined $opt_s ? $opt_s : 1500;
$Genome_Tag  = defined $opt_g ? $opt_g : "genome";

$Group_Num   = defined $opt_n ? $opt_n : 5;
$CPU_Num     = defined $opt_c ? $opt_c : 5;

$Min_Low_Per = defined $opt_d ? $opt_d : 0.2;

$Flank_Len   = defined $opt_f ? $opt_f : 60;
$TIR_Len     = defined $opt_t ? $opt_t : 10;
$Max_TSD     = defined $opt_M ? $opt_M : 10;
$Loose_Num   = defined $opt_l ? $opt_l : 1;
$Dust_Per    = defined $opt_p ? $opt_p : 0.2;

$Exe_Min_Len = defined $opt_L ? $opt_L : 90;	# for generating exemplars. minimal shared len(%)
$Exe_Min_Ide = defined $opt_I ? $opt_I : 85;	# for generating exemplars. minimal shared identity(%)

$Min_Copy    = defined $opt_m ? $opt_m : 3;		# minimal copy number detected by PHI
$Fixed_TSD   = defined $opt_T ? $opt_T : "TA_";

$Truncated   = defined $opt_C ? $opt_C : 0;		# a parameter for MITE_scouter_part3.pl

$Sample_Per  = defined $opt_P ? $opt_P : 1;	# increase this number (at most 1) if you are pacient enough
$Max_LowCom  = defined $opt_A ? $opt_A : 90;

$Step_code   = defined $opt_S ? $opt_S : 1;		# There are 8 steps. Input eg: 12 or 345678(if you finished step 1 and 2)
$Help        = defined $opt_h ? $opt_h : "";

#-----------------------------------------------------
if($Step_code =~ /1/) {
	usuage() if((!$Genomic_Seq)||($Help));
}else{
	usuage() if($Help);
}

# -----------------------------------------------------
$Output1 = $Genome_Tag."_Step1.fa";

$Output2 = $Genome_Tag."_Step2.fa";

$Output3 = $Genome_Tag."_Step3.fa";

$Output4 = $Genome_Tag."_Step4.fa";

$Output5 = $Genome_Tag."_Step5.fa";

$Output6 = $Genome_Tag."_Step6.fa";

$Output7 = $Genome_Tag."_Step7.fa";

$Output8 = $Genome_Tag."_Step8";

$Format_Index = "$script_path/reads_indexer.pl"; #shujun

#--------------------- split genome sequnces --------------
if($Step_code =~ /1/) {

	#------------------- format database if necessary --------
	if($Formated != 1) {
		print "formating database ...\n";
		system("perl $script_path/blast_formatdb_index.pl -d $Genomic_Seq -i $Format_Index \n"); #shujun
	}

	#------- cut genomic sequeneces into fragments --------
	system("perl $script_path/fasta_windows_maker.pl -i $Genomic_Seq -w $Window_Len -s $Step_Len -P $Sample_Per -o $Output1\n"); #shujun
}

#-------------------- find raw candidates -------------------
if($Step_code =~ /2/) {
	#------- split sequence fragments into groups ---------
	system("perl $script_path/fasta_spliter.pl -i $Output1 -n $Group_Num -o $Genome_Tag\n"); #shujun
	#------- identify putative MITEs from each group ------
	@Groups = glob("$Genome_Tag.*");
	$Sub_Group_Num = 1;
	foreach(@Groups) {
		$Group_Fasta_File = $_;
		print "Group file: $Group_Fasta_File\n";
		print "-i $Group_Fasta_File -f $Flank_Len -t $TIR_Len -m $Max_TSD -l $Loose_Num -p $Min_Low_Per -g $Genome_Tag \n";
		$Sub_Tag = $Genome_Tag."_".$Sub_Group_Num;
		system("perl $script_path/MITE_Hunter_worker1.pl -i $Group_Fasta_File -f $Flank_Len -t $TIR_Len -m $Max_TSD -p $Min_Low_Per -g $Sub_Tag &\n"); #shujun
		$Sub_Group_Num ++;
	}

	$Sleep = 1;
	while($Sleep == 1) {
		$Sleep = 0;
		foreach(@Groups) {
			if(!(-e "$_.done")) {
				$Sleep = 1;
			}
		}
		sleep(60) if($Sleep == 1);
	}
}

#-------------- filter out sequences mainly composed with low complexty ----------
if($Step_code =~ /3/) {
	system("cat *.TSD.* > $Output2\n");
	system("mdust $Output2 > $Output2.dusted\n");
	system("perl $script_path/low_complexity_filter.pl -i $Output2 -d $Output2.dusted -p $Dust_Per -l $Max_LowCom -o $Output3\n"); #shujun
}

#-------------------------- BLAST to filter candidates ------------------- 
if($Step_code =~ /4/) {
	#------------------------- classify candidates into groups based on their length -----------------------
	@Len_Files = ();
	@Len_Seq_Files = ();
	open(TF, "$Output3")||die"$!\n";
	while(<TF>) {
		chomp;
		if(/^>(\S+)/) {
			$Name = $1;
		}else{
			$Seq = $_;
			$Len = int((length($Seq) - $Flank_Len * 2)/100);
			$Temp_Output = $Genome_Tag."_raw_NC_".($Len+1)*100;
			if(-e "$Temp_Output") {
				print ($Len ">$Name\n$Seq\n");
			}else{
				open($Len, ">$Temp_Output")||die"$!\n";
				push(@Len_Files, $Len);
				push(@Len_Seq_Files, $Temp_Output);
				print ($Len ">$Name\n$Seq\n");
			}
		}
	}
	close(TF);
	
	foreach(@Len_Files) {
		close($_);
	}

	#------------------------------ self blast followed by MITE_scouter_part2.pl to find putative MITEs ----------------------------------
	foreach(@Len_Seq_Files) {
		$File = $_;
		system("formatdb -i $File -o F -p F\n");
		system("blastall -i $File -d $File -e 1e-10 -p blastn -o $File.self -m 8 -v 60 -b 60 -a $CPU_Num -F F -G 4 -E 2 -q -3 -r 2\n");
		system("perl $script_path/MITE_Hunter_worker2.pl -q $File -b $File.self -c 2 -o $File.pre\n"); #shujun
	}
	system("cat *.pre.fa > $Output4");
}

#-------------------------- make exemplars ---------------------
if($Step_code =~ /5/) {
	#------------------------------ self blast followed by examplar_maker.pl to group putative MITEs ---------------------------------
	$Pre_Num = 0;
	open(OF4, "$Output4")||die"$!\n";
	while(<OF4>) {
		if(/^>/) {
			$Pre_Num ++;
		}
	}
	close(OF4);

	if($Pre_Num > 30000) {
		#------------------------- classify candidates into groups based on their length -----------------------
		@Len_Files = ();
		@Len_Seq_Files = ();
		open(TF, "$Output4")||die"$!\n";
		while(<TF>) {
			chomp;
			if(/^>(\S+)/) {
				$Name = $1;
			}else{
				$Seq = $_;
				$Len = int((length($Seq) - 120)/100);
				$Temp_Output = "temp_".$Genome_Tag."_".($Len+1)*100;
				if(-e "$Temp_Output") {
					print ($Len ">$Name\n$Seq\n");
				}else{
					open($Len, ">$Temp_Output")||die"$!\n";
					push(@Len_Files, $Len);
					push(@Len_Seq_Files, $Temp_Output);
					print ($Len ">$Name\n$Seq\n");
				}
			}
		}
		close(TF);
	
		foreach(@Len_Files) {
			close($_);
		}

		#------------------------------ self blast ----------------------------------
		foreach(@Len_Seq_Files) {
			$File = $_;
			system("formatdb -i $File -o F -p F\n");
			system("blastall -i $File -d $File -e 1e-10 -p blastn -o $File.subself -m 8 -v 2000 -b 2000 -G 4 -E 2 -q -3 -r 2 -a $CPU_Num -F F\n");
		}
		system("cat *.subself > $Output4.self\n");	
		system("rm temp*\n");
	}else{
		system("formatdb -i $Output4 -o F -p F\n");
		system("blastall -i $Output4 -d $Output4 -e 1e-10 -p blastn -o $Output4.self -m 8 -v 2000 -b 2000 -G 4 -E 2 -q -3 -r 2 -a $CPU_Num -F F\n");
	}
	system("perl $script_path/examplar_maker.pl -i $Output4 -b $Output4.self -p 0 -n 2 -o $Output5 -m $Exe_Min_Len -s $Exe_Min_Ide\n"); #shujun
}

#-------------------------- check exemplars --------------------
if($Step_code =~ /6/) {
	$Max_Gap = $Max_LowCom + 10;
	#------------------------------ check whether the putative MITEs are real, get consensus sequences and group ----------------------
	system("perl $script_path/MITE_Hunter_worker3.pl -i $Output5 -D $Genomic_Seq -m $Min_Copy -I $Genome_Tag -f $Fixed_TSD -c $Truncated -o $Output6 -G $Max_Gap\n"); #shujun
}

#-------------------- further group and find exemplars again -------------------
if($Step_code =~ /7/) {
	#------------------------------ self blast followed by examplar_maker.pl to group putative MITEs ---------------------------------
	system("formatdb -i $Output6 -o F -p F\n");
	system("blastall -i $Output6 -d $Output6 -e 1e-10 -p blastn -o $Output6.self -m 8 -G 4 -E 2 -q -3 -r 2 -a $CPU_Num -F F\n");
	system("perl $script_path/examplar_maker.pl -i $Output6 -b $Output6.self -p 0 -n 2 -o $Output7 -m $Exe_Min_Len -s $Exe_Min_Ide -l 1\n"); #shujun
}

#-------------------- group into subfamilies -------------------
if($Step_code =~ /8/) {
	system("formatdb -i $Output7 -o F -p F\n");
	system("blastall -i $Output7 -d $Output7 -e 1e-10 -p blastn -o $Output7.self -v 500 -b 500 -m 9 -G 4 -E 2 -q -3 -r 2 -a $CPU_Num\n");
	system("perl $script_path/MITE_Hunter_worker4.pl -i $Output7 -b $Output7.self -o $Output8\n"); #shujun
}

# -----------------------------------------------------
sub usuage {
    print "\n","Hi,", ' need some help? @_@',"\n";   
    print STDERR <<"    _EOT_";

    Usage: program <options> <specification file> <default>
    
	\$Genomic_Seq = defined \$opt_i ? \$opt_i : "";
	\$Formated    = defined \$opt_F ? \$opt_F : 0;

	\$Window_Len  = defined \$opt_w ? \$opt_w : 2000;
	\$Step_Len    = defined \$opt_s ? \$opt_s : 1500;
	\$Genome_Tag  = defined \$opt_g ? \$opt_g : "genome";

	\$Group_Num   = defined \$opt_n ? \$opt_n : 5;

	\$Min_Low_Per = defined \$opt_d ? \$opt_d : 0.2; # sequnces with more than $Min_Low_Per percent low complexity will be filtered

	\$Flank_Len   = defined \$opt_f ? \$opt_f : 60;
	\$TIR_Len     = defined \$opt_t ? \$opt_t : 10;
	\$Max_TSD     = defined \$opt_M ? \$opt_M : 10;
	\$Loose_Num   = defined \$opt_l ? \$opt_l : 1;
	\$Dust_Per    = defined \$opt_p ? \$opt_p : 0.2;

	\$Exe_Min_Len = defined \$opt_L ? \$opt_L : 90;	# for generating exemplars. minimal shared len(%)
	\$Exe_Min_Ide = defined \$opt_I ? \$opt_I : 85;	# for generating exemplars. minimal shared identity(%)
	
	\$Min_Copy    = defined \$opt_m ? \$opt_m : 3;	# minimal copy number detected by PHI
	\$Fixed_TSD   = defined \$opt_T ? \$opt_T : "TA_"; # seperated by _ like TA_TAA 
	
	\$Truncated   = defined \$opt_C ? \$opt_C : 0;
	\$Sample_Per  = defined \$opt_P ? \$opt_P : 1;	# increase this number (at most 1) if you are pacient enough
	\$Max_LowCom  = defined \$opt_A ? \$opt_A : 90;

	\$Step        = defined \$opt_S ? \$opt_S : 1;	# There are 8 steps. Input eg: 12 or 345678(if you finished step 1 and 2)
	\$Help        = defined \$opt_h ? \$opt_h : "";
           
    _EOT_
    exit(1);
}
