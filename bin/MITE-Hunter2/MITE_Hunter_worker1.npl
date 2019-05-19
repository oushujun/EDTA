#!/usr/bin/perl -w
# 8/19/2010
# output a tag file to wake sleep
# 8/13/2010
# 7/25/2010
# cancel TA only for 2 bp TSD 
# get putative MITE-like sequences from given genomic sequences
# -----------------------------------------------------
use Getopt::Std;
# -----------------------------------------------------
getopts("i:f:t:m:p:g:T:h:");

$Input       = defined $opt_i ? $opt_i : "";
$Flank_Len   = defined $opt_f ? $opt_f : 60;
$TIR_Len     = defined $opt_t ? $opt_t : 10;
$Max_TSD     = defined $opt_m ? $opt_m : 10;
$Min_Low_Per = defined $opt_p ? $opt_p : 0.2;
$ID_Tag      = defined $opt_g ? $opt_g : "scouter";
$Fixed_TSDs  = defined $opt_T ? $opt_T : "TA_";
$Help        = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));

# -----------------------------------------------------
if($Fixed_TSDs ne "none") {
	@Fixed_Seqs = split(/_/, $Fixed_TSDs);
	foreach(@Fixed_Seqs) {
		push(@Fixed_Lens, length($_));
	}
}

#-----------------------------------------
print "Reading data...";
open(IF, "$Input")||die "$!";
while(<IF>) {
	chomp;
	if(/^>(\S+)/) {
		$Name = $1;
	}else{
		$Name_Seq{$Name} .= $_;
	}
}
close(IF);

print "Searching for MITEs ...\n";

$MITE_Num = 1;
%MITEs_Identified = ();

$Num = 0;
foreach(keys(%Name_Seq)) {
	if($Num % 1000 == 0) {
		print "$Num\n";
	}
	$Num ++;

	$Name = $_;
	$Seq  = $Name_Seq{$Name};
	$Len  = length($Seq);

	for($i = $Flank_Len; $i < $Len - $Flank_Len; $i ++) {
		$TIR = substr($Seq, $i, $TIR_Len);

		next if($TIR =~ /N/);
	
		next if(low_complex_checker($TIR) == 1);     # filter out low complexy sequences 

		# ----------------------------------------------------	get all putative TIRs
		$TIR_Rev = DNA_reverser($TIR);

		@Loose_TIR_Rev = loose_maker($TIR_Rev);

		foreach(@Loose_TIR_Rev) {
			$TIR_Rev = $_;
#			next if(low_complex_checker($TIR_Rev) == 1);     # filter out low complexy sequences 
	
			$Found = $i + 40;	# minimal length of candidate

			$Found = index($Seq, $TIR_Rev, $Found + 1);

			if($Found != -1) {
				($TSD, $TSD_Len) = TSD_Finder($Seq, $i, $Found + $TIR_Len, $Max_TSD);
				if($TSD ne "none") {
#					print "$TIR\t$TIR_Rev\n";

					$Seq_With_Flank = substr($Seq, $i - $Flank_Len, $Found + $TIR_Len + 2*$Flank_Len - $i);

					next if($Seq_With_Flank =~ /N/);

					if(defined($MITEs_Identified{$Seq_With_Flank." ".$TSD_Len})) {
						next;
					}else{
						$MITEs_Identified{$Seq_With_Flank." ".$TSD_Len} = 1;
					}

					$Putative_MITE_Name = $ID_Tag."_".$MITE_Num;
					if(defined($TSD_Len_Found{$TSD_Len})) {
						print ($TSD_Len ">$Putative_MITE_Name Contig:$Name TIR:$TIR $TIR_Rev TSD:$TSD\n$Seq_With_Flank\n");
						$MITE_Num ++;						
					}else{
						$TSD_Len_Found{$TSD_Len} = 1;
						open($TSD_Len, ">$Input.TSD.$TSD_Len")||die "$!";
						print ($TSD_Len ">$Putative_MITE_Name Contig:$Name TIR:$TIR $TIR_Rev TSD:$TSD\n$Seq_With_Flank\n");
						$MITE_Num ++;
					}
					$i += 5;	# jump to find next possiable MITE (even if it miss the correct ones, but it is still within 5 bps, which can be still detected)
				}
			}
		}
	}
}

foreach(keys(%TSD_Len_Found)) {
	close($_);
}

open(TF, ">$Input.done")||die"$!\n";
print (TF "Wake up!\n");
close(TF);

# -----------------------------------------------------
sub usuage {
    print "\n","Hi,", ' need some help? @_@',"\n";   
    print STDERR <<"    _EOT_";

    Usage: program <options> <specification file> <default>
    
	\$Input       = defined \$opt_i ? \$opt_i : "";
	\$Flank_Len   = defined \$opt_f ? \$opt_f : 60;
	\$TIR_Len     = defined \$opt_t ? \$opt_t : 10;
	\$Max_TSD     = defined \$opt_m ? \$opt_m : 10;
	\$Min_Low_Per = defined \$opt_p ? \$opt_p : 0.2;
	\$ID_Tag      = defined \$opt_g ? \$opt_g : "scouter";
	\$Fixed_TSDs  = defined \$opt_T ? \$opt_T : "TA_";
	\$Help        = defined \$opt_h ? \$opt_h : "";
           
    _EOT_
    exit(1);
}

#-----------------------------------------------------
sub DNA_reverser {
    my($Seq) = @_;
	$Seq = reverse $Seq;
	$Seq =~ tr/ACGTacgt/TGCAtgca/;
    return($Seq);
}

#-----------------------------------------------------
sub loose_maker {
	my($Left, $Mid, $Right, $Len, @Loose_Seqs, $i, %Loose_Seqs, @BPs);
	my($Seq) = @_;
	@Bps = ("A", "T", "C", "G");
	$Len = length($Seq);
	@Loose_Seqs = ();
	push(@Loose_Seqs, $Seq);

	@Seq_BPs = split(//, $Seq);

	for($i = 0; $i < $Len; $i ++) {
		if($i == 0) {
			$Right = substr($Seq, 1);
			foreach(@Bps) {
				if($_ ne $Seq_BPs[$i]) {
					push(@Loose_Seqs, $_.$Right);
				}
			}
		}elsif($i == $Len - 1) {
			$Left  = substr($Seq, 0, $Len - 1);
			foreach(@Bps) {
				if($_ ne $Seq_BPs[$i]) {
					push(@Loose_Seqs, $Left.$_);
				}
			}
		}else{
			$Left  = substr($Seq, 0, $i);
			$Right = substr($Seq, $i + 1);
			foreach(@Bps) {
				if($_ ne $Seq_BPs[$i]) {
					push(@Loose_Seqs, $Left.$_.$Right);
				}
			}
		}
	}
	return(@Loose_Seqs);
}

#-----------------------------------------------------
sub low_complex_checker {
	my(@BPs, $AT, $CG, $Low);
	my($Seq) = @_;
	@BPs = split(//, $Seq);
	$AT = 0;
	$CG = 0;
	foreach(@BPs) {
		$AT ++ if(($_ eq "A")||($_ eq "T"));
		$CG ++ if(($_ eq "C")||($_ eq "G"));
	}

	$Low = 0;
	if($AT * $CG == 0) {
		$Low = 1;
	}elsif(($Seq =~ /ATATATAT/)||($Seq =~ /TATATATA/)) {
		$Low = 1;
	}elsif(($Seq =~ /AAAAAAAA/)||($Seq =~ /TTTTTTTT/)) {
		$Low = 1;
	}elsif(($Seq =~ /CGCGCGCG/)||($Seq =~ /GCGCGCGC/)) {
		$Low = 1;
	}elsif(($Seq =~ /CCCCCCCC/)||($Seq =~ /GGGGGGGG/)) {
		$Low = 1;
	}elsif($CG/($AT + $CG) <= $Min_Low_Per) {
		$Low = 1;
	}
	return($Low);
}

#----------------------------------------------------
sub TSD_Finder {
	my($Seq, $Start, $Stop, $Max_TSD) = @_;
	my($i, $TSD_Seq, $TSD_Len, $TSD_L, $TSD_R);
	$TSD_Seq = "none";
	$TSD_Len = "0";
	for($i = $Max_TSD; $i >= 2; $i --) {
		$TSD_L = substr($Seq, $Start - $i, $i);
		$TSD_R = substr($Seq, $Stop, $i);

		if($TSD_L eq $TSD_R) {
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
			$TSD_Seq = $TSD_L;
			$TSD_Len = $i;
			last;
		}
	}
	return($TSD_Seq, $TSD_Len);
}