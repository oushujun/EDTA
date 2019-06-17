#!/usr/bin/perl -w
# 11/21/2010
# fix the formatdb
# 11/12/2009 a t c g -> A T C G
# if want to check the gene annoation information from the corrosponding genbank file, set -g 1 
# get rid of the \n at the end of each line
# altix
#-----------------------------------------------------
use Getopt::Std;
use FindBin; #shujun
#-----------------------------------------------------
my $script_path = $FindBin::Bin; #shujun
getopts("d:g:i:h:");

$Input   = defined $opt_d ? $opt_d : "";
$GenBank = defined $opt_g ? $opt_g : 0;
$Indexer = defined $opt_i ? $opt_i : "$script_path/reads_indexer.pl"; #shujun
$Help    = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------
if(($GenBank != 0)&&(-e "$Input.gb")){
	open(GF, "$Input.gb")||die"$!\n";
	while(<GF>) {
		chomp;
		$Line = $_;
		if(/ACCESSION\s+(\S+)\s/) {
			$Accession = $1;
			$Accession_Nums{$Accession} = 1;
		}
	}
	close(GF);
}else{
	$GenBank = 0;
}

open(IF, "$Input")||die"$!\n";
open(OF, ">$Input.temp")||die"$!\n";
$First = 1;
while(<IF>) {
	if(/^>/) {
		if($GenBank == 0) {
			if($First == 1) {
				print(OF "$_");
				$First = 0;
			}else{
				print(OF "\n$_");
			}
		}else{
			$Line = $_;
			if($First == 1) {
				$Acc_Found = 0;
				foreach(keys(%Accession_Nums)) {
					if($Line =~ /$_/) {
						print(OF ">$_\n");						
						$Acc_Found = 1;
						last;
					}
				}
				if($Acc_Found == 0) {
					print "Warnning: the information of $_ cannot be found in the GenBank file!";
				}

				$First = 0;
			}else{
				$Acc_Found = 0;
				foreach(keys(%Accession_Nums)) {
					if($Line =~ /$_/) {
						print(OF "\n>$_\n");
						$Acc_Found = 1;
						last;						
					}
				}
				if($Acc_Found == 0) {
					print "Warnning: the information of $_ cannot be found in the GenBank file!";
				}
			}
		}
	}else{
		$Line = $_;
		$Line =~ s/\s//g;
		$Line = uc($Line);
		print(OF "$Line");
	}
}
close(IF);
close(OF);

system "mv $Input.temp $Input\n";

system "formatdb -i $Input -o F -p F\n";

system "perl $Indexer -i $Input\n";

#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :tblastn_copy_finder.pl <options> <specification> <default>

	\$Input   = defined \$opt_d ? \$opt_d : "";
	\$GenBank = defined \$opt_g ? \$opt_g : 0;
	\$Indexer = defined \$opt_i ? \$opt_i : "$script_path/reads_indexer.pl";
	\$Help    = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}
