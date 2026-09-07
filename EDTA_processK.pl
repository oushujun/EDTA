#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;

#####################################################################
##### Perform EDTA basic and advance filtering on TE candidates #####
##### Shujun Ou (shujun.ou.1@gmail.com, 12/28/2023)             #####
#####################################################################

## Input:
#	$genome.SINE.raw.fa
#	$genome.LINE.raw.fa
#	$genome.LTR.raw.fa
#	$genome.LTR.intact.raw.fa
#	$genome.TIR.intact.raw.fa
#	$genome.Helitron.intact.raw.fa

## Output:
#	$genome.EDTA.fa.stg1

my $usage = "\nPerform EDTA basic and advance filtering for raw TE candidates and generate the stage 1 library
	perl EDTA_processK.pl [options]
		-genome	[File]	The genome FASTA
		-ltr	[File]	The raw LTR library FASTA
		-ltrint	[File]	The intact LTR library FASTA
		-sine	[File]	The raw SINE library FASTA
		-line	[File]	The raw LINE library FASTA
		-tir	[File]	The raw TIR library FASTA
		-helitron	[File]	The raw Helitron library FASTA
		-mindiff_ltr	[float]	The minimum fold difference in richness between LTRs and contaminants (default: 1)
		-mindiff_tir	[float]	The minimum fold difference in richness between TIRs and contaminants (default: 1)
		-mindiff_hel	[float]	The minimum fold difference in richness between Helitrons and contaminants (default: 1.5)
		-repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
		-blast [path]	The directory containing Blastn (default: read from ENV)
		-threads|-t	[int]	Number of theads to run this script
		-help|-h	Display this help info
\n";

# user input
my $genome = '';
my $LTRraw = '';
my $LTRintact = '';
my $SINEraw = '';
my $LINEraw = '';
my $TIRraw = '';
my $HELraw = '';
my $err = '';

# minimum richness difference between $TE1 and $TE2 for a sequence to be considered as REAL to $TE1.
# Smaller number is more inclusive during purging, hence higher false positives
my $mindiff_LTR = 1;
my $mindiff_TIR = 1;
my $mindiff_HEL = 1.5;

my $threads = 4;
my $script_path = $FindBin::Bin;
my $TE_purifier = "$script_path/bin/TE_purifier.pl";
my $rename_TE = "$script_path/bin/rename_TE.pl";
my $cleanup_tandem = "$script_path/bin/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/bin/cleanup_nested.pl";
my $cleanup_proteins = "$script_path/bin/cleanup_proteins.pl";
my $repeatmasker = '';
my $blast = '';

# read parameters
my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$LTRraw = $ARGV[$k+1] if /^-ltr$/i and $ARGV[$k+1] !~ /^-/;
	$LTRintact = $ARGV[$k+1] if /^-ltrint$/i and $ARGV[$k+1] !~ /^-/;
	$SINEraw = $ARGV[$k+1] if /^-sine$/i and $ARGV[$k+1] !~ /^-/;
	$LINEraw = $ARGV[$k+1] if /^-line$/i and $ARGV[$k+1] !~ /^-/;
	$TIRraw = $ARGV[$k+1] if /^-tir/i and $ARGV[$k+1] !~ /^-/;
	$HELraw = $ARGV[$k+1] if /^-helitron/i and $ARGV[$k+1] !~ /^-/;
	$mindiff_LTR = $ARGV[$k+1] if /^-mindiff_ltr/i and $ARGV[$k+1] !~ /^-/;
	$mindiff_TIR = $ARGV[$k+1] if /^-mindiff_tir/i and $ARGV[$k+1] !~ /^-/;
	$mindiff_HEL = $ARGV[$k+1] if /^-mindiff_hel/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker/i and $ARGV[$k+1] !~ /^-/;
	$blast = $ARGV[$k+1] if /^-blast/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-help$|^-h$/i;
	$k++;
        }

# check files and dependencies
die "Genome file $genome not exists!\n$usage" unless -s $genome;
die "LTR raw library file $LTRraw not exists!\n$usage" unless -s $LTRraw;
die "Intact LTR file $LTRintact not exists!\n$usage" unless -s $LTRintact;
#die "LINE raw library file $LINEraw not exists!\n$usage" unless -e $LINE; # allow empty file
#die "SINE raw library file $SINEraw not exists!\n$usage" unless -e $SINE; # allow empty file
die "TIR raw library file $TIRraw not exists!\n$usage" unless -e $TIRraw;
print STDERR "Warning: The TIR raw library $TIRraw is empty (0 bp). Continuing with an empty TIR library.\n" unless -s $TIRraw;
die "Helitron raw library file $HELraw not exists!\n$usage" unless -e $HELraw;
print STDERR "Warning: The Helitron raw library $HELraw is empty (0 bp). Continuing with an empty Helitron library.\n" unless -s $HELraw;
die "The script TE_purifier.pl is not found in $TE_purifier!\n" unless -s $TE_purifier;
die "The script rename_TE.pl is not found in $rename_TE!\n" unless -s $rename_TE;
die "The script cleanup_tandem.pl is not found in $cleanup_tandem!\n" unless -s $cleanup_tandem;
die "The script cleanup_nested.pl is not found in $cleanup_nested!\n" unless -s $cleanup_nested;
die "The script cleanup_proteins.pl is not found in $cleanup_proteins!\n" unless -s $cleanup_proteins;

# make a softlink to the genome
my $genome_file = basename($genome);
`ln -s $genome $genome_file` unless -e $genome_file;
$genome = $genome_file;
my $LTR = "$genome.LTR.raw.fa";
my $LTRint = "$genome.LTR.intact.raw.fa";
my $SINE = "$genome.SINE.raw.fa";
my $LINE = "$genome.LINE.raw.fa";
my $TIR = "$genome.TIR.intact.raw.fa";
my $HEL = "$genome.Helitron.intact.raw.fa";

# Make working directories
`mkdir $genome.EDTA.combine` unless -e "$genome.EDTA.combine" && -d "$genome.EDTA.combine";

# enter the combine folder for EDTA processing
chdir "$genome.EDTA.combine" or die "Cannot enter $genome.EDTA.combine: $!\n";

# --- Resume support (added 2026-08-24) --------------------------------------
# EDTA_processK.pl originally had no restart logic: every invocation redid the
# whole filtering stage from scratch. On an 11 Gb genome that stage takes far
# longer than one wall-clock window, so a timeout used to throw away every
# completed step. Each block below now drops a dot-stamp on completion and is
# skipped if its stamp is present. Stamps are named ".<step>.done" -- the
# leading dot keeps them clear of the `rm $genome*` / `rm *.ndb` cleanups in
# this script and in EDTA.pl, both of which glob on non-dot names.
# Remove ./.step*.done (or pass --overwrite 1 to EDTA.pl) to force a rerun.
sub done { return (-e ".$_[0].done") ? 1 : 0; }
# A stamp is written only after the step's expected products are verified: files in
# $nonempty must have content, files in $may_empty must exist (TIR/Helitron/SINE products
# are legitimately 0 bp in genomes lacking those TEs). mark() dies without stamping otherwise.
sub mark {
	my ($step, $nonempty, $may_empty) = @_;
	$nonempty = [] unless defined $nonempty;
	$may_empty = [] unless defined $may_empty;
	my @missing = ((grep { ! -s $_ } @$nonempty), (grep { ! -e $_ } @$may_empty));
	die "Step $step did not produce expected output(s), stamp not written: @missing\n" if @missing;
	`touch ".$step.done"`;
	}
# ----------------------------------------------------------------------------
`cp ../$LTRraw $LTR`;
`cp ../$LTRintact $LTRint`;
`cp ../$SINEraw $SINE`;
`cp ../$LINEraw $LINE`;
`cp ../$TIRraw $TIR`;
`cp ../$HELraw $HEL`;
die "Failed to stage raw TE libraries into $genome.EDTA.combine (missing inputs for: "
	. join(", ", grep { !-e $_ } ($LTR, $LTRint, $SINE, $LINE, $TIR, $HEL)) . ").\n"
	. "EDTA_processK.pl must be run from the run directory, as EDTA.pl does.\n"
	unless -e $LTR and -e $LTRint and -e $SINE and -e $LINE and -e $TIR and -e $HEL;


##################################
######  define subroutines  ######
##################################

# purify $TE2 contaminants in $TE1
# This function better works for redundant libraries
sub Purifier() {
	my ($TE1, $TE2, $mindiff) = ($_[0], $_[1], $_[2]);
	# mark contaminents with lowercase letters based on relative richness
	if (-s $TE1 and -s $TE2){
		`perl $TE_purifier -TE1 $TE1 -TE2 $TE2 -t $threads -mindiff $mindiff`;
		} else {
		`cp $TE1 $TE1-$TE2.fa`; # empty $TE1 or $TE2: no purging possible, keep the file chain flowing
		}
	# remove lowercase sequences
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 1 -minrm 1 -trf 0 -f $TE1-$TE2.fa > $TE1.HQ`;
	}

# cat with existence checks: die listing missing mandatory inputs (LTR/LINE classes),
# warn and skip missing optional ones (TIR/Helitron/SINE classes may be absent or 0 bp)
sub cat_checked {
	my ($dst, $mandatory, $optional) = @_;
	my @missing = grep { ! -e $_ } @$mandatory;
	die "Missing mandatory input(s) for $dst: @missing\n" if @missing;
	my @inputs = @$mandatory;
	foreach my $file (@$optional){
		if (-e $file){ push @inputs, $file; }
		else { print STDERR "Warning: optional input $file not found, skipped in $dst\n"; }
		}
	`cat @inputs > $dst.tmp.$$ && mv $dst.tmp.$$ $dst` if @inputs;
	die "Failed to generate $dst: cat/mv exited with status $?\n" if @inputs and $? != 0;
	`touch $dst` unless @inputs;
	}


#################################
###### Advance filtering ######
#################################

## Purge contaminants in redundant libraries
# purify raw LTR (clean LTR library is better than dirty intact LTR for purging LTRs from other TEs)
unless (&done("step1_LTR_vs_TIR")){
	&Purifier("$LTR", "$TIR", $mindiff_LTR);
	&mark("step1_LTR_vs_TIR", ["$LTR.HQ"]);
	}
unless (&done("step2_LTR_vs_HEL")){
	&Purifier("$LTR.HQ", "$HEL", $mindiff_LTR);
	`mv $LTR.HQ.HQ $LTR.HQ`;
	&mark("step2_LTR_vs_HEL", ["$LTR.HQ"]);
	}

# purify Helitron
unless (&done("step3_HEL_vs_TIR")){
	&Purifier("$HEL", "$TIR", $mindiff_HEL);
	&mark("step3_HEL_vs_TIR", [], ["$HEL.HQ"]);
	}
unless (&done("step4_HEL_vs_LTR")){
	&Purifier("$HEL.HQ", "$LTR", $mindiff_LTR);
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 0 -minrm 1 -trf 0 -f $HEL.HQ-$LTR.fa > $HEL.int.cln`; # more relaxed in filtering intact helitrons
	`mv $HEL.HQ.HQ $HEL.cln`;
	&mark("step4_HEL_vs_LTR", [], ["$HEL.int.cln", "$HEL.cln"]);
	}

# purify TIR
unless (&done("step5_TIR_vs_LTR")){
	&Purifier("$TIR", "$LTR", $mindiff_TIR);
	&mark("step5_TIR_vs_LTR", [], ["$TIR.HQ"]);
	}
unless (&done("step6_TIR_vs_HEL")){
	&Purifier("$TIR.HQ", "$HEL", $mindiff_TIR);
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 0 -minrm 1 -trf 0 -f $TIR.HQ-$HEL.fa > $TIR.int.cln`; # more relaxed in filtering intact TIRs
	`mv $TIR.HQ.HQ $TIR.cln`;
	&mark("step6_TIR_vs_HEL", [], ["$TIR.int.cln", "$TIR.cln"]);
	}

# purify intact LTR from TIRs. Including Helitron is too damaging for now.
unless (&done("step7_LTRint_vs_TIRcln")){
	&Purifier("$LTRint", "$TIR.cln", 10); # 10 is permissive
	`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 0 -minrm 1 -trf 0 -f $LTRint-$TIR.cln.fa > $LTRint.cln`;
	&mark("step7_LTRint_vs_TIRcln", ["$LTRint.cln"]);
	}
#&Purifier("$LTRint.HQ", "$HEL.cln", 10); # 10 is permissive
#`perl $cleanup_tandem -misschar l -Nscreen 1 -nc 50000 -nr 0.8 -minlen 80 -cleanN 1 -cleanT 0 -minrm 1 -trf 0 -f $LTRint.HQ-$HEL.cln.fa > $LTRint.cln`; # more relaxed in filtering intact LTRs

## Purge contaminants in non-redundant libraries
# clean LINEs in LTRs
if (&done("step8_LINE_in_LTR")){
	# skip: $LTR.cln already produced
	} elsif (-s "$LINE"){
	$err = `${repeatmasker}RepeatMasker -e ncbi -pa $threads -q -no_is -nolow -div 40 -lib $LINE $LTR 2>&1`;
	die "RepeatMasker failed in step8_LINE_in_LTR (exit status $?):\n$err\n" if $? != 0 and $err !~ /No repetitive sequences were detected/i;
	if ($err !~ /done/) {
        	`rm -f $LTR.masked; cp $LTR $LTR.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: No sequences were masked/si;
	        print STDERR "\n$err\n";
        	}
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $LTR.masked > $LTR.cln`;
	&mark("step8_LINE_in_LTR", ["$LTR.cln"]);
	} else {
		`cp $LTR $LTR.cln`;
		&mark("step8_LINE_in_LTR", ["$LTR.cln"]);
	}

# clean LINEs and LTRs in SINEs
if (&done("step9_LINE_LTR_in_SINE")){
	# skip: $SINE.cln already produced
	} elsif (-s "$SINE"){
	&cat_checked("$genome.LINE_LTR.raw.fa", ["$LTR.cln", "$LINE"], []);
	$err = `${repeatmasker}RepeatMasker -e ncbi -pa $threads -q -no_is -nolow -div 40 -lib $genome.LINE_LTR.raw.fa $SINE 2>&1`;
	die "RepeatMasker failed in step9_LINE_LTR_in_SINE (exit status $?):\n$err\n" if $? != 0 and $err !~ /No repetitive sequences were detected/i;
	if ($err !~ /done/) {
        	`rm -f $SINE.masked; cp $SINE $SINE.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: No sequences were masked/si;
	        print STDERR "\n$err\n";
        	}
	`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -f $SINE.masked > $SINE.cln`;
	&mark("step9_LINE_LTR_in_SINE", [], ["$SINE.cln"]);
	} else {
		-e $SINE ? `cp $SINE $SINE.cln` : `touch $SINE.cln`; # empty or absent SINE library contributes nothing
		&mark("step9_LINE_LTR_in_SINE", [], ["$SINE.cln"]);
	}


## clean LTRs and nonLTRs in TIRs and Helitrons
unless (&done("step10_mask_TIR_HEL")){
	my @tirhel;
	foreach my $file ("$TIR.cln", "$HEL.cln"){
		if (-e $file){ push @tirhel, $file; }
		else { print STDERR "Warning: optional input $file not found, skipped in $genome.TIR.Helitron.fa.stg1.raw\n"; }
		}
	if (@tirhel){
		`cat @tirhel | perl -nle 's/>/\\n>/g unless /^>/; print \$_' > $genome.TIR.Helitron.fa.stg1.raw`;
		} else {
		`touch $genome.TIR.Helitron.fa.stg1.raw`; # empty TIR+Helitron library is a legitimate outcome
		}
	&cat_checked("$genome.LTR.SINE.LINE.fa", ["$LTR.HQ", "$LINE"], ["$SINE.cln"]);
	if (-s "$genome.TIR.Helitron.fa.stg1.raw"){
		$err = `${repeatmasker}RepeatMasker -e ncbi -pa $threads -q -no_is -nolow -div 40 -lib $genome.LTR.SINE.LINE.fa $genome.TIR.Helitron.fa.stg1.raw 2>&1`;
		die "RepeatMasker failed in step10_mask_TIR_HEL (exit status $?):\n$err\n" if $? != 0 and $err !~ /No repetitive sequences were detected/i;
		if ($err !~ /done/) {
			`rm -f $genome.TIR.Helitron.fa.stg1.raw.masked; cp $genome.TIR.Helitron.fa.stg1.raw $genome.TIR.Helitron.fa.stg1.raw.masked` if $err =~ s/^.*(No repetitive sequences were detected.*)\s+$/Warning: No sequences were masked/si;
			print STDERR "\n$err\n";
			}
		`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.Helitron.fa.stg1.raw.masked > $genome.TIR.Helitron.fa.stg1.raw.cln`;
		} else {
		`cp $genome.TIR.Helitron.fa.stg1.raw $genome.TIR.Helitron.fa.stg1.raw.cln`; # RepeatMasker cannot take an empty query
		}
	&mark("step10_mask_TIR_HEL", ["$genome.LTR.SINE.LINE.fa"], ["$genome.TIR.Helitron.fa.stg1.raw.cln"]);
	}


## cluster TIRs and Helitrons and make stg1 raw library
unless (&done("step11_cleanup_nested")){
	if (-s "$genome.TIR.Helitron.fa.stg1.raw.cln"){
		`perl $cleanup_nested -in $genome.TIR.Helitron.fa.stg1.raw.cln -threads $threads -minlen 80 -cov 0.95 -blastplus $blast`;
		} else {
		`cp $genome.TIR.Helitron.fa.stg1.raw.cln $genome.TIR.Helitron.fa.stg1.raw.cln.cln`; # cleanup_nested dies on empty input
		}
	&mark("step11_cleanup_nested", [], ["$genome.TIR.Helitron.fa.stg1.raw.cln.cln"]);
	}
&cat_checked("$genome.EDTA.fa.stg1", ["$LTR.cln", "$LINE"], ["$SINE.cln", "$genome.TIR.Helitron.fa.stg1.raw.cln.cln"]);

## generate clean intact TEs
&cat_checked("$genome.EDTA.intact.fa.cln", ["$LTRint.cln", "$LINE"], ["$SINE.cln", "$TIR.int.cln", "$HEL.int.cln"]);

## clean up the folder
`rm *.ndb *.not *.ntf *.nto *.cat.gz *.cat *.masked *.ori.out *.nhr *.nin *.nsq *.njs 2>/dev/null`;

chdir '..';
