Generic Repeat Finder (GRF) is a C++ program package for detecting terminal inverted repeats (TIRs), terminal direct repeats (TDRs), interspersed repeats, miniature inverted repeat transposable elements (MITEs), and long terminal repeat (LTR) transposons in genomes.

If you have any question or comment, please contact Dr. Chun Liang (liangc@miamioh.edu).

This software is an open-source tool that follows specifications from the website (http://creativecommons.org/licenses/by-nc-sa/3.0/)

###############################################################################

Content

[1] Installation of GRF
[2] Installation of third-party software
[3] How to use GRF
[4] How to filter out tandem repeats in outputs using phobos
[5] Modified LTR_FINDER
[6] Test data

###############################################################################

[1] Installation of GRF

Download and unzip our source code and executable programs - "grf.XXX.tar.gz".

# tar zxf grf.XXX.tar.gz

The source code is under "grf.XXX/src/".

If you want to compile the source code yourself, 
(1) go to "grf.XXX/src/", and run "make" (easy way)
 
or (2) use g++ to compile with flag: "-std=c++11 -fopenmp". e.g.,
# cd grf.XXX/src/grf-main
# g++ main.cpp DetectMITE.h DetectMITE.cpp DetectIR.h DetectIR.cpp functions.h functions.cpp -std=c++11 -fopenmp -O3 -o grf-main

The executable programs are under "grf.XXX/bin/". Optionally, you can add the executables to your $PATH for convenience. e.g., 		

# echo 'export PATH=$PATH:<your_path>/grf.XXX/bin/' >> ~/.bashrc

###############################################################################

[2] Installation of third-party software

[2.1] CD-HIT

CD-HIT (https://code.google.com/archive/p/cdhit/downloads) is required for MITE detection. We tested version 4.6.1.

# wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/cdhit/cd-hit-v4.6.1-2012-08-27.tgz
# tar zxf cd-hit-v4.6.1-2012-08-27.tgz
# cd cd-hit-v4.6.1-2012-08-27
# make

Optionally, you can add the executables to your $PATH for convenience. e.g., 		

# echo 'export PATH=$PATH:<your_path>/cd-hit-v4.6.1-2012-08-27' >> ~/.bashrc

[2.2] g++ and gcc (version 4.9+)

Below is the commands to install version 4.9:

# sudo add-apt-repository ppa:ubuntu-toolchain-r/test
# sudo apt-get update; sudo apt-get install gcc-4.9 g++-4.9
# sudo update-alternatives --remove-all gcc
# sudo update-alternatives --remove-all g++
# sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 20
# sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 20
# sudo update-alternatives --config gcc
# sudo update-alternatives --config g++
# gcc -v
# g++ -v

###############################################################################

[3] How to use GRF

Note: All the commands below assume that the path to "grf.XXX/bin/" and "cd-hit-v4.6.1-2012-08-27/" have been added to your $PATH.

[optional] OpenMP schedule environment variable can be set at runtime by command:

# export OMP_SCHEDULE=[strategy]

Details can be found at https://software.intel.com/en-us/articles/openmp-loop-scheduling

The program "grf-main" performs genome-wide TIR, TDR, and MITE candidate detection for input genome sequences (see [3.1]-[3.3]).

Usage: grf-main [options]
[mandatory parameters]
-i <string>  Input genome sequence file in FASTA format.
-c <int>  Choice of analysis: 0: TIR detection; 1: MITE candidate detection; 2: TDR detection.
-o <string>  Output directory.

[optional parameters]
-t <int>  Number of threads used in this program; default = 1.
-f <int> Format for outputs; 0: FASTA format; 1: only IDs; default = 0.
--max_mismatch <int>  Maximum number of mismatches allowed in the terminal repeats; set -1 to be unlimited; default = -1.
--max_indel <int>  Maximum number of indels allowed in the terminal repeats; set -1 to be unlimited; default = -1.
-p <int>  Maximum percentage of unpaired nucleotides in the terminal repeats; set -1 to be unlimited; default = 10.
-r <float>  Maximum length ratio of spacer/total sequence; set -1 to be unlimited; default = -1.
-s <int>  Length of the seed region; default = 10; must >= 5 and <= '--min_tr'.
--seed_mismatch <int>  Maximum mismatch number in the seed region; default = 1.
--min_space <int>  Minimum distance between two seed regions; for TIRs/TDRs, default = 0; for MITEs, default = 30.
--max_space <int>  Maximum distance between two seed regions; for TIRs/TDRs, default = 980; for MITEs, default = 780.

If indel is enabled, in the extension of seed regions, the following scoring matrix for alignment is used.
--match <int>  Award score (positive number) for 1 match; default = 1.
--mismatch <int>  Penalty score (positive number) for 1 mismatch; default = 1.
--indel <int>  Penalty score (positive number) for 1 indel; default = 2.
--block <int>  Block size during alignment; default = 100.
--block_ratio <float>  For the best alignment in the current block, if the length of aligned sequences <= block_ratio * block_size, the alignment procedure will stop and the end position of the best alignment will be returned. Otherwise, a new block will be created and the alignment will continue from the current end position; default = 0.8.
For MITE detection,
--min_tsd <int>  Minimum length of TSDs; default = 2.
--max_tsd <int>  Maximum length of TSDs; default = 10.

To restrict terminal repeat (TR) and spacer length of detected TIRs, TDRs, and MITE candidates, set the following options:
Note: minimum TR length has been set in the mandatory parameter '--min_tr'.
--max_tr <int>  Maximum TR length.
--min_spacer_len <int>  Minimum spacer length.
--max_spacer_len <int>  Maximum spacer length.

[3.1] TIR detection

Use program "grf-main" with option "-c 0". e.g.,

# grf-main -i genome.fa -o . -c 0 --min_tr 10 -t 16

Here, the input file is "genome.fa"; the output directory is current directory; minimum terminal repeat length = 10; 16 threads are used.

To specify the maximum number of mismatches and indels in TRs, set '--max_mismatch' and '--max_indel' options. To change maximum percentage of unpaired bases, set '-p' option. e.g.,

# grf-main -i genome.fa -o . -c 0 --min_tr 10 -t 16 --max_indel 4 --max_mismatch 4 -p 20

Here, at most 4 indels and 4 mismatches are allowed in TRs, and at most 20% unpaired bases are allowed in TRs (all these criteria must be satisfied at the same time).

Outputs:

"perfect.fasta": TIRs with perfect stems and no spacer in the middle.
"perfect.spacer.fasta": TIRs with perfect stems and a spacer in the middle.
"imperfect.fasta": TIRs with imperfect stems.

The identifier of each TIR sequence (e.g., ">6:2367817:2367964:10m") has the format: ">ChromosomeName:GenomicStartPosition:GenomicStopPosition:TIR_pairing".

For TIR_pairing:
	"m" means matches in base pairing (A and T; C and G);
	"M" means mismatches in base pairing;
	"I" means insertions in base pairing;
	"D" means deletions in base pairing;

[3.2] TDR detection	  
			  
Use program "grf-main" with option "-c 2". e.g.,

# grf-main -i genome.fa -o . -c 2 --min_tr 10 -t 16

Here, the input file is "genome.fa"; the output directory is current directory; minimum terminal repeat length = 10; 16 threads are used.

Outputs:

"perfect.fasta": TDRs with perfect stems and no spacer in the middle.
"perfect.spacer.fasta": TDRs with perfect stems and a spacer in the middle.
"imperfect.fasta": TDRs with imperfect stems.

The identifier of each TDR sequence (e.g., ">6:2367817:2367964:10m") has the format: ">ChromosomeName:GenomicStartPosition:GenomicStopPosition:TDR_alignment".

For TDR alignment:
	"m" means matches in base alignment (A and A; T and T; C and C; G and G);
	"M" means mismatches in base alignment;
	"I" means insertions in base alignment;
	"D" means deletions in base alignment;

[3.3] MITE candidate detection

Use program "grf-main" with option "-c 1". e.g.,

# grf-main -i genome.fa -o . -c 1 --min_tr 10 -t 16

Here, the input file is "genome.fa"; the output directory is current directory; minimum terminal repeat length = 10; 16 threads are used.

Outputs:

"candidate.fasta": MITE candidate sequences.

The identifier of each MITE sequence (e.g., ">6:2367817:2367964:10m:AT") has the format: ">ChromosomeName:GenomicStartPosition:GenomicStopPosition:TIR_pairing:TSD_sequence".

[3.4] MITE family clustering and filtration

(1) Use CD-HIT to cluster similar MITE sequences. e.g.,

# cd-hit-est -i candidate.fasta -o clusteredCandidate.fasta -c 0.80 -n 5 -d 0 -T 16 -aL 0.99 -s 0.8 -M 0 > cd-hit-est.out

Detailed help information of "cd-hit-est" can be found by command:

# cd-hit-est -h
	
(2) Use "grf-cluster" to filter results generated from "cd-hit-est" (i.e., "clusteredCandidate.fasta"):

Usage: grf-cluster [options]
[mandatory parameters]
-i <string>  'cd-hit-est' output file: "*.clstr".
-g <string>  Input genome sequence file in FASTA format.
-o <string>  Output directory.

[optional parameters]
-t <int>  Number of threads used in this program; default = 1.
-f <int>  Length of flanking sequences of MITEs; default = 50.
-c <int>  Minimum copy number of MITEs in the genome; default = 3.
In comparing flanking sequences of MITE candidates, the scoring matrix for alignment is as below:
--match <int>  Award score (positive number) for 1 match; default = 1.
--mismatch <int>  Penalty score (positive number) for 1 mismatch; default = 1.
--indel <int>  Penalty score (positive number) for 1 indel; default  = 2.

e.g., 

# grf-cluster -i clusteredCandidate.fasta.clstr -g genome.fa -o .

Outputs:

"mite.fasta": representative sequences of each MITE family.
"miteSet.fasta": all MITE sequences in each family; each family is separated by dashes.

The identifier of the representative sequence of each MITE family (e.g., ">6:2367817:2367964:10m:AT:10") has the format: ">ChromosomeName:GenomicStartPosition:GenomicStopPosition:TIR_pairing:TSD_seqeunce:CopyNumber".

[3.5] Find nested TIRs/MITEs

Usage: grf-nest <input_fasta> <genome_fasta> <output_fasta>

e.g.,

# grf-nest miteSet.fasta genome.fa miteSet.overlap.fasta

The output file includes groups of nested TIR/MITE sequences; each group is separated by dashes; the first sequence in each group is the longest sequence, and the rest sequences in the group are within the spacer of the largest sequence.

[3.6] Interspersed repeat detection

Usage: grf-intersperse [options]
[mandatory parameters]
-i <string>  Input genome sequence file in FASTA format.
-o <string>  Output directory.

[optional parameters]
-t <int>  Number of threads used in this program; default = 1.
-f <int>  Format for outputs; 0: consensus sequence with position and sequence of each repeat copy; 1: consensus sequence with position of each repeat copy only; default = 0.
-c <int>  Minimum copy number of seeds in the genome; default = 3.
-s <int>  Length of the seed region; default = 20; must >= 10.
-n <int>  Maximum number of undetermined bases in the extension of the seed region (either direction); default = 1.
-p <int>  Minimum identity percentage for a determined base in the consensus sequence; default = 80, which means that any determined base in the consensus sequence must be identical to at least 80% of the bases in the same position from all repeat copies.
-m <int>  Maximum number of mismatches in repeat copies compared with the determined bases of the consensus sequence; repeat copies with mismatches > this value will be removed from the group; default = 2.

e.g.,

# grf-intersperse -i genome.fa -o .

Outputs:

"interspersed_repeat.out": detected interspersed repeats in groups. e.g.,

--------------------------------------------------
>Mt:276777:276835:-
GAGAGGTCCAACGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG
>Mt:86303:86361:+
GAGAGGTCCAAGGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG
>Mt:59244:59302:+
GAGAGGTCCAACGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG

The output represents a group of interspersed repeats with the identifiers (">chromosome:start:end:strand") and sequences.

[3.7] Show Dot-Bracket Notation (DBN) structures of TIRs/MITEs

Usage: grf-dbn <input_fasta> <output_dbn>

e.g.,

# grf-dbn miteSet.fasta miteSet.dbn

Outputs: DBN format of TIR/MITE sequences.

[3.8] Show alignments of TIRs/TDRs/MITEs

Usage: grf-alignment <type> <input_fasta> <output>
type <int>  input sequence type; 1: inverted repeats (or MITEs); 2: direct repeats.

e.g.,

# grf-alignment 1 miteSet.fasta miteSet.alignment

Here are some output examples:

(1) inverted repeats:

>Mt:31842:31861:7m1M2m
CTCTCAGTTTAATCTGAGAG
CTCTCAGTTT
||||||| ||
GAGAGTCTAA

>Mt:328656:328697:4m1M5m1I5m
CTGTAGAGAATGAAGAGGGGCCTAGGATCTTCTTCTCAACAG
CTGTAGAGAATGAAGA
|||| ||||| |||||
GACAACTCTT-CTTCT

(2) direct repeats:

>Mt:64084:64103:5m1M4m
AGGCTGCGGGAGGCTACGGG
AGGCTGCGGG
||||| ||||
AGGCTACGGG

>Mt:66921:66953:10m1I1m1M4m
AGTCTTCTTATCTTAACTTAATGAGTCTTCTTACGTAAC
AGTCTTCTTATCTTAAC
|||||||||| | ||||
AGTCTTCTTA-CGTAAC

[3.9] Filter TIRs/TDRs/MITEs according to given spacer and terminal repeat (TR) lengths

Usage: grf-filter <min_TR_len> <max_TR_len> <min_spacer_len> <max_spacer_len> <input_fasta> <output>

e.g.,

# grf-filter 20 100 10 50 imperfect.fasta imperfect.filtered.fasta

[3.10] Show the alignments and consensus sequences of interspersed repeats

Usage: grf-alignment2 <input> <output>
The input file must be the output file of grf-intersperse (i.e., interspersed_repeat.out) and it must include repeat sequences.

e.g.,

# grf-alignment2 interspersed_repeat.out interspersed_repeat.alignment

Example outputs:

--------------------------------------------------
GAGAGGTCCAACGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG
>Mt:276777:276835:-
GAGAGGTCCAACGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG
mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
>Mt:86303:86361:+
GAGAGGTCCAAGGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG
mmmmmmmmmmmMmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
>Mt:59244:59302:+
GAGAGGTCCAACGTAATTTATTACTCTTATAAAAGAGGGAACTCGACTGAAAGGAGAGG
mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

The output represents a group of interspersed repeats. The first line is the consensus sequence. The rest lines are identifiers of individual repeats (">chromosome:start:end:strand"), their sequences, and their alignments against the consensus sequences. "M" means mismatch; "m" means match.

###############################################################################

[4] How to filter out tandem repeats in outputs using phobos [optional]

[4.1] Installation of phobos

Go to http://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm and download the Linux version of phobos (tested version: 3.3.12).

Extract package:

# tar zxf phobos-v3.3.12-linux.tar.gz

The executable program is "phobos-v3.3.12-linux/bin/phobos-linux-gcc4.1.2".

Optionally, you can add the executable to your $PATH for convenience. e.g., 		

# echo 'export PATH=$PATH:<your_path>/ phobos-v3.3.12-linux/bin' >> ~/.bashrc

To get help information:

# phobos-linux-gcc4.1.2 --help

If error occurs, try to fix it by:

# sudo apt-get install lib32stdc++6

[4.2] Filter out TIRs/TDRs with tandem repeats in TRs

(1) Assuming the file containing TIRs/TDRs in FASTA format generated by "bin/grf-main" is named "input.fa", extract only the TRs (without spacer) using "script/extract_tr.py" in GRF package; the output is named "tr.fa".

# python script/extract_tr.py input.fa tr.fa

Here are example TR sequences in the output:

>Pt:8695:9694:5m1M4m:l
AATAGTTAAA
>Pt:8695:9694:5m1M4m:r
TTTATCTATT

The identifier has the format: ">chrom:start:end:cigar:direction". "l" means the TR is in the 5' end; "r" means the TR is in the 3' end.

(2) Find tandem repeats in TRs using phobos; the output is named "tr.out".

# phobos-linux-gcc4.1.2 --outputFormat 0 tr.fa > tr.out

(3) Filter out TIRs/TDRs with tandem repeats in TRs using "script/filter.py"; the output is named "output.fa".

# python script/filter.py input.fa tr.out output.fa

The output file only includes the TIRs/TDRs without tandem repeats in TRs.

[4.3] Filter out interspersed repeats containing tandem repeats

(1) Covert the interspersed repeat output of GRF (i.e., "interspersed_repeat.out") to FASTA format using "script/convert.py"; the output is named "repeat.fa".

# python script/convert.py interspersed_repeat.out repeat.fa

Here are example sequences in the output:

>1:15739646:15739666:+:1
TGAGGGTGGAGGGGGGGGGGG
>5:13146195:13146215:+:1
TGAGGGTGGAGGGGGGGGGGG
>5:12297814:12297834:+:1
TGAGGGTGGAGGGGGGGGGGG
>5:11365103:11365123:+:1

The identifier has the format ">chrom:start:end:strand:family_number".

(2) Find tandem repeats in interspersed repeats using phobos; the output is named "repeat.out".

# phobos-linux-gcc4.1.2 --outputFormat 0 repeat.fa > repeat.out

(3) Filter out interspersed repeats containing tandem repeats using "script/filter2.py"; the output is named "output.fa".

# python script/filter2.py repeat.fa repeat.out output.fa

(4) Format output file using "script/format.py"; the new output file is named "interspersed_repeat.filtered.out", and it will have the same format with "interspersed_repeat.out".

# python script/format.py output.fa interspersed_repeat.filtered.out

###############################################################################

[5] Modified LTR_FINDER for LTR transposon detection

We modified the source code of LTR_FINDER (https://code.google.com/archive/p/ltr-finder/) to make it accept TDRs generated by "bin/grf-main" as input and perform LTR transposon detection.

The modified source code is in "ltr_finder_modified_src/". The executable program is "bin/ltr_finder".

Notice: all the TDR outputs of GRF (i.e., "perfect.fasta", "perfect.spaer.fasta", and "imperfect.fasta") need to be combined into one file before running the program.

# cat perfect.fasta perfect.spaer.fasta imperfect.fasta > repeat.fasta

Usage: ltr_finder [options] <genome_fasta_file> <TDR_output_of_GRF>

Detailed options can be found by command:

# ltr_finder -h

Optionally, LTR_FINDER can use program "ps_scan" and file "prosite.dat" to find protein domains, and they can be downloaded from the following websites:

ftp://ftp.expasy.org/databases/prosite/ps_scan/ps_scan_linux_x86_elf.tar.gz

ftp://ftp.expasy.org/databases/prosite/prosite.dat

Notice: "prosite.dat" must be put in the folder where "ps_scan" is installed.

LTR_FINDER can also use tRNA database for a specific species in the detection, which can be found in the original LTR_FINDER package ("LTR_FINDER.x86_64-1.0.6/tRNAdb/").

e.g., assuming "ps_scan" is installed in "ps_scan"; Arabidopsis tRNA database is named "LTR_FINDER.x86_64-1.0.6/tRNAdb/Athal-tRNAs.fa"; to detect LTR transposons in Arabidopsis:

# ltr_finder genome.fa repeat.fasta -a ps_scan -s LTR_FINDER.x86_64-1.0.6/tRNAdb/Athal-tRNAs.fa > repeat.out

Here, "repeat.out" is the LTR_FINDER output file including the detected LTR transposons.

###############################################################################


[6] Test data

The folder "tests" includes the genome FASTA file of Saccharomyces_cerevisiae.

To detect TIR in the genome with 8 threads:

# grf-main -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o tir -c 0 --min_tr 10 -t 8

To detect TDR in the genome with 8 threads:

# grf-main -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o tdr -c 2 --min_tr 10 -t 8

To detect MITE candidates in the genome with 8 threads:

# grf-main -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o mite -c 1 --min_tr 10 -t 8

To detect interspersed repeats in the genome with 8 threads:

# grf-intersperse -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o interspersed -t 8
