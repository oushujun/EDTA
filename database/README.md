## Databases
`rice7.0.0.liban` is the manually curated rice TE database. It includes LTR retrotransposons, non-LTR retrotransposons (SINEs and LINEs), DNA transposons (autonomous TIR elements, nonautonomous TIR elements, MITEs, and Helitrons), the rice centromeric repeat (Os1304#Centro/tandem) and the satellite repeat (Os2182#Satellite/rice). In addition, the LTR, SINE, LINE, TIR, and Helitron sublibraries were also provided. This library is obtained from the [riceTElib](https://github.com/oushujun/riceTElib).

`maizeTE11082026` is the Maize TE Consortium (MTEC) curated TE library (version 11/08/2026, 1,517 sequences). It combines three historical versions of the MTEC library, with sequence names converted to the RepeatMasker naming scheme (`name#class/superfamily`) and genic sequences removed using the B73v4 CDS; the rDNA spacer (AF013103.1), the subtelomere 4-12-1 (CL569186.1), and the consensus knob180, TR-1, and CentC sequences provided by Jianing Liu were added. The 2026 updates added 40 curated TEs and strand orientation: each header carries a tab-separated tag (`strand=+` or `strand=?`, with the evidence), and the 844 stranded sequences are given on the (+) strand. This library is obtained from the [MTEC](https://github.com/oushujun/MTEC) repository, which documents the curation steps and the strand evidence of each entry.

`Rice_MSU7.fasta.std7.0.0.out` is the whole-genome TE annotation of [the rice genome (MSUv7)](http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con) based on the `rice7.0.0.liban` library. This command was used to generate the file:
	RepeatMasker -pa 100 -div 40 -cutoff 225 -lib rice7.0.0.liban Rice_MSU7.fasta

`HelitronScanner.training.set.fa` is 100 bp head and tail sequences from curated Helitrons obtained from the [HelitronScanner](https://sourceforge.net/projects/helitronscanner/) package. Sequences were first masked using the `rice6.9.5.liban` library to filter out non-Helitron sequences, then sequence names were renamed to fit the `RepeatMasker` sequence name length requirement.

`Tpases020812DNA` is the DNA TE transposase database obtained from the [LTR_retriever](https://github.com/oushujun/LTR_retriever) package.

`Tpases020812LINE` is the LINE retrotransposase database obtained from the [LTR_retriever](https://github.com/oushujun/LTR_retriever) package.

`alluniRefprexp082813` is the unique plant protein transposase database obtained from the [LTR_retriever](https://github.com/oushujun/LTR_retriever) package. For non-plant species, you may use [uniport_sprot databases](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/) instead.

`dummy060817.fa` just a random sequence copied from the `LTR_retriever` package. It is used to test run `RepeatMasker`.
