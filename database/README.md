## Databases
`rice6.9.5.liban` is the manually curated rice TE database. It includes LTR retrotransposons, non-LTR retrotransposons (SINEs and LINEs), DNA transposons (autonomous TIR elements, nonautonomous TIR elements, MITEs, and Helitrons), the rice centromeric repeat (Os1304#Centro/tandem) and the satellite repeat (Os2182#Satellite/rice). An additional of 13 SINE sequences obtained from the [SINE_scan](https://github.com/maohlzj/SINE_Scan) were also added to this library.

`Rice_MSU7.fasta.std6.9.5.out` is the whole-genome TE annotation of [the rice genome (MSUv7)](http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con) based on the `rice6.9.5.liban` library. The parameters were used to generate this file:
	RepeatMasker -pa 36 -no_is -norna -nolow -div 40 -lib rice6.9.5.liban -cutoff 225 Rice_MSU7.fasta.std6.9.5

`HelitronScanner.training.set.fa` is 100 bp head and tail sequences from curated Helitrons obtained from the [HelitronScanner](https://sourceforge.net/projects/helitronscanner/) package. Sequences were first masked using the `rice6.9.5.liban` library to filter out non-Helitron sequences, then sequence names were renamed to fit the `RepeatMasker` sequence name length requirement.

`Tpases020812DNA` is the DNA TE transposase database obtained from the [LTR_retriever](https://github.com/oushujun/LTR_retriever) package.

`Tpases020812LINE` is the LINE retrotransposase database obtained from the [LTR_retriever](https://github.com/oushujun/LTR_retriever) package.

`alluniRefprexp082813` is the unique plant protein transposase database obtained from the [LTR_retriever](https://github.com/oushujun/LTR_retriever) package. For non-plant species, you may use [uniport_sprot databases](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/) instead.

`dummy060817.fa` just a random sequence copied from the `LTR_retriever` package. It is used to test run `RepeatMasker`.
