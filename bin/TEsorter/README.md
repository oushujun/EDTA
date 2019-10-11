# TEsorter
It is coded for [LTR_retriever](https://github.com/oushujun/LTR_retriever) to classify long terminal repeat retrotransposons (LTR-RTs) at first. It can also be used to classify any other TE sequences, including Class I and Class II elements which are covered by the [REXdb](http://repeatexplorer.org/?page_id=918) database.
  
For more details of methods and benchmarking results, see the [preprint paper](https://doi.org/10.1101/800177).

### Installation ###
Dependencies:
+    [python 2.7](https://www.python.org/)  
	+   [biopython](https://biopython.org/): quickly install by `pip install biopython`  
	+   [parallel python](https://www.parallelpython.com/): quickly install by `pip install pp`
+    [hmmscan 3.1x or 3.2x](http://hmmer.org/): be compatible with HMMER3/f database format.
+   [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 
```
git clone https://github.com/zhangrengang/TEsorter
```
Note: do not move or hard link `TEsorter.py` alone to anywhere else, as it rely on `database/` and `bin/`. You can add the directory to `PATH` or soft link `TEsorter.py` to `PATH`.

### Quick Start ###
```
git clone https://github.com/zhangrengang/TEsorter
cd TEsorter
sh build_database.sh 
cd test

# run
python ../TEsorter.py rice6.9.5.liban
```
By default, the newly released [REXdb](http://repeatexplorer.org/?page_id=918) ([viridiplantae_v3.0 + metazoa_v3](https://bitbucket.org/petrnovak/re_databases)) database is used, which is more sensitive and more common and thus is recommended. 
  
Classical [GyDB](http://gydb.org/) can also be used:
```
python ../TEsorter.py rice6.9.5.liban -db gydb
```
To speed up, use more processors [default=4]:
```
python ../TEsorter.py rice6.9.5.liban -p 20
```
To improve sensitivity, reduce the criteria (coverage and E-value):
```
python ../TEsorter.py rice6.9.5.liban -p 20 -cov 10 -eval 1e-2
```
To improve specificity, increase the criteria and disable the pass2 mode:
```
python ../TEsorter.py rice6.9.5.liban -p 20 -cov 30 -eval 1e-5 -dp2
```
To improve sensitivity of pass-2, reduce the rule:
```
python ../TEsorter.py rice6.9.5.liban -p 20 -rule 70-30-80
```
To classify TE polyprotein sequences ([an example](http://www.repeatmasker.org/RMDownload.html)):
```
python ../TEsorter.py RepeatPeps.lib -st prot -p 20
```
### Outputs ###
```
rice6.9.5.liban.rexdb.domtbl        HMMScan raw output
rice6.9.5.liban.rexdb.dom.faa       protein sequences of domain, which can be used for phylogenetic analysis.
rice6.9.5.liban.rexdb.dom.tsv       inner domains of LTRs, which might be used to filter domains based on their scores and coverages.
rice6.9.5.liban.rexdb.dom.gff3      domain annotations in `gff3` format
rice6.9.5.liban.rexdb.cls.tsv       TEs/LTRs classifications
    Column 1: raw id
    Column 2: Order, e.g. LTR
    Column 3: Superfamily, e.g. Copia
    Column 4: Clade, e.g. SIRE
    Column 5: Complete, "yes" means one LTR Copia/Gypsy element with full GAG-POL domains.
    Column 6: Strand, + or - or ?
    Column 7: Domains, e.g. GAG|SIRE PROT|SIRE INT|SIRE RT|SIRE RH|SIRE; `none` for pass-2 classifications
rice6.9.5.liban.rexdb.cls.lib       fasta library for RepeatMakser
rice6.9.5.liban.rexdb.cls.pep       the same sequences as `rice6.9.5.liban.rexdb.dom.faa`, but id is changed with classifications.
```

### Usage ###
```
$ python TEsorter.py  -h
usage: TEsorter.py [-h] [-v] [-db {rexdb,rexdb-plant,rexdb-metazoa,gydb}]
                   [-st {nucl,prot}] [-pre PREFIX] [-fw] [-p PROCESSORS]
                   [-tmp TMP_DIR] [-cov MIN_COVERAGE] [-eval MAX_EVALUE]
                   [-dp2] [-rule PASS2_RULE] [-nolib] [-norc] [-nocln]
                   sequence

positional arguments:
  sequence              input TE sequences in fasta format [required]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -db {rexdb,rexdb-plant,rexdb-metazoa,gydb}, --hmm-database {rexdb,rexdb-plant,rexdb-metazoa,gydb}
                        the database used [default=rexdb]
  -st {nucl,prot}, --seq-type {nucl,prot}
                        'nucl' for DNA or 'prot' for protein [default=nucl]
  -pre PREFIX, --prefix PREFIX
                        output prefix [default='{-s}.{-db}']
  -fw, --force-write-hmmscan
                        if False, will use the existed hmmscan outfile and
                        skip hmmscan [default=False]
  -p PROCESSORS, --processors PROCESSORS
                        processors to use [default=4]
  -tmp TMP_DIR, --tmp-dir TMP_DIR
                        directory for temporary files [default=./tmp]
  -cov MIN_COVERAGE, --min-coverage MIN_COVERAGE
                        mininum coverage for protein domains in HMMScan output
                        [default=20]
  -eval MAX_EVALUE, --max-evalue MAX_EVALUE
                        maxinum E-value for protein domains in HMMScan output
                        [default=0.001]
  -dp2, --disable-pass2
                        do not further classify the unclassified sequences
                        [default=False for `nucl`, True for `prot`]
  -rule PASS2_RULE, --pass2-rule PASS2_RULE
                        classifying rule [identity-coverage-length] in pass-2
                        based on simliarity [default=80-80-80]
  -nolib, --no-library  do not generate a library file for RepeatMasker
                        [default=False]
  -norc, --no-reverse   do not reverse complement sequences if they are
                        detected in minus strand [default=False]
  -nocln, --no-cleanup  do not clean up the temporary directory
                        [default=False]
```

### Limitations ###
1. For each domain (e.g. RT), only the best hit with the highest score will output, which means: 1) if frame is shifted, only one part can be annotated; 2) for example, if two or more RT domains are present in one query sequence, only one of these RT domains will be annotated.
2. Many LTR-RTs cannot be classified due to no hit, which might be because: 1) the database is still incompleted; 2) some LTR-RTs may have too many mutations such as frame shifts and stop gains or have lost protein domains; 3) some LTR-RTs may be false positive. For the test data set ([rice6.9.5.liban](https://raw.githubusercontent.com/oushujun/EDTA/master/database/rice6.9.5.liban)), ~84% LTR-RTs (_INT sequences) are classified.
3. Non-autonomous TEs that lack protein domains, some un-active autonomous TEs that have lost their protein domains and any other elements that contain none protein domains, are excepted to be un-classified.

### Further phylogenetic analyses ###
You may want to use the RT domains to analysis relationships among retrotransposons (LTR, LINE, DIRS, etc.). Here is an example (with [mafft](https://mafft.cbrc.jp/alignment/software/) and [iqtree](http://www.iqtree.org/) installed):
```
# to extract RT domain sequences
cat rice6.9.5.liban.rexdb.dom.tsv | grep -P "\-RT\t" | python ../bin/get_record.py -i rice6.9.5.liban.rexdb.dom.faa -o rice6.9.5.liban.rexdb.dom.RT.faa -t fasta

# to align with MAFFT or other tools
mafft --auto rice6.9.5.liban.rexdb.dom.RT.faa > rice6.9.5.liban.rexdb.dom.RT.aln

# to reconduct the phylogenetic tree with IQTREE or other tools
iqtree -s rice6.9.5.liban.rexdb.dom.RT.aln -bb 1000 -nt AUTO 

# Finally, visualize and edit the tree 'rice6.9.5.liban.rexdb.RT.faa.aln.treefile' with FigTree or other tools.
```
The alignments can also be generated by:
```
python ../bin/concatenate_domains.py rice6.9.5.liban.rexdb.cls.pep RT > rice6.9.5.liban.rexdb.cls.pep.RT.aln
```
The alignments of LTR-RTs full domains can be generated by (align and concatenate):
```
python ../bin/concatenate_domains.py rice6.9.5.liban.rexdb.cls.pep GAG PROT RH RT INT > rice6.9.5.liban.rexdb.cls.pep.full.aln
```
The alignments of Class I INT and Class II TPase (DDE-transposases) can be generated by:
```
python ../bin/concatenate_domains.py rice6.9.5.liban.rexdb.cls.pep INT > rice6.9.5.liban.rexdb.cls.pep.INT.aln
python ../bin/concatenate_domains.py rice6.9.5.liban.rexdb.cls.pep TPase > rice6.9.5.liban.rexdb.cls.pep.TPase.aln
cat rice6.9.5.liban.rexdb.cls.pep.INT.aln rice6.9.5.liban.rexdb.cls.pep.TPase.aln > rice6.9.5.liban.rexdb.cls.pep.INT_TPase.faa
mafft --auto rice6.9.5.liban.rexdb.cls.pep.INT_TPase.faa > rice6.9.5.liban.rexdb.cls.pep.INT_TPase.aln
```
Note: the domain names between rexdb and gydb are different: PROT (rexdb) = AP (gydb), RH (rexdb) = RNaseH (gydb). You should use the actual domain name.

### Extracting TE sequences from genome for TEsorter ###
Here are examples to extract TE sequences from outputs of wide-used softwares, when you have only genome sequences.

1. extract all TE sequences from [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) output:
```
# run RepeatMasker, which will generate a *.out file.
RepeatMasker [options] genome.fa

# extract sequences
python [path_to_TEsorter]/bin/RepeatMasker.py out2seqs genome.fa.out genome.fa > whole_genome_te.fa

# classify
python [path_to_TEsorter]/TEsorter.py whole_genome_te.fa [options]
```

2. extract all intact LTR-RTs sequences from [LTR_retriever](https://github.com/oushujun/LTR_retriever) outputs:
```
# run LTR_retriever, which generate two *.pass.list files.
LTR_retriever -genome genome.fa [options]

# extract sequences
python [path_to_TEsorter]/bin/LTR_retriever.py get_full_seqs genome.fa > intact_ltr.fa

# classify
python [path_to_TEsorter]/TEsorter.py intact_ltr.fa [options]
```
