## ~ ~ ~ Run LTR_HARVEST in parallel ~ ~ ~
This is a Perl wrapper for LTR_harvest modified from [LTR_FINDER_parallel](https://github.com/oushujun/LTR_FINDER_parallel). All rights reserved to the original authors.


### Prerequisite: [GenomeTools](https://github.com/genometools/genometools)
`conda create -n LTR_HARVEST_parallel`  
`conda install -n LTR_HARVEST_parallel -c bioconda -y genometools-genometools`  
`conda activate LTR_HARVEST_parallel`


### Installation of LTR_HARVEST_parallel: No need. Just download and run.

	Usage: perl LTR_HARVEST_parallel -seq [file] -size [int] -threads [int]  
	Options:
		-seq    [file]  Specify the sequence file.
		-size   [int]   Specify the size you want to split the genome sequence.
				Please make it large enough to avoid spliting too many LTR elements. Default 5000000 (bp).  			 
		-time   [int]   Specify the maximum time to run a subregion (a thread).
				This helps to skip simple repeat regions that take a substantial of time to run. Default: 1500 (seconds).
				Suggestion: 300 for -size 1000000. Increase -time when -size increased.  
		-try1   [0|1]   If a region requires more time than the specified -time (timeout), decide:  
					0, discard the entire region.
					1, further split to 50 Kb regions to salvage LTR candidates (default);
		-next           Only summarize the results for previous jobs without rerunning LTR_HARVEST (for -v).
		-verbose|-v     Retain LTR_HARVEST outputs for each sequence piece.
		-threads|-t     [int]   Indicate how many CPU/threads you want to run LTR_HARVEST.
		-check_dependencies Check if dependencies are fullfiled and quit
		-help|-h        Display this help information.


### Input
Genome file in multi-FASTA format.


### Output
GFF3, LTRharvest (STDOUT) formats of predicted LTR candidates.
