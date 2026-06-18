########
# EDTA #
########

## Test genome (genome.fa)
A ~1.18 Mb toy genome with TWO sequences that together exercise all five of EDTA's
de-novo structural detectors:

  Chr2_rice_Nipponbare_MSU7   1,000,000 bp   real rice (Oryza sativa) chr2 segment,
                                             unchanged. Yields de-novo LTR, TIR, Helitron.
  ChrSyn                        183,094 bp   synthetic TE panel adding the families rice
                                             chr2 lacks at de-novo-detectable copy number
                                             and structure:
                                               ~20 LINE copies   (from rice Os0651, L1; RT+EN)
                                               ~10 SINE copies   (from rice Os3708 / OsSN3, tRNA)
                                                1 intact LTR retrotransposon (Copia, from Os0008)
                                                1 intact Helitron (from Os0906)
                                             ...separated by unique complex spacer DNA, plus
                                             three nested insertions that also exercise EDTA's
                                             cross-type contamination purge (no "No sequences
                                             were masked" notices):
                                               - LINE fragment nested in the LTR internal region
                                               - chimeric SINE carrying a LINE-derived fragment
                                               - SINE nested in the Helitron body

Chr2 is byte-identical to the previous test genome, so genome.cds.fa and
genome.exclude.bed (coordinates referencing chr2 only) remain valid.

Ground truth for every engineered element (coordinates, family, divergence, TSD):
  genome.synthetic.manifest.tsv
Deterministic, reproducible build (same seed -> same genome):
  python3 build_test_genome.py --refs <dir-with-source-FASTAs> --out . --seed 20260616

## IMPORTANT — conda environment for SINE detection
De-novo SINE detection (AnnoSINE -> nhmmer) needs libopenblas. An EDTA env built with the
MKL BLAS backend (e.g. an MKL-pinned yml) lacks libopenblas.so.0, so nhmmer fails to load
and finds ZERO SINEs for ANY genome -- and because EDTA runs AnnoSINE with `2>/dev/null`
(EDTA_raw.pl), this surfaces only as the misleading "The SINE result file has 0 bp!".
Use an env built from the openblas-based EDTA_2.3.yml, or patch an existing env:
  mamba install -n <env> libopenblas

## Test with the conda version
conda activate EDTA   # an env with libopenblas (see note above)
nohup /usr/bin/time -v perl ../EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 > EDTA.test &


## Test with the Docker version
nohup docker run -v $PWD:/in -w /in oushujun/edta:2.0.0 EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 > EDTA.test &


## Test with the Singularity version
nohup singularity exec EDTA.sif EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 > EDTA.test &

###########
# panEDTA #
###########
nohup sh ../panEDTA.sh -g genome.cds.list -c genome.cds.fa -l ../database/athrep.updated.nonredun.fasta -t 20 -f 3 &

## Expected result
All five structural detectors (LTR, SINE, LINE, TIR, Helitron) produce non-empty raw
libraries, and the whole-genome annotation (genome.fa.mod.EDTA.TEanno.sum) reports all
five classes. The run is clean: all four notices the previous toy genome emitted are
gone -- "The SINE result file has 0 bp!" and "The LINE result file has 0 bp!" (SINE/LINE
now detected de novo), "LOC list ...ltrTE.veryfalse is empty" (LTR_retriever), and
"No sequences were masked" (the three internal cross-type purges now each find the
engineered nested contamination).
