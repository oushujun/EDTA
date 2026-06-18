########
# EDTA #
########

## Test with the conda version
conda activate EDTA
nohup /usr/bin/time -v perl ../EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 > EDTA.test &

## Test with the Docker version
nohup docker run -v $PWD:/in -w /in oushujun/edta:2.0.0 EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 > EDTA.test &

## Test with the Singularity version
nohup singularity exec EDTA.sif EDTA.pl --genome genome.fa --cds genome.cds.fa --curatedlib ../database/rice7.0.0.liban --exclude genome.exclude.bed --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10 > EDTA.test &

## Expected result
All five structural detectors (LTR, SINE, LINE, TIR, Helitron) produce non-empty raw
libraries, and the whole-genome annotation (genome.fa.mod.EDTA.TEanno.sum) reports all
five classes.

###########
# panEDTA #
###########
nohup sh ../panEDTA.sh -g genome.cds.list -c genome.cds.fa -l ../database/athrep.updated.nonredun.fasta -t 20 -f 3 &

