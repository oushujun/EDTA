genome=Rice_MSU7.fasta
threads=36

### Use stg1 to clean stg0.HQ, then generate stg0.HQ2

# remove mite and helitron in LTR candidates
cat $genome.TIR.fa.stg1 $genome.Helitron.fa.stg1 > $genome.TIR.Helitron.fa.stg1
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR.Helitron.fa.stg1 $genome.LTR.fa.stg0.HQ
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 100 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.stg0.HQ.masked > $genome.LTR.fa.stg0.HQ2

# remove LTR and helitron in TIR candidates
#cat $genome.LTR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ > $genome.LTR.Helitron.fa.stg0.HQ
#RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.Helitron.fa.stg0.HQ $genome.TIR.fa.stg0
#perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg1

# remove LTR in TIR candidates
#cat $genome.LTR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ > $genome.LTR.Helitron.fa.stg0.HQ
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.stg1 $genome.TIR.fa.stg0.HQ
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 100 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.HQ.masked > $genome.TIR.fa.stg0.HQ2

# remove LTR and TIR in Helitron candidates
cat $genome.LTR.fa.stg1 $genome.TIR.fa.stg1 > $genome.LTR.TIR.fa.stg1
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.fa.stg1 $genome.Helitron.fa.stg0.HQ
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 100 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.HQ.masked > $genome.Helitron.fa.stg0.HQ2


### use stg0.HQ2 to clean stg0

# remove mite and helitron in LTR candidates
cat $genome.TIR.fa.stg0.HQ2 $genome.Helitron.fa.stg0.HQ2 > $genome.TIR.Helitron.fa.stg0.HQ2
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.TIR.Helitron.fa.stg0.HQ2 $genome.LTR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.LTR.fa.stg0.masked > $genome.LTR.fa.stg2

# remove LTR and helitron in TIR candidates
#cat $genome.LTR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ > $genome.LTR.Helitron.fa.stg0.HQ
#RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.Helitron.fa.stg0.HQ $genome.TIR.fa.stg0
#perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 1 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg1

# remove LTR in TIR candidates
#cat $genome.LTR.fa.stg0.HQ $genome.Helitron.fa.stg0.HQ > $genome.LTR.Helitron.fa.stg0.HQ
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.fa.stg0.HQ2 $genome.TIR.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.TIR.fa.stg0.masked > $genome.TIR.fa.stg2

# remove LTR and TIR in Helitron candidates
cat $genome.LTR.fa.stg0.HQ2 $genome.TIR.fa.stg0.HQ2 > $genome.LTR.TIR.fa.stg0.HQ2
RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $genome.LTR.TIR.fa.stg0.HQ2 $genome.Helitron.fa.stg0
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 0 -cleanN 1 -cleanT 1 -f $genome.Helitron.fa.stg0.masked > $genome.Helitron.fa.stg2

# aggregate clean sublibraries and cluster
cat $genome.LTR.fa.stg2 $genome.TIR.fa.stg2 $genome.Helitron.fa.stg2 > $genome.LTR.TIR.Helitron.fa.stg2
