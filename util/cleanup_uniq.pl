#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "";

# files and parameters
my $genome = $ARGV[0];
my $repeat = $ARGV[1]; #the Red produced .rpt file. To get this file: Red -gnm ./genome/ -rpt ./output/
my $candidate = $ARGV[2]; #The file that you want to have non-repeat sequence removed. MSUID format is required.
my $min_rep_len = 80; #regions in $repeat shorter than this will be removed
my $max_uniq_len = 5000; #non-repeat regions in $candidate longer than this will be removed
my $min_can_len = 80; #after removal of non-repeats, sequence in $candidate shorter than this will be removed
my $threads = 16;

# dependencies
my $combine = "~/las/git_bin/EDTA/util/combine_overlap.pl";
my $substract = "~/las/git_bin/EDTA/util/substract.pl";
my $call_seq = "~/las/git_bin/EDTA/util/call_seq_by_list.pl";
my $RepatMasker = '';

if (1){
# combine repeat regions with filter
`perl $combine $repeat $repeat.cmb`;
`perl -i -nle 'my (\$start, \$end)=(split)[1,2]; my \$len=abs(\$start-\$end)+1; next if \$len<$min_rep_len; print \$_' $repeat.cmb`;

# get candidate regions
`grep \\> $candidate | perl -nle 's/>//; my \$id=(split)[0]; \$id=~s/:/\\t/; \$id=~s/\\.\\./\\t/; print \$id' > $candidate.list`;
`perl $combine $candidate.list $candidate.list.cmb`;

# get non-repeat regions with filter
`perl $substract $candidate.list.cmb $repeat.cmb`;
`perl -nle 'my (\$chr, \$start, \$end)=(split)[0,1,2]; my \$len=abs(\$start-\$end)+1; next if \$len<$max_uniq_len; print "\$chr:\$start..\$end\\t\$chr:\$start..\$end"' $candidate.list.cmb-$repeat.cmb > $candidate.list.uniq`;
`perl $call_seq $candidate.list.uniq -C $genome > $candidate.list.uniq.max$max_uniq_len.fa`;
#`perl $call_seq $candidate.list.uniq -C $genome > $candidate.list.uniq.fa`;
}

# mask uniq regions in $candidate
`${RepatMasker}RepeatMasker -e ncbi -pa $threads -qq -no_is -norna -nolow -div 0.0001 -xsmall -small -lib $candidate.list.uniq.max$max_uniq_len.fa $candidate`;
`perl ~/las/git_bin/EDTA/util/cleanup_tandem.pl -nr 1 -minlen $min_can_len -misschar l -cleanN 1 -cleanT 1 -minrm $max_uniq_len -trf 0 -f $candidate.masked > $candidate.max$max_uniq_len.cln`;

# quick test
`${RepatMasker}RepeatMasker -e ncbi -pa 16 -q -no_is -norna -nolow -div 40 -lib $candidate.max$max_uniq_len.cln -cutoff 225 rice6.9.5.liban.TIR`;
#`${RepatMasker}RepeatMasker -e ncbi -pa 16 -q -no_is -norna -nolow -div 40 -lib $candidate.max$max_uniq_len.cln -cutoff 225 rice6.9.5.liban.Helitron`;
