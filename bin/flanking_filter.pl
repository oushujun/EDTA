#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
use File::Temp qw(tempdir);
use File::Spec;

my $usage = "\nFilter TE fasta candidates by the repeatiness of their flanking/terminal regions.
	perl flanking_filter.pl -genome genome.fa -query candidate.fa [options]
		-genome	[file]	The multifasta file that used to generate the -query
		-query	[file]	The candidate TE sequence to be filtered by this script
		-extlen	[int]	The length of extended flanking sequence in -query. Default: 30 (bp)
		-tgt_out	[int]	Output taget site with [int] length on each terminal. Default: 15 (bp)
		-miniden	[int]	Minimum identity for flanking sequence alignment. Default: 80 (%)
		-mincov	[float]	Minimum coverage for flanking sequence alignment that counts as full match. Default: 0.8
		-maxct	[int]	Maximum allowed copy number for flanking sequence for a true element. Default: 1.
		-blastplus	[path]	Path to the blastn program. Default: read from \$ENV
		-t|-threads	[int]	Number of threads to run this program. Default: 4
	  Performance/robustness options (results are identical to the per-query version):
		-chunk_size	[int]	Max queries pooled into one blastn (one DB scan amortized over the chunk).
					Default: adaptive ceil(njobs/threads), capped at 1000.
		-route_maxmult	[int]	Queries whose most frequent k-mer (k=word_size) occurs > this many
					times are blasted as singletons so a low-complexity 'detonator' that
					exceeds the timeout isolates itself instead of stalling a whole chunk. Default: 3
		-timeout	[int]	Per-chunk blastn KILL timeout in seconds. Default: 120
	  Checkpoint / resume options:
		-batch_size	[int]	Candidates decided and appended to the tabout per batch (the tabout
					checkpoint granularity). Default: 1000
		-overwrite	[0|1]	0 (default): resume - salvage finished candidates from an existing
					tabout and finished blast counts from the .ckpt sidecar (both optional).
					1: ignore/truncate any existing tabout+ckpt and start fresh.
		-h|-help	Display this help messege and exit.
\n";

my $genome = '';
my $query = '';
my $ext_len = 30;
my $tgt_out = 15;
my $min_iden = 80;
my $min_cov = 0.8;
my $max_ct = 1;
my $max_ct_flank = 50000;
my $blastplus = '';
my $threads = 4;
my $chunk_size = 0;
my $route_maxmult = 3;
my $timeout = 120;
my $word_size = 7;
my $batch_size = 1000;
my $overwrite = 0;          # 0 = resume from existing tabout+ckpt; 1 = start fresh

my $k=0;
foreach (@ARGV){
	$genome = $ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
	$query = $ARGV[$k+1] if /^-query$/i and $ARGV[$k+1] !~ /^-/;
	$ext_len = $ARGV[$k+1] if /^-extlen$/i and $ARGV[$k+1] !~ /^-/;
	$tgt_out = $ARGV[$k+1] if /^-tgt_out$/i and $ARGV[$k+1] !~ /^-/;
	$min_iden = $ARGV[$k+1] if /^-miniden$/i and $ARGV[$k+1] !~ /^-/;
	$min_cov = $ARGV[$k+1] if /^-mincov$/i and $ARGV[$k+1] !~ /^-/;
	$max_ct = $ARGV[$k+1] if /^-maxct$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^-blastplus$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-t$|^-threads$/i and $ARGV[$k+1] !~ /^-/;
	$chunk_size = $ARGV[$k+1] if /^-chunk_size$/i and $ARGV[$k+1] !~ /^-/;
	$route_maxmult = $ARGV[$k+1] if /^-route_maxmult$/i and $ARGV[$k+1] !~ /^-/;
	$timeout = $ARGV[$k+1] if /^-timeout$/i and $ARGV[$k+1] !~ /^-/;
	$batch_size = $ARGV[$k+1] if /^-batch_size$/i and $ARGV[$k+1] !~ /^-/;
	$overwrite = $ARGV[$k+1] if /^-overwrite$/i and $ARGV[$k+1] !~ /^-/;
	die $usage if /^-h$|^-help$/i;
	$k++;
	}

die "The genome file $genome is not found!\n$usage" unless -e $genome or -s "$genome.nsq";  # FASTA, or a prebuilt DB (resume)
die "The query file $query is not found!\n$usage" unless -e $query;

# resolve blastplus path from ENV if not provided
if ($blastplus eq ''){
	chomp ($blastplus=`command -v blastn 2>/dev/null`);
	$blastplus =~ s/blastn$//;
	}
$blastplus="$blastplus/" if $blastplus ne '' and $blastplus !~ /\/$/;

## make blast db for $genome (skip if it already exists - avoids rebuilding a large DB on
## resume, and works when only the DB, not the source FASTA, is present)
`${blastplus}makeblastdb -in $genome -out $genome -dbtype nucl 2> /dev/null` unless -s "$genome.nsq";

my $tabout = "$query.cov${min_cov}iden$min_iden.tabout";
my $ckpt   = "$tabout.ckpt";
my $header = "#Decision\t5'count\t3'count\tflank_count\tChr\tStart\tEnd\tLOC\tflanking\t5'flank\t5'seq\t3'seq\t3'flank\n";

## Store sequence information (all candidates, in input order)
open Query, "<$query" or die $usage;
my @FA;
$/ = "\n>";
while (<Query>){
	chomp;
	s/>//g;
	my ($name, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	$seq = uc $seq;
	push @FA, [$name, $seq];
	}
$/ = "\n";
close Query;

## ---- resume: salvage finished candidates (tabout) and finished blast counts (ckpt) ----
my %done_cand;   # LOC => 1 : candidate already has a complete tabout row -> skip entirely
my %CKPT;        # "LOC\tendtype" => count : a blast already resolved -> reuse, do not re-run
if ($overwrite){
	unlink $tabout, $ckpt;
	}
else {
	if (-s $tabout){
		open my $tf, "<", $tabout or die "Cannot read $tabout: $!\n";
		while (<$tf>){
			next if /^#/;
			chomp;
			my @f = split /\t/;
			next unless @f >= 13;                  # skip a truncated/short (crash) row
			$done_cand{$f[7]} = 1 if defined $f[7] and $f[7] ne '';   # LOC column
			}
		close $tf;
		}
	if (-s $ckpt){
		open my $cf, "<", $ckpt or die "Cannot read $ckpt: $!\n";
		while (<$cf>){
			chomp;
			my ($loc, $et, $ct) = split /\t/;
			next unless defined $ct and $et and ($ct eq 'NA' or $ct =~ /^\d+$/);   # skip malformed/truncated line
			$CKPT{"$loc\t$et"} = $ct;
			}
		close $cf;
		}
	}

## open tabout: append when resuming onto existing content, else create with a header
my $need_header = !(-s $tabout);
open my $OUT, ">>", $tabout or die "Cannot open $tabout: $!\n";
print $OUT $header if $need_header;
$OUT->autoflush(1);

## shared state
my %COUNT :shared;          # job id => copy count (int) or 'NA' (blast did not complete)
my $tmpbase = (-d "/dev/shm" && -w "/dev/shm") ? "/dev/shm" : File::Spec->tmpdir;
my $tmpdir = tempdir("flankf.XXXXXXXX", DIR => $tmpbase, CLEANUP => 1);

## ---- build the work list: candidates not already finished in the tabout ----
my @cand;
for (my $i=0; $i<=$#FA; $i++){
	my ($name, $seq) = @{$FA[$i]}[0,1];
	next unless defined $name;
	my ($chr, $str, $end) = ('','','');
	($chr, $str, $end) = ($1, $2, $3) if $name =~ /^(.*):([0-9]+)\.\.([0-9]+)/;
	my $loc = "$chr:$str..$end";
	next if $chr ne '' and $done_cand{$loc};          # salvaged from a previous run

	my $flank5 = substr $seq, 0, $ext_len;
	my $flank3 = substr $seq, -$ext_len;
	my $seq5 = substr $seq, $ext_len, 30;
	my $seq3 = substr $seq, -($ext_len+30), 30;
	my $tgt_ste = uc((substr $flank5, -$tgt_out)."*".(substr $flank3, 0, $tgt_out));

	# drop candidates whose flanking sequences are simple repeat (unchanged: not printed)
	my ($ssr5, $ssr3) = ('NA', 'NA');
	$ssr5 = &count_base($flank5) if length $flank5 > 0;
	$ssr3 = &count_base($flank3) if length $flank3 > 0;
	next if $ssr5 eq 'true' or $ssr3 eq 'true';

	push @cand, {
		idx=>$i, chr=>$chr, str=>$str, end=>$end, loc=>$loc,
		flank5=>$flank5, seq5=>$seq5, seq3=>$seq3, flank3=>$flank3, tgt=>$tgt_ste,
		end5=>"$flank5$seq5", end3=>"$seq3$flank3", flank=>"$flank5$flank3",
		};
	}

## ---- process the work list in batches; each batch appends its rows to the tabout ----
my $job_id = 0;
for (my $b=0; $b<=$#cand; $b+=$batch_size){
	my $e = $b + $batch_size - 1; $e = $#cand if $e > $#cand;
	my @batch = @cand[$b..$e];

	# stage 1: end5/end3 (reuse ckpt counts where present, else queue for blast)
	my @end_jobs;
	for my $c (@batch){
		for my $side (['e5','end5'], ['e3','end3']){
			my ($tag, $et) = @$side;
			my $id = $job_id++;
			$c->{$tag} = $id;
			my $key = "$c->{loc}\t$et";
			if (exists $CKPT{$key}){ lock(%COUNT); $COUNT{$id} = $CKPT{$key}; }
			else { push @end_jobs, [$id, $c->{$et}, length($c->{$et}), &maxmult($c->{$et}), $key]; }
			}
		}
	&run_jobs(\@end_jobs);

	# stage 2: flank (only for double-repetitive candidates)
	my @flank_jobs;
	for my $c (@batch){
		my $e5 = $COUNT{$c->{e5}}; my $e3 = $COUNT{$c->{e3}};
		next unless $e5 ne 'NA' and $e5 > $max_ct and $e3 ne 'NA' and $e3 > $max_ct;
		my $id = $job_id++;
		$c->{fid} = $id;
		my $key = "$c->{loc}\tflank";
		if (exists $CKPT{$key}){ lock(%COUNT); $COUNT{$id} = $CKPT{$key}; }
		else { push @flank_jobs, [$id, $c->{flank}, length($c->{flank}), &maxmult($c->{flank}), $key]; }
		}
	&run_jobs(\@flank_jobs) if @flank_jobs;

	# stage 3: decide (verbatim logic) + append this batch's rows to the tabout
	for my $c (@batch){
		my $end5_count = $COUNT{$c->{e5}};
		my $end3_count = $COUNT{$c->{e3}};
		my $decision = "true";
		my $end5_repeat = "false";
		($end5_repeat = "true", $decision = "false") if $end5_count ne 'NA' and $end5_count > $max_ct;
		my $end3_repeat = "false";
		($end3_repeat = "true", $decision = "false") if $end3_count ne 'NA' and $end3_count > $max_ct;

		my $flank_count = "NA";
		if ($end5_repeat eq "true" and $end3_repeat eq "true"){
			$flank_count = defined $c->{fid} ? $COUNT{$c->{fid}} : 0;
			$flank_count = 0 if $flank_count eq 'NA';
			if ($flank_count >= 1){
				if (($end5_count + $end3_count)/(2*$flank_count) < 10000 and $end5_count < $max_ct_flank and $end3_count < $max_ct_flank){
					$decision = "true";
					}
				}
			}
		$decision = "false" if $end5_count eq 'NA' or $end3_count eq 'NA';

		print $OUT "$decision\t$end5_count\t$end3_count\t$flank_count\t$c->{chr}\t$c->{str}\t$c->{end}\t$c->{loc}\t$c->{tgt}\t$c->{flank5}\t$c->{seq5}\t$c->{seq3}\t$c->{flank3}\n";
		}
	{ lock(%COUNT); %COUNT = (); }    # free this batch's counts
	}
close $OUT;

## ---- completion: derive pass.fa from the full tabout, then clean up ----
my %pass_loc;
open my $tf, "<", $tabout or die "Cannot read $tabout: $!\n";
while (<$tf>){
	next if /^#/;
	my @f = split /\t/;
	next unless @f >= 13;
	$pass_loc{$f[7]} = 1 if $f[0] eq 'true';
	}
close $tf;
my %TE_cln;
for my $r (@FA){
	my ($name, $seq) = @$r;
	next unless defined $name;
	my ($chr,$str,$end) = ('','','');
	($chr,$str,$end) = ($1,$2,$3) if $name =~ /^(.*):([0-9]+)\.\.([0-9]+)/;
	my $loc = "$chr:$str..$end";
	next unless $pass_loc{$loc};
	$TE_cln{"$chr:$str..$end"} = substr $seq, $ext_len, -$ext_len;
	}
open Seq, ">$query.pass.fa";
foreach my $id (sort {$a cmp $b} keys %TE_cln){ print Seq ">$id\n$TE_cln{$id}\n"; }
close Seq;

unlink $ckpt;                                           # resume no longer needed
`rm $genome.nhr $genome.nin $genome.nsq 2> /dev/null`;  # remove database


## ===================== engine =====================
sub run_jobs {
	my ($jobs) = @_;
	return unless @$jobs;
	my (@big, @singleton);
	for my $j (@$jobs){
		if ($j->[3] > $route_maxmult){ push @singleton, $j; } else { push @big, $j; }
		}
	my $cs = $chunk_size;
	if ($cs <= 0){
		$cs = int((scalar(@$jobs) + $threads - 1)/$threads);
		$cs = 1 if $cs < 1;
		$cs = 1000 if $cs > 1000;
		}
	my @chunks;
	for (my $i=0; $i<=$#big; $i+=$cs){
		my $j = $i+$cs-1; $j = $#big if $j > $#big;
		push @chunks, [ @big[$i..$j] ];
		}
	push @chunks, [$_] for @singleton;
	return unless @chunks;

	my $cq = Thread::Queue->new();
	foreach my $ch (@chunks){
		# chunk data carries [id, seq, qlen, ckptkey] per job
		my @data = map { shared_clone([$_->[0], $_->[1], $_->[2], $_->[4]]) } @$ch;
		$cq->enqueue(shared_clone(\@data));
		}
	$cq->end();

	my $nw = (scalar(@chunks) < $threads) ? scalar(@chunks) : $threads;
	my @workers = map { threads->create(\&chunk_worker, $cq) } (1..$nw);
	$_->join() foreach @workers;
	}

sub chunk_worker {
	my ($cq) = @_;
	my $tid = threads->tid();
	# per-worker append handle to the ckpt; O_APPEND makes small line writes atomic across workers
	open(my $ckfh, '>>', $ckpt) or die "ERROR: cannot append $ckpt: $!\n";
	while (defined(my $chunk = $cq->dequeue())){
		&process_chunk($chunk, $tid, $ckfh);
		}
	close $ckfh;
	}

## Run one chunk (a list of [id,seq,qlen,ckptkey]) as one single-threaded blastn, streaming the output
## through a per-id counter. Each resolved job is recorded to %COUNT and durably checkpointed to .ckpt.
## KILL (timeout/OOM) on a singleton -> NA; on a larger chunk -> bisect so the detonator isolates while
## its neighbours still get counted. A transient (non-KILL) failure is retried up to 3 times.
sub process_chunk {
	my ($jobs, $tid, $ckfh) = @_;
	my $n = scalar @$jobs;
	return unless $n;
	my $qfile = "$tmpdir/c$tid.fa";
	my (%qlen, %key);
	open(my $qfh, '>', $qfile) or die "ERROR: cannot write $qfile: $!\n";
	for my $j (@$jobs){ print $qfh ">$j->[0]\n$j->[1]\n"; $qlen{$j->[0]} = $j->[2]; $key{$j->[0]} = $j->[3]; }
	close $qfh;

	my $cmd = "timeout -s KILL ${timeout}s ${blastplus}blastn -query $qfile -db $genome -outfmt 6 -word_size $word_size -evalue 1e-5 -dust no 2> /dev/null";
	my (%seen, %pass);
	my $rc = -1;
	for (my $attempt=0; $attempt<3; $attempt++){
		%seen = (); %pass = ();
		my $bh;
		unless (open($bh, '-|', $cmd)){ $rc = -1; next; }
		while (my $line = <$bh>){
			my ($q, undef, $iden, $len) = split /\t/, $line;
			next unless defined $len;
			$seen{$q} = 1;
			next unless $iden =~ /^[0-9]/ and $len =~ /^[0-9]/;   # skip a truncated/non-numeric line (e.g. the last line of a KILLed huge-output chunk); count-equivalent to before, without the warning
			$pass{$q}++ if $iden >= $min_iden and $len >= $qlen{$q} * $min_cov;
			}
		close $bh;
		$rc = $?;
		last if $rc == 0;
		last if ($rc >> 8) == 137;    # KILLed -> do not retry
		}

	if ($rc == 0){
		my $rec = '';
		{ lock(%COUNT);
		  for my $j (@$jobs){ my $id = $j->[0]; my $c = $seen{$id} ? ($pass{$id} // 0) : 'NA'; $COUNT{$id} = $c; $rec .= "$key{$id}\t$c\n"; }
		}
		syswrite($ckfh, $rec);        # durable per-chunk checkpoint (atomic append)
		}
	elsif ($n == 1){
		{ lock(%COUNT); $COUNT{$jobs->[0][0]} = 'NA'; }
		syswrite($ckfh, "$key{$jobs->[0][0]}\tNA\n");
		}
	else {
		my $mid = int($n/2);
		&process_chunk([ @{$jobs}[0 .. $mid-1] ], $tid, $ckfh);
		&process_chunk([ @{$jobs}[$mid .. $n-1] ], $tid, $ckfh);
		}
	}

## max k-mer multiplicity (k = word size) within a query: routing metric only, never removes a query.
sub maxmult {
	my ($s) = @_;
	my $L = length $s;
	return 0 if $L < $word_size;
	my %c; my $max = 0;
	for (my $i=0; $i <= $L-$word_size; $i++){
		my $v = ++$c{ substr($s, $i, $word_size) };
		$max = $v if $v > $max;
		}
	return $max;
	}

# determine if the given sequence is simple repeat (unchanged from the per-query version)
sub count_base {
	my $seq = lc $_[0];
	my $repeat = "false";
	$seq =~ s/[nx]+//gi;
	if (length $seq > 0){
		my @base = ('a', 't', 'c', 'g');
		my @count = map { $_ = () = ($seq =~ /$_/gi) } @base;
		@count = (sort { $b<=>$a } @count);
		my $dominant_base = ($count[0] + $count[1])/length $seq;
		$repeat = "true" if $dominant_base >= 0.85;
		} else {
		$repeat = "true";
		}
	return $repeat;
	}
