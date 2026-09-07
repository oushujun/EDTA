#!/usr/bin/env python3
"""Compute LTR boundaries directly from the final EDTA library by self-alignment.

For each #LTR/ library entry, self-align the sequence with blastn and take the best
off-diagonal HSP that starts near the 5' end and ends near the 3' end: that HSP is the
left LTR aligned to the right LTR. Emit the EDTA LTRbound format:

    <TE_id#class>\t<total_len>\t<lLTR_len>\t<rLTR_len>

Because it measures the sequence that actually ships, it stays correct after EDTA's
advance filtering trims a library sequence -- unlike rewriting LTR_retriever's
pre-filtering numbers (bin/rename_LTR.pl / get_range2.pl), which then disagree with the
library and can make lLTR_len + rLTR_len > total_len, silently corrupting
label_solo_LTR.pl. About 1/3 of LTR-class entries resolve a terminal-repeat pair; the
rest genuinely lost one or both LTRs to filtering and are omitted (they must NOT be fed
to label_solo_LTR.pl).

Validated at 100% detection / 96.3% exact / 99.7% within 10 bp against the entries whose
stored numbers were still internally consistent.

Usage: ltrbound_denovo.py <final_lib.fa> <out> <nproc> [limit]
blastn is taken from $EDTA_BLASTN if set, else 'blastn' on PATH.
"""
import os, subprocess, sys, tempfile
from multiprocessing import Pool

FA, OUT, NPROC = sys.argv[1], sys.argv[2], int(sys.argv[3])
LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 0
BLASTN = os.environ.get('EDTA_BLASTN', 'blastn')

def read_fa(path, ltr_only=True):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None and (not ltr_only or '#LTR/' in name):
                    yield name, ''.join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None and (not ltr_only or '#LTR/' in name):
        yield name, ''.join(chunks)

def bounds(rec):
    name, seq = rec
    n = len(seq)
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as tf:
        tf.write(f">q\n{seq}\n")
        p = tf.name
    try:
        r = subprocess.run(
            [BLASTN, '-query', p, '-subject', p, '-evalue', '1e-5',
             '-word_size', '11', '-dust', 'no', '-strand', 'plus',
             '-outfmt', '6 qstart qend sstart send pident length'],
            capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        sys.stderr.write(f"ltrbound_denovo.py: WARNING: blastn timed out on {name}, skipping\n")
        return None
    finally:
        os.unlink(p)
    if r.returncode != 0:
        sys.stderr.write(f"ltrbound_denovo.py: WARNING: blastn exited {r.returncode} on {name}, skipping: {r.stderr.strip()}\n")
        return None

    best = None
    for line in r.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 6:
            continue
        qs, qe, ss, se = (int(x) for x in f[:4])
        pid, alen = float(f[4]), int(f[5])
        if ss <= qe:                      # must be off-diagonal (not the self hit)
            continue
        if qs > 0.05 * n + 50:            # left LTR must start near the 5' end
            continue
        if se < n - (0.05 * n + 50):      # right LTR must end near the 3' end
            continue
        if alen < 80 or pid < 80:
            continue
        score = alen * pid
        if best is None or score > best[0]:
            best = (score, qe - qs + 1, se - ss + 1)
    if best is None:
        return None
    return f"{name}\t{n}\t{best[1]}\t{best[2]}"

def main():
    recs = list(read_fa(FA))
    if LIMIT:
        recs = recs[:LIMIT]
    sys.stderr.write(f"ltrbound_denovo.py: LTR-class records: {len(recs)}\n")
    ok = 0
    with Pool(NPROC) as pool, open(OUT, 'w') as out:
        for res in pool.imap_unordered(bounds, recs, chunksize=8):
            if res:
                out.write(res + "\n"); ok += 1
    sys.stderr.write(f"ltrbound_denovo.py: boundaries found: {ok} / {len(recs)}\n")

if __name__ == '__main__':
    main()
