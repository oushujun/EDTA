#!/usr/bin/env python3
"""
build_test_genome.py - deterministically build the "ChrSyn" synthetic TE contig
that augments EDTA's toy test genome (rice Chr2) so all five structural
detectors fire and the SINE/LINE "0 bp" + "No sequences were masked" warnings
disappear.

ChrSyn carries, separated by unique complex spacer DNA:
  * ~20 dispersed LINE copies of a master verified (TEsorter) to carry RT+EN
        -> RepeatModeler2 models a family TEsorter classifies LINE (warning #3)
  * ~10 dispersed SINE copies (tRNA head + LINE-3'UTR tail + poly-A + TSD)
        -> AnnoSINE_v2 detects them (warning #2); the LINE tail also makes the
           SINE raw lib homologous to the LINE lib (EDTA_processK.pl:182)
  * 1 intact LTR retrotransposon (Copia) with a LINE 5'UTR fragment nested in
        its internal region -> LTR raw lib homologous to LINE lib (:170)
  * 1 intact Helitron with a full SINE copy nested in its body
        -> TIR/Helitron raw lib homologous to SINE lib (:200)

All three masking-purge steps therefore find hits => no "No sequences were masked".

The rice Chr2 sequence is NOT touched here; the caller concatenates Chr2 + ChrSyn.

Deterministic: seeded RNG only; no wall-clock / urandom. Same seed => same genome.
"""
import argparse
import os
import random

BASES = "ACGT"
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def first_seq(path):
    d = read_fasta(path)
    if not d:
        raise SystemExit(f"ERROR: no sequence in {path}")
    k = next(iter(d))
    return k, d[k]


def revcomp(s):
    return s.translate(_COMP)[::-1]


def rand_dna(n, rng, gc=0.5):
    """Complex random DNA, ~gc GC content, no homopolymer run > 5."""
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]  # A,C,G,T
    out, run = [], 0
    for _ in range(n):
        b = rng.choices(BASES, weights=weights)[0]
        if out and b == out[-1]:
            run += 1
            if run >= 5:
                b = rng.choice([x for x in BASES if x != b])
                run = 0
        else:
            run = 0
        out.append(b)
    return "".join(out)


def mutate(seq, rate, rng):
    """i.i.d. point substitutions only (length preserved => consensus == master)."""
    s = list(seq)
    for i, c in enumerate(s):
        if c in BASES and rng.random() < rate:
            s[i] = rng.choice([x for x in BASES if x != c])
    return "".join(s)


def make_tsd(n, rng):
    """Mixed-base TSD, never pure A/T (AnnoSINE rejects pure-A/T TSDs)."""
    while True:
        t = "".join(rng.choice(BASES) for _ in range(n))
        if any(c in "GC" for c in t):
            return t


def utr3_window(master, length, rng):
    """An ~`length` bp window from the 3' UTR, skipping any terminal A/T homopolymer
    so the fragment is complex (not poly-A) yet still LINE-homologous."""
    end = len(master)
    while end > length and master[end - 1] in "AT":
        end -= 1
    start = max(0, end - length)
    return master[start:end]


# ---------------------------------------------------------------------------
# element builders -> each returns dict(seq, name, type, family, strand,
#                     divergence, tsd, nested=[{rel_start,rel_end,...}])
# ---------------------------------------------------------------------------
def build_line_copies(master, family, k, div, rng, tail_lo=25, tail_hi=40):
    els = []
    for i in range(k):
        body = mutate(master, div, rng) + ("A" * rng.randint(tail_lo, tail_hi))
        strand = rng.choice("+-")
        core = revcomp(body) if strand == "-" else body
        tsd = make_tsd(rng.randint(12, 15), rng)
        els.append(dict(seq=tsd + core + tsd, name=f"LINE_{i:02d}", type="LINE",
                        family=family, strand=strand, divergence=div, tsd=tsd,
                        nested=[]))
    return els


def build_sine_copies(master, family, line_frag, k, div, rng, n_chimeric=4, ins_pos=90,
                      tail_lo=12, tail_hi=20):
    """Pure SINE = TSD + Os3708 + polyA + TSD. The poly-A tail and TSD sit ADJACENT to
    the SINE body end so AnnoSINE's TSD search (which keys off the detected SINE 3'
    boundary) succeeds -- the earlier design appended a LINE 3' tail that pushed the TSD
    ~100 bp past the boundary and AnnoSINE dropped every copy.

    Chimeric copies (first n_chimeric) insert a short non-coding LINE fragment INTO the
    body (after the head, before the natural tail), giving the SINE library LINE homology
    (-> EDTA_processK.pl:182) WITHOUT displacing the 3' TSD."""
    els = []
    for i in range(k):
        m = mutate(master, div, rng)
        chim = i < n_chimeric
        if chim:
            frag = mutate(line_frag, div, rng)
            body = m[:ins_pos] + frag + m[ins_pos:]
        else:
            body = m
        core = body + ("A" * rng.randint(tail_lo, tail_hi))
        strand = rng.choice("+-")
        oriented = revcomp(core) if strand == "-" else core
        tsd = make_tsd(rng.randint(12, 14), rng)
        nested = ([{"type": "LINE_frag", "note": "LINE frag in SINE body (->:182)"}]
                  if chim else [])
        els.append(dict(seq=tsd + oriented + tsd,
                        name=f"SINE_{i:02d}{'c' if chim else ''}", type="SINE",
                        family=family, strand=strand, divergence=div, tsd=tsd,
                        nested=nested))
    return els


def build_one_sine(master, family, div, rng):
    """A single pure SINE copy (Helitron-nested cargo -> :200 SINE homology)."""
    body = mutate(master, div, rng) + ("A" * 16)
    tsd = make_tsd(13, rng)
    return tsd + body + tsd


def build_ltr_host(ltr, intr, line5utr, family, rng):
    """LTR + (INT with LINE-5'UTR fragment nested mid-internal) + LTR, flanked by 5bp TSD."""
    mid = len(intr) // 2
    intr_n = intr[:mid] + line5utr + intr[mid:]
    elem = ltr + intr_n + ltr  # two identical LTRs => 100% LTR-LTR identity
    tsd = make_tsd(5, rng)
    seq = tsd + elem + tsd
    nstart = len(tsd) + len(ltr) + mid          # 0-based start of nested frag in `seq`
    return dict(seq=seq, name="LTRRT_host", type="LTR", family=family, strand="+",
                divergence=0.0, tsd=tsd,
                nested=[{"type": "LINE_frag", "rel_start": nstart,
                         "rel_end": nstart + len(line5utr),
                         "note": "LINE-5UTR nested in LTR internal (->:170)"}])


def build_helitron_host(hel_seq, sine_cargo, family, rng, req_left="A", req_right="T"):
    """Full real Helitron (Os0906) with a SINE nested in its body, away from the TC head
    and the CTRR-hairpin tail (the termini HelitronScanner's LCV models need intact).
    Using the whole real element gives reliable detection (a reconstructed head+tail did
    not); the body-nested SINE then sits inside the detected element -> :200 hit."""
    # HelitronScanner detects only the 3' ~3.6 kb (anchors on the CTRR tail, extends back
    # to a head signal ~3.6 kb upstream). Insert the SINE ~1.7 kb before the 3' terminus
    # so it lands INSIDE that detected window, clear of the hairpin.
    ins = len(hel_seq) - 1700
    seq = hel_seq[:ins] + sine_cargo + hel_seq[ins:]
    return dict(seq=seq, name="Helitron_host", type="Helitron", family=family,
                strand="+", divergence=0.0, tsd="none",
                req_left=req_left, req_right=req_right,
                nested=[{"type": "SINE_frag", "rel_start": ins,
                         "rel_end": ins + len(sine_cargo),
                         "note": "SINE nested in Helitron body (->:200)"}])


def assemble(elements, rng, sp_lo=300, sp_hi=500):
    """Concatenate elements with unique complex spacers. 1-based inclusive coords."""
    parts, rows, pos = [], [], 0

    def add_spacer(force_last=None, force_first=None):
        nonlocal pos
        sp = rand_dna(rng.randint(sp_lo, sp_hi), rng)
        if force_last:
            sp = sp[:-1] + force_last
        if force_first:
            sp = force_first + sp[1:]
        parts.append(sp)
        pos += len(sp)

    add_spacer()
    for el in elements:
        # Helitron needs host 'A' immediately 5' (fix tail of preceding spacer)
        if el.get("req_left"):
            parts[-1] = parts[-1][:-1] + el["req_left"]
        start = pos + 1                      # 1-based inclusive
        parts.append(el["seq"]); pos += len(el["seq"])
        end = pos
        rows.append(dict(start=start, end=end, **{k: el[k] for k in
                    ("name", "type", "family", "strand", "divergence", "tsd")}))
        for nf in el.get("nested", []):
            if "rel_start" in nf:
                rows.append(dict(start=start + nf["rel_start"], end=start + nf["rel_end"] - 1,
                                 name=el["name"] + ":" + nf["type"], type=nf["type"],
                                 family="(nested)", strand=el["strand"],
                                 divergence="", tsd="", note=nf.get("note", "")))
        add_spacer(force_first=el.get("req_right"))
    return "".join(parts), rows


def write_fasta(path, name, seq, width=60):
    with open(path, "w") as o:
        o.write(f">{name}\n")
        for i in range(0, len(seq), width):
            o.write(seq[i:i + width] + "\n")


def write_manifest(path, contig, rows):
    cols = ["contig", "start", "end", "name", "type", "family", "strand",
            "divergence", "tsd", "note"]
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r.get(c, "")) if c != "contig" else contig
                              for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", default="refs")
    ap.add_argument("--out", default="build")
    ap.add_argument("--seed", type=int, default=20260616)
    ap.add_argument("--contig", default="ChrSyn")
    ap.add_argument("--n-line", type=int, default=20)
    ap.add_argument("--n-sine", type=int, default=10)
    ap.add_argument("--line-div", type=float, default=0.05)
    ap.add_argument("--sine-div", type=float, default=0.03)
    ap.add_argument("--line5-len", type=int, default=400)
    ap.add_argument("--line3-len", type=int, default=80)
    # source files (created in Tasks 0-1)
    ap.add_argument("--line-master", default=None, help="default <refs>/master_line.fa")
    ap.add_argument("--sine-master", default=None, help="default <refs>/master_sine.fa")
    ap.add_argument("--ltr", default=None, help="default <refs>/ltr_os0008.fa")
    ap.add_argument("--helitron", default=None, help="default <refs>/helitron_os0906.fa")
    a = ap.parse_args()

    R = a.refs
    line_master_p = a.line_master or os.path.join(R, "master_line.fa")
    sine_master_p = a.sine_master or os.path.join(R, "master_sine.fa")
    ltr_p = a.ltr or os.path.join(R, "ltr_os0008.fa")
    hel_p = a.helitron or os.path.join(R, "helitron_os0906.fa")

    rng = random.Random(a.seed)

    line_fam, line_seq = first_seq(line_master_p)
    sine_fam, sine_seq = first_seq(sine_master_p)
    line5utr = line_seq[:a.line5_len]
    line_sine_frag = line_seq[150:230]   # 80 bp non-coding 5'UTR LINE fragment for SINE body

    ltr_d = read_fasta(ltr_p)
    ltr = next(v for k, v in ltr_d.items() if "_LTR" in k)
    intr = next(v for k, v in ltr_d.items() if "_INT" in k)

    hel_seq = first_seq(hel_p)[1]  # full real Helitron Os0906 (intact TC head + CTRR tail)

    # ---- build elements, then interleave types so copies stay dispersed ----
    lines = build_line_copies(line_seq, line_fam, a.n_line, a.line_div, rng)
    sines = build_sine_copies(sine_seq, sine_fam, line_sine_frag, a.n_sine, a.sine_div, rng)
    ltr_host = build_ltr_host(ltr, intr, line5utr, ltr_d_family(ltr_d), rng)
    hel_cargo = build_one_sine(sine_seq, sine_fam, a.sine_div, rng)
    hel_host = build_helitron_host(hel_seq, hel_cargo, "rice_Helitron_Os0906", rng)

    elements = interleave(lines, sines, [ltr_host, hel_host], rng)

    contig_seq, rows = assemble(elements, rng)

    os.makedirs(a.out, exist_ok=True)
    write_fasta(os.path.join(a.out, a.contig + ".fa"), a.contig, contig_seq)
    write_manifest(os.path.join(a.out, "manifest.tsv"), a.contig, rows)

    print(f"contig {a.contig}: {len(contig_seq):,} bp")
    print(f"elements: {a.n_line} LINE, {a.n_sine} SINE, 1 LTR(+LINE frag), "
          f"1 Helitron(+SINE cargo)")
    print(f"LINE master={line_fam} ({len(line_seq)} bp)  SINE master={sine_fam} "
          f"({len(sine_seq)} bp)")
    print(f"manifest rows: {len(rows)}  -> {os.path.join(a.out, 'manifest.tsv')}")


def ltr_d_family(ltr_d):
    # family label from the _INT header base id (best-effort)
    for k in ltr_d:
        if "_INT" in k:
            return "rice_Copia_" + k.split("_")[0]
    return "rice_Copia"


def interleave(lines, sines, hosts, rng):
    """Round-robin lines+sines (keeps same-family copies dispersed), then insert each
    host at an evenly-spread position so no two engineered hosts are adjacent. Adjacent
    LTR+Helitron hosts make HelitronScanner detect a single chimeric element spanning
    both, which TEsorter then reclassifies as LTR and drops from the Helitron lib."""
    base = []
    L, S = list(lines), list(sines)
    while L or S:
        if L:
            base.append(L.pop(0))
        if S:
            base.append(S.pop(0))
    out = list(base)
    k = len(hosts)
    n = len(out)
    for j, h in enumerate(hosts):
        pos = (j + 1) * n // (k + 1) + j   # evenly spread; +j accounts for prior inserts
        out.insert(pos, h)
    return out


if __name__ == "__main__":
    main()
