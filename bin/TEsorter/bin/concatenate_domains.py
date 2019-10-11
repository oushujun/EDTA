#!/bin/env python
import os
import sys
import re
import uuid
from Bio import SeqIO
from RunCmdsMP import run_cmd

def main(inSeq, domains, outSeq=sys.stdout, tmpdir='/tmp'):
	d_domain = {domain: [] for domain in domains}
	uid = uuid.uuid1()
	# intersect
	for rc in SeqIO.parse(inSeq, 'fasta'):
		domain = re.compile(r'gene=([^;\s]+)').search(rc.description).groups()[0]
		if domain in d_domain:
			raw_id = '#'.join(rc.id.split('#')[:-1])
			d_domain[domain] += [raw_id]
	
	i = 0
	for domain, rawids in d_domain.iteritems():
		i += 1
		if i == 1:
			intersect = set(rawids)
			continue
		intersect = intersect & set(rawids)
	print >>sys.stderr, '{} sequences contain {} domains'.format(len(intersect), domains)
	# open files
	d_file = {}
	files = []
	for i, domain in enumerate(domains):
		outfile = '{}/{}.{}.fa'.format(tmpdir, uid, i)
		fout = open(outfile, 'w')
		d_file[domain] = fout
		files += [outfile]

	# write
	for  rc in SeqIO.parse(inSeq, 'fasta'):
		domain = re.compile(r'gene=([^;\s]+)').search(rc.description).groups()[0]
		if not domain in d_file:
			continue
		raw_id = '#'.join(rc.id.split('#')[:-1])
		if raw_id not in intersect:
			continue
		fout = d_file[domain]
		rc.id = raw_id
		SeqIO.write(rc, fout, 'fasta')

	# close files
	for fout in d_file.values():
		fout.close()

	# align
	alnfiles = []
	for seqfile in files:
		alnfile = seqfile + '.aln'
		cmd = 'mafft --auto {} > {}'.format(seqfile, alnfile)
#		os.system(cmd)
		run_cmd(cmd, log=True)
		alnfiles += [alnfile]
	# concatenate
	catAln(alnfiles, outSeq)

def catAln(inALNs, outALN):
    d_seqs = {}
    lens = []
    for inALN in inALNs:
        for rc in SeqIO.parse(inALN, 'fasta'):
            sp = rc.id
            seq = str(rc.seq)
            try: d_seqs[sp] += [seq]
            except KeyError: d_seqs[sp] = [seq]
        lens += [len(seq)]
    description = 'genes:{} sites:{} blocks:{}'.format(len(lens), sum(lens), lens)
    for sp, seqs in d_seqs.items():
        seqs = ''.join(seqs)
        print >> outALN, '>{} {}\n{}'.format(sp, description, seqs)

if __name__ == '__main__':
	inSeq = sys.argv[1]
	domains = sys.argv[2:]
	main(inSeq=inSeq, domains=domains, outSeq=sys.stdout)
