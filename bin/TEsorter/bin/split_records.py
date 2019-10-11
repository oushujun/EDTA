#!/bin/env python
import sys
import os
import shutil
import argparse
from Bio import SeqIO

from small_tools import open_file as open

__version__ = '1.0'
def split_sam_by_chunk_size(inSam, prefix, chunk_size):
	if not isinstance(inSam, file):
		inSam = open(inSam)
	i = 0
	j = 0
	outfiles = []
	header = ''
	for line in inSam:
		if line.startswith('@'):
			header += line
			continue
		chunk_id = i/chunk_size + 1
		hname = 'f%s' % (chunk_id,)
		if hname in dir():
			pass
		else:
			last_hname = 'f%s' % (chunk_id-1,)
			if last_hname in dir():
				exec '%s.close()' % (last_hname, )
			outfile = '%s.%s.sam' % (prefix, chunk_id)
			outfiles += [outfile]
			exec '%s = open("%s", "w")' % (hname, outfile)
			exec '%s.write(header)' % (hname, )
			j += 1
			
		exec '%s.write(line)' % (hname, )
		i += 1
	# close files
	exec '%s.close()' % (hname, )
	return (i, j, chunk_size, outfiles)
	
def split_sam_by_chr(inSam, prefix):
	if not isinstance(inSam, file):
		inSam = open(inSam)
	header = ''
	d_chr = {}
	i = 0
	j = 0
	outfiles = []
	for line in inSam:
		if line.startswith('@'):
			header += line
			continue
		j += 1
		temp = line.split('\t')
		chrom = temp[2]
		if chrom in d_chr:
			hname = 'f%s' % (d_chr[chrom],)
		else:
			d_chr[chrom] = i + 1
			outfile = '%s.%s.sam' % (prefix, d_chr[chrom])
			outfiles += [outfile]
			hname = 'f%s' % (d_chr[chrom],)
			exec '%s = open("%s", "w")' % (hname, outfile)
			exec '%s.write(header)' % (hname, )
			i += 1
		exec '%s.write(line)' % (hname, )
	# close files
	for chrom in d_chr.keys():
		hname = 'f%s' % (d_chr[chrom],)
		exec '%s.close()' % (hname, )
	return (j, i, j/i, outfiles)

def split_sam_by_chunk_num(inSam, prefix, chunk_num):
	if not isinstance(inSam, file):
		inSam = open(inSam)
	# open files
	outfiles = []
	for chunk_id in range(chunk_num):
		chunk_id += 1
		outfile = '%s.%s.sam' % (prefix, chunk_id)
		outfiles += [outfile]
		hname = 'f%s' % (chunk_id,)
		exec '%s = open("%s", "w")' % (hname, outfile)
		
	i = 0
	header = ''
	for line in inSam:
		if line.startswith('@'):
			header += line
			continue
		if i == 0:
			for chunk_id in range(chunk_num):
				chunk_id += 1
				hname = 'f%s' % (chunk_id,)
				exec '%s.write(header)' % (hname, )
		chunk_id = i % chunk_num + 1
		hname = 'f%s' % (chunk_id,)
		exec '%s.write(line)' % (hname, )
		i += 1
	# close files
	for chunk_id in range(chunk_num):
		chunk_id += 1
		exec 'f%s.close()' % (chunk_id, )
	return (i, chunk_num, i/chunk_num, outfiles)

def split_paf_by_chr(inPaf, prefix, suffix=''):	
	if not isinstance(inPaf, file):
		inPaf = open(inPaf)
	d_chr = {}
	i = 0
	j = 0
	outfiles = []
	for line in inPaf:
		j += 1
		temp = line.split('\t')
		chrom = temp[5]
		if chrom in d_chr:
			hname = 'f%s' % (d_chr[chrom],)
		else:
			d_chr[chrom] = i + 1
			outfile = '%s.%s.paf%s' % (prefix, chrom, suffix)
			outfiles += [outfile]
			hname = 'f%s' % (d_chr[chrom],)
			exec '%s = open("%s", "w")' % (hname, outfile)
			i += 1
		exec '%s.write(line)' % (hname, )
	# close files
	for chrom in d_chr.keys():
		hname = 'f%s' % (d_chr[chrom],)
		exec '%s.close()' % (hname, )
	return (j, i, j/i, outfiles)
	
def split_fastx_by_chr(inFastx, prefix, seqfmt, suffix=''):
	if not isinstance(inFastx, file):
		inFastx = open(inFastx)
	i = 0
	j = 0
	outfiles = []
	for rc in parse_fastx(inFastx):
		j += 1
		last_hname = 'f%s' % (i-1,)
		if last_hname in dir():
			exec '%s.close()' % (last_hname, )
		chrom = rc.split()[0][1:]
		outfile = '%s.%s.%s%s' % (prefix, chrom, seqfmt, suffix)
		outfiles += [outfile]
		hname = 'f%s' % (i,)
		exec '%s = open("%s", "w")' % (hname, outfile)
		i += 1
		exec '%s.write(rc)' % (hname, )
	# close files
	exec '%s.close()' % (hname, )	
	return (j, i, j/i, outfiles)

def split_fastx_by_chunk_size(inFastx, prefix, chunk_size, seqfmt, suffix):
	if not isinstance(inFastx, file):
		inFastx = open(inFastx)
	i = 0
	j = 0
	outfiles = []
	header = ''
	for rc in parse_fastx(inFastx):
		chunk_id = i/chunk_size + 1
		hname = 'f%s' % (chunk_id,)
		if hname in dir():
			pass
		else:
			last_hname = 'f%s' % (chunk_id-1,)
			if last_hname in dir():
				exec '%s.close()' % (last_hname, )
			outfile = '%s.%s.%s%s' % (prefix, chunk_id, seqfmt, suffix)
			outfiles += [outfile]
			exec '%s = open("%s", "w")' % (hname, outfile)
			j += 1
		exec '%s.write(rc)' % (hname, )
		i += 1
	# close files
	exec '%s.close()' % (hname, )
	return (i, j, chunk_size, outfiles)

def split_fastx_by_chunk_num(inFastx, prefix, chunk_num, seqfmt, suffix):
	if not isinstance(inFastx, file):
		inFastx = open(inFastx)
	# open files
	outfiles = []
	for chunk_id in range(chunk_num):
		chunk_id += 1
		outfile = '%s.%s.%s%s' % (prefix, chunk_id, seqfmt, suffix)
		outfiles += [outfile]
		hname = 'f%s' % (chunk_id,)
		exec '%s = open("%s", "w")' % (hname, outfile)
		
	i = 0
	for rc in parse_fastx(inFastx):
		chunk_id = i % chunk_num + 1
		hname = 'f%s' % (chunk_id,)
		exec '%s.write(rc)' % (hname, )
		i += 1
	# close files
	for chunk_id in range(chunk_num):
		chunk_id += 1
		exec 'f%s.close()' % (chunk_id, )
	return (i, chunk_num, i/chunk_num, outfiles)

def split_fastx_by_size(inFastx, prefix, chunk_num, seqfmt, suffix, out_random=True):
	import binpacking
	d_seq = {}
	d_len = {}
	for rc in SeqIO.parse(inFastx, seqfmt):
		d_seq[rc.id] = rc
		d_len[rc.id] = len(rc.seq)
	bins = binpacking.to_constant_bin_number(d_len, chunk_num)
	i = 0
	j = 0
	outfiles = []
	if out_random:
		import random
		random.shuffle(bins)
	for d_bin in bins:
		chunk_id = i + 1
		out_file = '%s.%s.%s%s' % (prefix, chunk_id, seqfmt, suffix)
		outfiles += [out_file]
		f = open(out_file,'w')
		for id in d_bin.keys():
			j +=1
			SeqIO.write(d_seq[id], f, seqfmt)
		f.close()
		i += 1
	return (j, chunk_num, j/chunk_num, outfiles)

def parse_fastx(inFastx):
	i = 0
	lines = []
	for line in inFastx:
		i += 1
		if i == 1:
			if line[0] == '>':
				seqfmt = 'fasta'
			elif line[0] == '@':
				seqfmt = 'fastq'
			else:
				raise ValueError('Unknown sequence format, neither fasta nor fastq')
			lines.append(line)
			continue
		if (seqfmt == 'fasta' and line[0] == '>') or (seqfmt == 'fastq' and i % 4 == 1):
			yield ''.join(lines)
			lines = []
		lines.append(line)
	yield ''.join(lines)
	
def main():
	parser = argparse.ArgumentParser(version=__version__)
	parser.add_argument("-i","--input", action="store",type=str,
					dest="input", default=sys.stdin, 
					help="input [default=%(default)s]")
	parser.add_argument("--prefix", action="store",
					dest="prefix", default='chunk', 
					help="output prefix [default=%(default)s]")
	parser.add_argument("-n","--chunk-number", action="store",type=int,
					dest="chunk_num", default=None, 
					help="number of chunk [default=%(default)s]")
	parser.add_argument("-s","--chunk-size", action="store", type=int,
					dest="chunk_size", default=None, 
					help='size of chunk [default=%(default)s]')
	parser.add_argument("-f","--format", action="store", 
					dest="rcfmt", default='fasta', 
					choices=['fasta', 'fastq', 'fastx', 'sam', 'paf'], 
					help="record file format [default=%(default)s]")
	parser.add_argument("--gzip-output", action="store_true",
					dest="gzip_output", default=False, 
					help="if gzip output [default=%(default)s]")
	parser.add_argument("--by-size", action="store_true",
					dest="by_size", default=False, 
					help='split by size (binpacking, only for fastx format) [default=%(default)s]')
	parser.add_argument("--by-chrom", action="store_true",
					dest="by_chr", default=False,
					help='split by mapped chromsome [default=%(default)s]')
	parser.add_argument("-pfn", "--print-filenames", action="store_true",
					dest="print_filenames", default=False, 
					help='print filenames [default=%(default)s]')

	#parser.print_help()
	options = parser.parse_args()
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit()
	if options.chunk_num and options.chunk_size:
		parser.print_help(sys.stderr)
		print >> sys.stderr, 'chunk-number and chunk-size is not compatible.'
		sys.exit()
	elif not (options.chunk_num or options.chunk_size) and not options.by_chr:
		parser.print_help(sys.stderr)
		print >> sys.stderr, 'either chunk-number or chunk-size must be speicfied.'
		sys.exit()
	if options.gzip_output:
		suffix = '.gz'
		print >>sys.err, 'Warning: for large dataset, gzip output is very slow. You may want to diable it.'
	else:
		suffix = ''
	if not isinstance(options.input, file):
		options.input = open(options.input)
	
	#execute
	if options.rcfmt == 'sam':
		if options.chunk_num:
			stats = split_sam_by_chunk_num(options.input, options.prefix, options.chunk_num)
		elif options.chunk_size:
			stats = split_sam_by_chunk_size(options.input, options.prefix, options.chunk_size)
		elif options.by_chr:
			stats = split_sam_by_chr(options.input, options.prefix,)
	elif options.rcfmt in set(['fasta', 'fastq', 'fastx']):
		if options.chunk_num:
			if options.by_size:
				stats = split_fastx_by_size(options.input, options.prefix, 
				options.chunk_num, options.rcfmt, suffix)
			else:
				stats = split_fastx_by_chunk_num(options.input, options.prefix, 
				options.chunk_num, options.rcfmt, suffix)
		elif options.chunk_size:
			stats = split_fastx_by_chunk_size(options.input, options.prefix, 
			options.chunk_size, options.rcfmt, suffix)
		elif options.by_chr:
			stats = split_fastx_by_chr(options.input, options.prefix, options.rcfmt, suffix)
	elif options.rcfmt == 'paf':
		if options.by_chr:
			stats = split_paf_by_chr(options.input, options.prefix, suffix)
	n_records, n_chunks, per_chunk, outfiles = stats
	print >> sys.stderr, 'total %s records, splited into %s chunks, %s per chunk' % (n_records, n_chunks, per_chunk)
	if options.print_filenames:
		print >> sys.stdout, '\n'.join(outfiles)

if __name__ == '__main__':
	main()
