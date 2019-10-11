#!/bin/env python
# coding: utf-8
'''Author: Zhang, Ren-Gang and Wang, Zhao-Xuan
'''
import sys
import os
import re
import shutil
import glob
import argparse
from collections import Counter, OrderedDict
from Bio import SeqIO
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s -%(levelname)s- %(message)s')
logger = logging.getLogger(__name__)

bindir = os.path.dirname(os.path.realpath(__file__))
sys.path = [bindir + '/bin'] + sys.path

from translate_seq import six_frame_translate
# for multi-processing HMMScan
from RunCmdsMP import run_cmd, pp_run
from split_records import split_fastx_by_chunk_num
# for pass-2 blast classifying
from get_record import get_records

__version__ = '1.0'

DB = {
	'gydb' : bindir + '/database/GyDB2.hmm',
	'rexdb': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm',
	'rexdb-plant': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0.hmm',
	'rexdb-metazoa': bindir + '/database/REXdb_protein_database_metazoa_v3.hmm',
	}
BLASType = {
    'qseqid': str,
    'sseqid': str,
    'pident': float,
    'length': int,
    'mismatch': int,
    'gapopen': int,
    'qstart': int,
    'qend': int,
    'sstart': int,
    'send': int,
    'evalue': float,
    'bitscore': float,
    'qlen': int,
    'slen': int,
    'qcovs': float,
    'qcovhsp': float,
    'sstrand': str,
	}
	
	
def Args():
	parser = argparse.ArgumentParser(version=__version__)
	parser.add_argument("sequence", action="store",type=str,
					help="input TE sequences in fasta format [required]")
	parser.add_argument("-db","--hmm-database", action="store",type=str,
					default='rexdb', choices=DB.keys(),  
					help="the database used [default=%(default)s]")
	parser.add_argument("-st","--seq-type", action="store",type=str,
					default='nucl', choices=['nucl', 'prot'],  
					help="'nucl' for DNA or 'prot' for protein [default=%(default)s]")
	parser.add_argument("-pre", "--prefix", action="store",
					default=None, type=str,
					help="output prefix [default='{-s}.{-db}']")
	parser.add_argument("-fw", "--force-write-hmmscan", action="store_true",
					default=False, 
					help="if False, will use the existed hmmscan outfile and skip hmmscan [default=%(default)s]")
	parser.add_argument("-p", "--processors", action="store",
					default=4, type=int,
					help="processors to use [default=%(default)s]")
	parser.add_argument("-tmp", "--tmp-dir", action="store",
					default='./tmp', type=str,
					help="directory for temporary files [default=%(default)s]")
	parser.add_argument("-cov", "--min-coverage", action="store",
					default=20, type=float,
					help="mininum coverage for protein domains in HMMScan output [default=%(default)s]")
	parser.add_argument("-eval", "--max-evalue", action="store",
					default=1e-3, type=float,
					help="maxinum E-value for protein domains in HMMScan output [default=%(default)s]")				
	parser.add_argument("-dp2", "--disable-pass2", action="store_true",
					default=False, 
					help="do not further classify the unclassified sequences [default=%(default)s for `nucl`, True for `prot`]")
	parser.add_argument("-rule", "--pass2-rule", action="store",
					default='80-80-80', type=str,
					help="classifying rule [identity-coverage-length] in pass-2 based on similarity [default=%(default)s]")
	parser.add_argument("-nolib", "--no-library", action="store_true",
					default=False, 
					help="do not generate a library file for RepeatMasker [default=%(default)s]")
	parser.add_argument("-norc", "--no-reverse", action="store_true",
					default=False,
					help="do not reverse complement sequences if they are detected in minus strand [default=%(default)s]")
	parser.add_argument("-nocln", "--no-cleanup", action="store_true",
					default=False,
					help="do not clean up the temporary directory [default=%(default)s]")
	args = parser.parse_args()
	if args.prefix is None:
		args.prefix = '{}.{}'.format(os.path.basename(args.sequence), args.hmm_database)
	
	if args.seq_type == 'prot':
		args.disable_pass2 = True
		args.no_reverse = True
	if not args.disable_pass2:
		for key, par in zip(['p2_identity', 'p2_coverage', 'p2_length'], args.pass2_rule.split('-')):
			setattr(args, key, float(par))
	return args

def pipeline(args):
	logger.info( 'VARS: {}'.format(vars(args)))
	logger.info( 'checking dependencies:' )
	Dependency().check_hmmer(db=DB[args.hmm_database])
	if not args.disable_pass2:
		Dependency().check_blast()
	if not os.path.exists(args.tmp_dir):
		os.makedirs(args.tmp_dir)
	logger.info( 'Start classifying pipeline' )
	seq_num = len([1 for rc in SeqIO.parse(args.sequence, 'fasta')])
	logger.info('total {} sequences'.format(seq_num))	
	# search against DB and parse
	gff, geneSeq = LTRlibAnn(
			ltrlib = args.sequence, 
			hmmdb = args.hmm_database, 
			seqtype = args.seq_type,
			prefix = args.prefix,
			force_write_hmmscan = args.force_write_hmmscan, 
			processors = args.processors,
			tmpdir = args.tmp_dir,
			mincov = args.min_coverage,
			maxeval = args.max_evalue,
			)
			
	# classify	
	classify_out = args.prefix + '.cls.tsv'
	fc = open(classify_out, 'w')
	d_class = OrderedDict()
	for rc in Classifier(gff, db=args.hmm_database, fout=fc):
		d_class[rc.id] = rc
	fc.close()
	classfied_num = len(d_class)
	logger.info('{} sequences classified by HMM'.format(classfied_num))
	logger.info('see protein domain sequences in `{}` and annotation gff3 file in `{}`'.format(geneSeq, gff))

	# pass-2 classify
	if classfied_num == 0 and not args.disable_pass2:
			logger.warn('skipping pass-2 classification for zero classification in step-1')
			args.disable_pass2 = True
	if not args.disable_pass2:
		logger.info('classifying the unclassified sequences by searching against the classified ones')
		classified_seq = '{}/pass1_classified.fa'.format(args.tmp_dir)
		unclassified_seq = '{}/pass1_unclassified.fa'.format(args.tmp_dir)
		get_records(args.sequence, classified_seq, d_class.keys(), type='fasta', process='get')
		get_records(args.sequence, unclassified_seq, d_class.keys(), type='fasta', process='remove')

		logger.info('using the {} rule'.format(args.pass2_rule))
		d_class2 = classify_by_blast(classified_seq, unclassified_seq, 
						seqtype=args.seq_type, ncpu=args.processors,
						min_identtity=args.p2_identity, min_coverge=args.p2_coverage, min_length=args.p2_length,
						)
		fc = open(classify_out, 'a')
		for unclfed_id, clfed_id in d_class2.items():
			clfed = d_class[clfed_id]
			order, superfamily, clade = clfed.order, clfed.superfamily, 'unknown'
			line = [unclfed_id, order, superfamily, clade, 'none', '?', 'none']
			print >> fc, '\t'.join(line)
			# update
			d_class[unclfed_id] = CommonClassification(*line)
		fc.close()
		logger.info('{} sequences classified in pass 2'.format(len(d_class2)))
		logger.info('total {} sequences classified.'.format(len(d_class)))
	logger.info('see classified sequences in `{}`'.format(classify_out))

	# output library
	if not args.no_library:
		out_lib = args.prefix + '.cls.lib'
		logger.info( 'writing library for RepeatMasker in `{}`'.format(out_lib) )
		fout = open(out_lib, 'w')
		for rc in SeqIO.parse(args.sequence, 'fasta'):
			if rc.id in d_class:
				cl = d_class[rc.id]
				strand = cl.strand
				cl = fmt_cls(cl.order, cl.superfamily, cl.clade)
				if not args.no_reverse and strand == '-':
					rc.seq = rc.seq.reverse_complement()
			else:
				cl = 'Unknown'
			rc.id = rc.id.split('#')[0] + '#' + cl
			SeqIO.write(rc, fout, 'fasta')
		fout.close()

	# rename id of protein domains
	pep_lib = args.prefix + '.cls.pep'
	logger.info( 'writing classified protein domains in `{}`'.format(pep_lib) )
	fout = open(pep_lib, 'w')
	for rc in SeqIO.parse(geneSeq, 'fasta'):
		raw_id = '|'.join(rc.id.split('|')[:-1])
		assert raw_id in d_class
		cl = d_class[raw_id]
		cl = fmt_cls(cl.order, cl.superfamily, cl.clade)
		d_desc = dict([pair.split('=')for pair in rc.description.split()[-1].split(';')])
		gene, clade = d_desc['gene'], d_desc['clade']
		new_id = '{}#{}#{}|{}'.format(raw_id.split('#')[0], cl, gene, clade)
		rc.id = new_id
		SeqIO.write(rc, fout, 'fasta')
	fout.close()
	logger.info('Summary of classifications:')
	summary(d_class)
	logger.info( 'Pipeline done.' )

	# clean up
	if not args.no_cleanup:
		logger.info( 'cleaning the temporary directory {}'.format(args.tmp_dir) )
		shutil.rmtree(args.tmp_dir)
def summary(d_class):
	d_sum = {}
	for sid, clf in d_class.iteritems():
		key = (clf.order, clf.superfamily)
		d_sum[key] = [0, 0, [], 0] # #seqs, #seqs in clades, #clades, #full domains
	for sid, clf in d_class.iteritems():
		key = (clf.order, clf.superfamily)
		d_sum[key][0] += 1
		if clf.clade not in {'unknown', 'mixture'}:
			d_sum[key][1] += 1
			d_sum[key][2] += [clf.clade]
		if clf.completed == 'yes':
			d_sum[key][3] += 1
	out_order = ['LTR', 'pararetrovirus', 'DIRS', 'Penelope', 'LINE', 'TIR', 'Helitron', 'Maverick', 'mixture', 'Unknown']
	template = '{:<16}{:<16}{:>15}{:>15}{:>15}{:>15}'
	line = ['Order', 'Superfamily', '# of Sequences', '# of Clade Sequences', '# of Clades', '# of full Domains']
	
	print >> sys.stderr, template.format(*line)
	for (order, superfamliy), summary in \
			sorted(d_sum.items(), key=lambda x: (out_order.index(x[0][0]), x[0][1])):
		line = [order, superfamliy, summary[0], summary[1], len(set(summary[2])), summary[3]]
		line = map(str, line)
		line = template.format(*line)
		print >> sys.stderr, line
def fmt_cls(*args):
	values = []
	for arg in args:
		if arg == 'unknown' or arg in set(values):
			continue
		values += [arg]
	return '/'.join(values)
	
class CommonClassification(object):
	def __init__(self, id=None, order=None, superfamily=None, 
				clade=None, completed=None, strand=None, domains=None):
		self.id = id
		self.order = order
		self.superfamily = superfamily
		self.clade = clade
		self.completed = completed
		self.strand = strand
		self.domains = domains

def classify_by_blast(db_seq, qry_seq, blast_out=None, seqtype='nucl', ncpu=4, 
					  min_identtity=80, min_coverge=80, min_length=80):
	'''pass-2 classify'''
	blast_outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sstrand'"
	blast_out = blast(db_seq, qry_seq, seqtype=seqtype, blast_out=blast_out, blast_outfmt=blast_outfmt, ncpu=ncpu)
	with open(blast_out+'.best', 'w') as fb:
		d_best_hit = BlastOut(blast_out, blast_outfmt).filter_besthit(fout=fb)
	for qseqid, rc in d_best_hit.iteritems():
		if not (rc.pident >= min_identtity and rc.qcovs >= min_coverge and rc.length >= min_length):
			del d_best_hit[qseqid]
	d_class = OrderedDict([(qseqid, rc.sseqid) for qseqid, rc in d_best_hit.iteritems()])
	return d_class
	
def blast(db_seq, qry_seq, seqtype='nucl', blast_out=None, blast_outfmt=None, ncpu=4):
	if seqtype == 'nucl':
		blast_app = 'blastn'
	elif seqtype == 'prot':
		blast_app = 'blastp'
	else:
		raise ValueError('Unknown molecule type "{}" for blast'.format(seqtype))
	if blast_out is None:
		blast_out = qry_seq + '.blastout'
	if blast_outfmt is None:
		blast_outfmt = '6'
	cmd = 'makeblastdb -in {} -dbtype {}'.format(db_seq, seqtype)
	run_cmd(cmd, logger=logger)
		
	cmd = '{} -query {} -db {} -out {} -outfmt {} -num_threads {}'.format(blast_app, qry_seq, db_seq, blast_out, blast_outfmt, ncpu)
#	cmd += " " + blast_opts
	run_cmd(cmd, logger=logger)
	return blast_out

class BlastOut(object):
	def __init__(self, blast_out, outfmt=None):
		self.blast_out = blast_out
		if outfmt is not None:
			outfmt = outfmt.strip(''''"''')
		if outfmt is None:
			self.outfmt = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()
		elif outfmt[0] == '6':
			self.outfmt = outfmt.split()[1:]
		else:
			raise ValueError('Only support for blast outfmt 6 = tabular, but {} input'.format(outfmt))
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.blast_out):
			values = line.strip().split('\t')
			yield BlastOutRecord(self.outfmt, values)
	def filter_besthit(self, fout=sys.stdout):
		d_best_hit = OrderedDict()
		for rc in self.parse():
			if rc.qseqid in d_best_hit:
				if rc.bitscore > d_best_hit[rc.qseqid]:
					d_best_hit[rc.qseqid] = rc
			else:
				d_best_hit[rc.qseqid] = rc
		if fout is None:
			return d_best_hit
		for qseqid, rc in d_best_hit.iteritems():
			rc.write(fout)
		return d_best_hit
		
class BlastOutRecord(object):
	def __init__(self, outfmt, values):
		self.values = values
		for key, value in zip(outfmt, values):
			setattr(self, key, BLASType[key](value))
	def write(self, fout=sys.stdout):
		print >> fout, '\t'.join(self.values)
		
class Classifier(object):
	def __init__(self, gff, db='rexdb', fout=sys.stdout): # gff is sorted
		self.gff = gff
		self.db = db
		self.fout = fout
		if self.db.startswith('rexdb'):
			self.markers = {'GAG', 'PROT', 'INT', 'RT', 'RH'}
		elif self.db == 'gydb':
			self.markers = {'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV'}
	def __iter__(self):
		return self.classify()
	def parse(self):
		record = []
		last_lid = ''
		for line in open(self.gff):
			if line.startswith('#'):
				continue
			line = LTRgffLine(line)
			lid = line.ltrid
			if record and  not lid == last_lid:
				yield record
				record = []
			record.append(line)
			last_lid = lid
		yield record
	def classify(self, ):
		line = ['#TE', 'Order', 'Superfamily', 'Clade', 'Complete', 'Strand', 'Domains']
		print >> self.fout, '\t'.join(line)
		for rc in self.parse():
			rc_flt = rc
			strands = [line.strand for line in rc_flt]
			strands_num = len(set(strands))
			if strands_num >1 :
				strand = '?'
			elif strands_num == 1:
				strand = strands[0]
			else:
				continue
			if strand == '-':
				rc_flt.reverse()
			lid = rc_flt[0].ltrid
			domains = ' '.join(['{}|{}'.format(line.gene, line.clade)  for line in rc])
			genes  = [line.gene  for line in rc_flt]
			clades = [line.clade for line in rc_flt]
			names = [line.name for line in rc_flt]
			if self.db.startswith('rexdb'):
				order, superfamily, max_clade, coding = self.identify_rexdb(genes, names)
			elif self.db == 'gydb':
				order, superfamily, max_clade, coding = self.identify_gydb(genes, clades)
			line = [lid, order, superfamily, max_clade, coding, strand, domains]
			print >> self.fout, '\t'.join(line)
			yield CommonClassification(*line)
	def identify_rexdb(self, genes, clades):
		perfect_structure = {
#            ('LTR', 'Copia'): ['Ty1-GAG', 'Ty1-PROT', 'Ty1-INT', 'Ty1-RT', 'Ty1-RH'],
#            ('LTR', 'Gypsy'): ['Ty3-GAG', 'Ty3-PROT', 'Ty3-RT', 'Ty3-RH', 'Ty3-INT'],
			('LTR', 'Copia'): ['GAG', 'PROT', 'INT', 'RT', 'RH'],
			('LTR', 'Gypsy'): ['GAG', 'PROT', 'RT', 'RH', 'INT'],
			('LTR', 'Bel-Pao'): ['GAG', 'PROT', 'RT', 'RH', 'INT'],
			}
		clade_count = Counter(clades)
		max_clade = max(clade_count, key=lambda x: clade_count[x])
		order, superfamily = self._parse_rexdb(max_clade)
		if len(clade_count) == 1 or clade_count[max_clade] > 1:
			max_clade = max_clade.split('/')[-1]
		elif len(clade_count) > 1:
			max_clade = 'mixture'
			superfamlies = [self._parse_rexdb(clade)[1] for clade in clades]
			if len(Counter(superfamlies)) > 1:
				superfamily = 'mixture'
				orders = [self._parse_rexdb(clade)[0] for clade in clades]
				if len(Counter(orders)) > 1:
					order = 'mixture'
		try:
			ordered_genes = perfect_structure[(order, superfamily)]
			my_genes = [gene for gene in genes if gene in set(ordered_genes)]
			if ordered_genes == my_genes:
				coding = 'yes' # completed gene structure
			else:
				coding = 'no'
		except KeyError:
			coding = 'unknown'
		if superfamily not in {'Copia', 'Gypsy'}:
			max_clade = 'unknown'
		if max_clade.startswith('Ty'): # Ty3_gypsy, Ty1_copia, Ty1-outgroup in metazoa_v3
			max_clade = 'unknown'
		return order, superfamily, max_clade, coding
	def _parse_rexdb(self, clade): # full clade name
		if clade.startswith('Class_I/LTR/Ty1_copia'):
			order, superfamily = 'LTR', 'Copia'
		elif clade.startswith('Class_I/LTR/Ty3_gypsy'):
			order, superfamily = 'LTR', 'Gypsy'
		elif clade.startswith('Class_I/LTR/'): # LTR/Bel-Pao, LTR/Retrovirus
			order, superfamily = clade.split('/')[1:3]
		elif clade.startswith('Class_I/'): # LINE, pararetrovirus, Penelope, DIRS
			order, superfamily = clade.split('/')[1], 'unknown'
		elif clade.startswith('Class_II/'): # TIR/hAT, Helitro, Maverick
			try: order, superfamily = clade.split('/')[2:4]
			except ValueError: order, superfamily = clade.split('/')[2], 'unknown'
		elif clade.startswith('NA'): # "NA:Retrovirus-RH"
			order, superfamily = 'LTR', 'Retrovirus'
		else:
			logger.error('Unknown clade {}'.format(clade))
		return order, superfamily
	def identify_gydb(self, genes, clades):
		perfect_structure = {
			('LTR', 'Copia')         : ['GAG', 'AP', 'INT', 'RT', 'RNaseH'],
			('LTR', 'Gypsy')         : ['GAG', 'AP', 'RT', 'RNaseH', 'INT'],
			('LTR', 'Pao')           : ['GAG', 'AP', 'RT', 'RNaseH', 'INT'],
			('LTR', 'Retroviridae')  : ['GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ENV'],
			('LTR', 'Caulimoviridae'): ['GAG', 'AP', 'RT', 'RNaseH', ],
			}
		d_map = self.clade_map
		clade_count = Counter(clades)
		max_clade = max(clade_count, key=lambda x: clade_count[x])
		try: (order, superfamily) = d_map[max_clade]
		except KeyError: 
			(order, superfamily) = ('Unknown', 'unknown')
			logger.warn( 'unknown clade: {}'.format(max_clade) )
		try:
			ordered_genes = perfect_structure[(order, superfamily)]
			my_genes = [gene for gene in genes if gene in set(ordered_genes)]
			if ordered_genes == my_genes:
				coding = 'yes' # completed gene structure and the same order
			else:
				coding = 'no'
		except KeyError:
			coding = 'unknown'
		return order, superfamily, max_clade, coding
	def replace_annotation(self, rawseq, fout=sys.stdout, idmap = None):
		d_class = self.classification
		i = 0
		for rc in SeqIO.parse(rawseq, 'fasta'):
			if idmap is None:
				intid = rc.id.split('#')[0]
			else:
				try: intid = idmap[rc.id.split('#')[0]]
				except KeyError as e:	# this should be rare
					logger.warn( 'skipped KeyError: {}'.format(e) )
			if intid in d_class:
				neword, newfam = d_class[intid]
				re_org = self.re_orgnize(rc.id, neword, newfam)
				if re_org:
					i += 1
					rc.id = re_org
			SeqIO.write(rc, fout, 'fasta')
		logger.info( 'sequences re-classified'.format(i) )
	def re_orgnize(self, rawid, neword, newfam):
		rawid, rawcls = rawid.split('#')
		try: raword, rawfam = rawcls.split('/')[:2]
		except ValueError: raword, rawfam = rawcls, 'unknown'
		if raword.lower() == 'unknown' or (rawfam.lower() == 'unknown' and raword == neword):
			return '{}#{}/{}'.format(rawid, neword, newfam)
		else:
			return False
	@property
	def classification(self):
		return {rc.ltrid.split('#')[0]: (rc.order, rc.superfamily) for rc in self.classify()}
	@property
	def clade_map(self):
		return {rc.clade: (rc.order, rc.superfamily) for rc in CladeInfo()}
			
class CladeInfo():
	def __init__(self, infile=DB['gydb']+'.info'):
		self.infile = infile
		self.clade_map = {
			'Ty_(Pseudovirus)': 'pseudovirus',
			'Cer2-3': 'cer2-3',
			'412/Mdg1': '412_mdg1',
			'TF1-2': 'TF',
			'Micropia/Mdg3': 'micropia_mdg3',
			'CoDi-I': 'codi_I',
			'CoDi-II': 'codi_II',
			'17.6': '17_6',
			}
	def __iter__(self):
		return self.parse()
	def parse(self):
		i = 0
		for line in open(self.infile):
			i += 1
			temp = line.strip().split('\t')
			if i == 1:
				title = temp
				continue
			self.dict = dict(zip(title, temp))
			if self.dict['Clade'] == 'NA':
				self.clade = self.dict['Cluster_or_genus']
			else:
				self.clade = self.dict['Clade']
			self.superfamily = self.dict['Family'].split('/')[-1]
			if self.superfamily == 'Retroviridae':	# deltaretroviridae gammaretroviridae
				self.clade = self.dict['Cluster_or_genus'].replace('virus', 'viridae')
			if self.superfamily == 'Retrovirus':	# an exception
				self.superfamily = 'Retroviridae'
			self.order = 'LTR' if self.dict['System'] in {'LTR_retroelements', 'LTR_Retroelements', 'LTR_retroid_elements'} else self.dict['System']
			yield self
			if self.clade in self.clade_map:
				self.clade = self.clade_map[self.clade]
				yield self
				if self.clade == '412_mdg1':
					self.clade = '412-mdg1'  # 412-mdg1 and 412_mdg1
					yield self

			self.clade = self.clade.replace('-', '_') # A-clade V-clade C-clade
			yield self
			self.clade = self.clade.lower()
			yield self
		self.order, self.superfamily, self.clade, self.dict = ['LTR', 'Copia', 'ty1/copia', {}]  # AP_ty1/copia
		yield self
		order_map = { 	# some unknown clade
			'retroelement': 'LTR',
			'retroviridae': 'LTR',
			'B-type_betaretroviridae': 'LTR',
			'D-type_betaretroviridae': 'LTR',
			'caulimoviruses': 'LTR',
			'caulimoviridae_dom2': 'LTR',
			'errantiviridae': 'LTR',
			'retropepsins': 'LTR',
			'VPX_retroviridae': 'LTR',
			'cog5550': 'Unknown',
			'ddi': 'Unknown',
			'dtg_ilg_template': 'Unknown',
			'saspase': 'Unknown',
			'GIN1': 'Unknown',
			'shadow': 'Unknown',
			'all': 'Unknown',
			}
		for clade, order in order_map.items():
			self.order, self.superfamily, self.clade, self.dict = [order, 'unknown', clade, {}]  # CHR_retroelement
			yield self	
				
class GffLine(object):
	def __init__(self, line):
		temp = line.strip().split('\t')
		self.chr, self.source, self.type, self.start, self.end, self.score, self.strand, self.frame, self.attributes = temp
		self.start, self.end = int(self.start), int(self.end)
		try: self.score = float(self.score)
		except: pass
		try: self.frame = int(self.frame)
		except: pass
		self.attributes = self.parse(self.attributes)
	def parse(self, attributes):
		return dict(kv.split('=') for kv in attributes.split(';'))
class LTRgffLine(GffLine):
	def __init__(self, line):
		super(LTRgffLine, self).__init__(line)
		self.gene = self.attributes['gene']
		self.clade = self.attributes['clade']
		self.ltrid = '|'.join(self.attributes['ID'].split('|')[:-1])
		self.name = self.attributes['ID'].split('|')[-1].split(':')[0]

class HmmScan(object):
	def __init__(self, hmmout, hmmfmt='domtbl'):
		self.hmmout = hmmout
		self.hmmfmt = hmmfmt
	def __iter__(self):
		return self.parse()
	def parse(self):
		for line in open(self.hmmout):
			if line.startswith('#'):
				continue
			if self.hmmfmt == 'domtbl':
				yield HmmDomRecord(line)
class HmmDomRecord(object):
	def __init__(self, line):
		temp = line.strip().split()
		self.tname, self.tacc, self.tlen, self.qname, self.qacc, self.qlen, \
			self.evalue, self.score, self.bias, \
			self.domi, self.domn, self.cevalue, self.ievalue, self.domscore, self.dombias, \
			self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend, \
			self.acc \
				= temp[:22]
		self.tlen, self.qlen, self.domi, self.domn, \
			self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend = \
			map(int, [self.tlen, self.qlen, self.domi, self.domn, \
				self.hmmstart, self.hmmend, self.alnstart, self.alnend, self.envstart, self.envend])
		self.evalue, self.score, self.bias, self.cevalue, self.ievalue, self.domscore, self.dombias, self.acc = \
			map(float, [self.evalue, self.score, self.bias, self.cevalue, self.ievalue, self.domscore, self.dombias, self.acc])
		self.tdesc = ' '.join(temp[22:])
	@property
	def hmmcov(self):
		return round(1e2*(self.hmmend - self.hmmstart + 1) / self.tlen, 1)
		
def seq2dict(inSeq):
	return dict([(rc.id, rc) for rc in SeqIO.parse(inSeq, 'fasta')])
	
def parse_hmmname(hmmname, db='gydb'):
	db = db.lower()
	if db == 'gydb':
		temp = hmmname.split('_')
		gene, clade = temp[0], '_'.join(temp[1:])
	elif db.startswith('rexdb'):	# Class_I/LTR/Ty3_gypsy/chromovirus/Tekay:Ty3-RT
		gene = hmmname.split(':')[1] #.split('-')[1]
		clade = hmmname.split(':')[0].split('/')[-1]
	elif db.startswith('pfam'):
		gene = hmmname
		clade = hmmname
	return gene, clade
class HmmCluster(object):
	def __int__(self, hmmout, seqtype = 'nucl'): # only for nucl
		self.hmmout = hmmout
		self.seqtype = seqtype
	def cluster(self):
		d_cluster = OrderedDict()
		for rc in HmmScan(self.hmmout):
			suffix = rc.qname.split('|')[-1]
			qid = '|'.join(rc.qname.split('|')[:-1])
			strand, frame = parse_frame(rc.qname.split('|')[-1])
			rc.qid, rc.strand, rc.frame = qid, strand, frame
			key = (rc.tname, qid, strand)
			try: d_cluster[key] += [rc]
			except KeyError: d_cluster[key] = [rc]

		return d_cluster
	def extand(self, records, max_mal=5):
		if len(records) == 1:
			return records
		records = sorted(records, key=lambda x:x.hmmstart)
		best_idx, best_rc = self.maxscore(records)
		new_rcs = [best_rc]
		for rc in records[:best_idx][::-1]:		# <<- left extand
			right_rc = new_rcs[0]
			mal_pos = right_rc.hmmstart - rc.hmmend
			if abs(mal_pos) <= max_mal:
				if mal_pos <= 0:
					diff = 1-mal_pos
					rc.hmmend -= diff
					rc.alnend -= diff
				new_rcs = [rc] + new_rcs
		for rc in records[best_idx:]:			# ->> right extand
			left_rc = new_rcs[-1]
			mal_pos = rc.hmmstart - left_rc.hmmend
			if abs(mal_pos) <= max_mal:
				if mal_pos <= 0:
					diff = 1-mal_pos
					rc.hmmstart += diff
					rc.alnstart += diff
				new_rcs += [rc]
		return HmmClusterRecord(new_rcs)
	def maxscore(self, records):
		for i, rc in enumerate(records):
			if i == 0:
				best = (i, rc)
				continue
			if rc.score > best[1].score:
				best = (i, rc)
		return best
			
class HmmClusterRecord(object):
	def __init__(self, records):
		self.records = records
		self.score = sum([rc.score for rc in records])
		self.hmmstart = records[0].hmmstart
		self.hmmend = records[-1].hmmend
		self.alnstart = records[0].alnstart
		self.alnend = records[-1].alnend
		self.tlen = records[0].tlen
		self.hmmcov = round(1e2*(self.hmmend - self.hmmstart + 1) / self.tlen, 1)
		self.evalue = multi(*[rc.evalue for rc in records])
def multi(*n):
	result = 1
	for i in n:
		result = result * i
	return result
	
def hmm2best(inSeq, inHmmouts, prefix=None, db='rexdb', seqtype='nucl', mincov=20, maxeval=1e-3):
	if prefix is None:
		prefix = inSeq
	d_besthit = {}
	for inHmmout in inHmmouts:
		for rc in HmmScan(inHmmout):
			suffix = rc.qname.split('|')[-1]
			if seqtype == 'nucl' and (suffix.startswith('aa') or suffix.startswith('rev_aa')):
				qid = '|'.join(rc.qname.split('|')[:-1])
			else:
				qid = rc.qname
			domain,clade = parse_hmmname(rc.tname, db=db)
			if db.startswith('rexdb'):
				cdomain = domain.split('-')[1]
				if cdomain == 'aRH':
					cdomain = 'RH'
				if cdomain == 'TPase':
					cdomain = 'INT'
				key = (qid, cdomain)
				if key in d_besthit:
					best_rc = d_besthit[key]
					if rc.score > best_rc.score:
						best_domain, _ = parse_hmmname(best_rc.tname, db=db)
						if domain == best_domain:
							d_besthit[key] = rc
						elif rc.envstart <= best_rc.envend and rc.envend >= best_rc.envstart: # overlap
							d_besthit[key] = rc
				else:
					d_besthit[key] = rc
			else: # gydb
				key = (qid, domain)
				if key in d_besthit:
					if rc.score > d_besthit[key].score:
						d_besthit[key] = rc
				else:
					d_besthit[key] = rc
	d_seqs = seq2dict(inSeq)
	lines = []
	for (qid, domain), rc in d_besthit.items():
		if rc.hmmcov < mincov or rc.evalue > maxeval:
			continue
		rawid = qid
		gene,clade = parse_hmmname(rc.tname, db=db)
		if db.startswith('rexdb'):
			domain = gene.split('-')[1]
		gid = '{}|{}'.format(qid, rc.tname)
		gseq = d_seqs[rc.qname].seq[rc.envstart-1:rc.envend]
		if seqtype == 'nucl':
			strand, frame = parse_frame(rc.qname.split('|')[-1])
			if strand == '+':
				nuc_start = rc.envstart * 3 - 2  + frame
				nuc_end = rc.envend* 3 + frame
			elif strand == '-':
				nuc_start = rc.qlen*3 - (rc.envend* 3 + frame) + 1
				nuc_end = rc.qlen*3 - (rc.envstart* 3 + frame) + 1
			else:
				nuc_start = rc.envstart
				nuc_end = rc.envend
		elif seqtype == 'prot':
			strand, frame = '+', '.'
			nuc_start, nuc_end, = rc.envstart, rc.envend
		match = re.compile(r'(\S+?):(\d+)\.\.(\d+)').match(qid)
		if match:
			qid, ltrstart, ltrend = match.groups()
			ltrstart = int(ltrstart)
			nuc_start = ltrstart + nuc_start - 1
			nuc_end = ltrstart + nuc_end -1
		attr = 'ID={};gene={};clade={};evalue={};coverage={};probability={}'.format(gid, domain, clade, rc.evalue, rc.hmmcov, rc.acc)
		gffline = [qid, 'TE_classifier', 'CDS', nuc_start, nuc_end, rc.score, strand, frame, attr, rc.evalue, rc.hmmcov, rc.acc, rawid, gid, gseq]
		lines.append(gffline)
	gff, seq, tsv = '{}.dom.gff3'.format(prefix), '{}.dom.faa'.format(prefix), '{}.dom.tsv'.format(prefix)
	fgff = open(gff, 'w')
	fseq = open(seq, 'w')
	ftsv = open(tsv, 'w')
	print >> ftsv, '\t'.join(['#id', 'length', 'evalue', 'coverge', 'probability'])
	for line in sorted(lines, key=lambda x: (x[0], x[-3], x[3])):
		gffline = line[:9]
		gffline = map(str, gffline)
		print >> fgff, '\t'.join(gffline)
		gid, gseq = line[-2:]
		gdesc = line[8]
		print >> fseq, '>{} {}\n{}'.format(gid, gdesc, gseq)
		evalue, hmmcov, acc = line[-6:-3]
		line = [gid, len(gseq), evalue, hmmcov, acc]
		print >> ftsv, '\t'.join(map(str, line))
	fgff.close()
	fseq.close()
	return gff, seq
def translate(inSeq, prefix=None):
	if prefix is None:
		prefix = inSeq
	outSeq = prefix + '.aa'
	with open(outSeq, 'w') as fp:
		six_frame_translate(inSeq, fp)
	return outSeq
def hmmscan(inSeq, hmmdb='rexdb.hmm', hmmout=None, ncpu=4):
	if hmmout is None:
		hmmout = prefix + '.domtbl'
	cmd = 'hmmscan --notextw -E 0.01 --domE 0.01 --noali --cpu {} --domtblout {} {} {} > /dev/null'.format(ncpu, hmmout, hmmdb, inSeq)
	run_cmd(cmd, logger=logger)
	return hmmout
def hmmscan_pp(inSeq, hmmdb='rexdb.hmm', hmmout=None, tmpdir='./tmp', processors=4):
	chunk_prefix = '{}/{}'.format(tmpdir, 'chunk_aaseq')
	_, _, _, chunk_files = split_fastx_by_chunk_num(
			inSeq, prefix=chunk_prefix, chunk_num=processors, seqfmt='fasta', suffix='')
	domtbl_files = [chunk_file + '.domtbl' for chunk_file in chunk_files]
	cmds = [ 
		'hmmscan --notextw -E 0.01 --domE 0.01 --noali --domtblout {} {} {}'.format(domtbl_file, hmmdb, chunk_file) \
			for chunk_file, domtbl_file in zip(chunk_files, domtbl_files)]
	jobs = pp_run(cmds, processors=processors)
	for cmd, (stdout, stderr, status) in zip(cmds, jobs):
		if not status == 0:
			logger.warn( "exit code {} for CMD '{}'".format(status, cmd) )
			logger.warn('\n\tSTDOUT:\n{0}\n\tSTDERR:\n{1}\n\n'.format(stdout, stderr))
	# cat files
	if hmmout is None:
		hmmout = prefix + '.domtbl'
	with open(hmmout, 'w') as f:
		for domtbl_file in domtbl_files:
			for line in open(domtbl_file):
				f.write(line)
	return hmmout

def LTRlibAnn(ltrlib, hmmdb='rexdb', seqtype='nucl', prefix=None, 
			force_write_hmmscan=False,
			processors=4, tmpdir='./tmp', 
			mincov=20, maxeval=1e-3):
	if prefix is None:
		prefix = '{}.{}'.format(ltrlib, hmmdb)
	
	if seqtype == 'nucl':
		logger.info( 'translating `{}` in six frames'.format(ltrlib) )
		tmp_prefix = '{}/translated'.format(tmpdir)
		aaSeq = translate(ltrlib, prefix=tmp_prefix)
	elif seqtype == 'prot':
		aaSeq = ltrlib
	
	logger.info( 'HMM scanning against `{}`'.format(DB[hmmdb]) )
	domtbl = prefix + '.domtbl'
	if not (os.path.exists(domtbl) and os.path.getsize(domtbl) >0) or force_write_hmmscan:
		if processors > 1:
			hmmscan_pp(aaSeq, hmmdb=DB[hmmdb], hmmout=domtbl, tmpdir=tmpdir, processors=processors)
		else:
			hmmscan(aaSeq, hmmdb=DB[hmmdb], hmmout=domtbl)
	else:
		logger.info( 'use existed non-empty `{}` and skip hmmscan'.format(domtbl) )
	logger.info( 'generating gene anntations' )
	gff, geneSeq = hmm2best(aaSeq, [domtbl], db=hmmdb, prefix=prefix, seqtype=seqtype, mincov=mincov, maxeval=maxeval)
	return gff, geneSeq
def replaceCls(ltrlib, seqtype='nucl', db='rexdb'):
	gff = ltrlib + '.' + db + '.gff3'
	if not os.path.exists(gff):
		gff, geneSeq, aaSeq = LTRlibAnn(ltrlib, seqtype=seqtype, hmmdb=db)
	annout = '{}.anno'.format(gff)
	newlib = '{}.reclassified'.format(ltrlib)
	fann = open(annout, 'w')
	flib = open(newlib, 'w')
	Classifier(gff, fout=fann).replace_annotation(ltrlib, fout=flib)
	fann.close()
	flib.close()
def parse_frame(string):
	if string.startswith('rev'):
		strand = '-'
	elif string.startswith('aa'):
		strand = '+'
	else:
		return '.', '.' #None,None
	frame = int(string[-1]) -1
	return strand, frame

class Dependency(object):
	def __init__(self,):
		pass
	def check(self):
		pass
	def check_hmmer(self, db, program='hmmscan'):
		dp_version = self.get_hmm_version(db)[:3]
		if self.check_presence(program):
			version0 = self.check_hmmer_verion(program)
			version = version0[:3]
			if version >= dp_version:
				logger.info('hmmer\t{}\tOK'.format(version0))
			elif version < dp_version:
				logger.warn('hmmer version {} is too low. Please update to {} from http://hmmer.org/download.html'.format(version, dp_version))
			else:
				logger.info('hmmer version {} is too high. You may use the version {}. However, I update the database first.'.format(version, dp_version))
				self.update_hmmer(db)
		else:
			logger.error('hmmer>={} not found'.format(dp_version))
	def get_hmm_version(self, db):
		line = open(db).readline()
		version = re.compile(r'HMMER\S+ \[([\w\.]+)').search(line).groups()[0]
		return version
	def update_hmmer(self, db):
		from small_tools import backup_file
		bk_db, db = backup_file(db)
		for suffix in ['.h3f', '.h3i', '.h3m', '.h3p']:
			backup_file(bk_db + suffix)
		cmd = 'hmmconvert {} > {}'.format(bk_db, db) 
		out, err, status0 = run_cmd(cmd, logger=logger)
		cmd = 'hmmpress {}'.format(db)
		out, err, status1 = run_cmd(cmd, logger=logger)
		if status0 + status1 == 0:
			logger.info('HMM converted. it will continue')
		else:
			logger.error('HMM failed to convert. exit')
			sys.exit(1)
	def check_blast(self, program='blastn'):
		if self.check_presence(program):
			version = self.check_blast_version(program)
			logger.info('{}\t{}\tOK'.format(program, version))
		else:
			logger.error('{} not found'.format(program))
	def check_presence(self, program):
		cmd = 'which {}'.format(program)
		out, err, status = run_cmd(cmd)
		if status == 0:
			return True
		else:
			return False
	def check_hmmer_verion(self, program):
		cmd = '{} -h'.format(program)
		out, err, status = run_cmd(cmd)
		version = re.compile(r'HMMER (\S+)').search(out).groups()[0]
		return version
	def check_blast_version(self, program):
		cmd = '{} -version'.format(program)
		out, err, status = run_cmd(cmd)
		version = re.compile(r'blast\S* ([\d\.\+]+)').search(out).groups()[0]
		return version
def main():
	subcmd = sys.argv[1]
	if subcmd == 'LTRlibAnn':   # hmmscan + HmmBest
		ltrlib = sys.argv[2]	# input is LTR library (fasta)
		try: 
			hmmdb = sys.argv[3] # rexdb, gydb, pfam, etc.
			try: seqtype = sys.argv[4]
			except IndexError: seqtype = 'nucl'
			LTRlibAnn(ltrlib, hmmdb=hmmdb, seqtype=seqtype)
		except IndexError:
			LTRlibAnn(ltrlib)
	elif subcmd == 'HmmBest':
		inSeq = sys.argv[2]          # input: LTR library (translated protein)
		prefix = inSeq
		inHmmouts = sys.argv[3:]     # input: hmmscan output (inSeq search against hmmdb)
		hmm2best(inSeq, inHmmouts, prefix)
	elif subcmd == 'Classifier':
		gff = sys.argv[2]	    # input: gff3 output by LTRlibAnn or HmmBest
		try: db = sys.argv[3]	# rexdb or gydb
		except IndexError: db = 'rexdb'
		for line in Classifier(gff, db=db):
			continue
	elif subcmd == 'replaceCls':	# LTRlibAnn + Classifier
		ltrlib = sys.argv[2]	    # input: LTR library (nucl fasta) 
		replaceCls(ltrlib)
	elif subcmd == 'replaceClsLR':
		genome = sys.argv[2]		# input: genome input for LTR_retriever pipeline
		Retriever(genome).re_classify()
	else:
		raise ValueError('Unknown command: {}'.format(subcmd))

if __name__ == '__main__':
	#main()
	pipeline(Args())
