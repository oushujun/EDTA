#!/bin/env python
# coding: utf-8
'''# Author: zrg1989@qq.com
# Version: 0.1
'''
import sys
import glob
import os, re
from math import log
from collections import Counter, OrderedDict
from Bio import SeqIO

bindir = os.path.dirname(os.path.realpath(__file__))
DB = {
	'gydb' : bindir + '/database/GyDB2.hmm',
	'rexdb': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm',
	'rexdb-plant': bindir + '/database/REXdb_protein_database_viridiplantae_v3.0.hmm',
	'rexdb-metazoa': bindir + '/database/REXdb_protein_database_metazoa_v3.hmm',
	}

class IntactRecord():
	def __init__(self, title, temp):
		self.dict = dict([(key, value) for key, value in zip(title, temp)])
#		print self.dict
		self.dict['Identity'] = float(self.dict['Identity'])
		self.dict['Insertion_Time'] = int(self.dict['Insertion_Time'])
		for key, value in self.dict.items():
			try: exec 'self.{} = value'.format(key)
			except SyntaxError: pass
		self.chr, self.start, self.end = re.compile(r'(\S+?):(\d+)\.\.(\d+)').match(self.LTR_loc).groups()
		self.start, self.end = int(self.start), int(self.end)
		
class CandidateRecord():
	def __init__(self, title, temp, mu=1.3e-8):
		if len(temp) < len(title):
			temp += [None] * (len(title)-len(temp))
		self.dict = dict([(key, value) for key, value in zip(title, temp)])
		for key, value in self.dict.items():
			try: self.dict[key] = int(value)
			except: continue
		try: self.dict['similarity'] = float(self.dict['similarity'])
		except ValueError: self.dict['similarity'] = 0
		try: self.dict['ageya'] = int(self.dict['ageya'])
		except ValueError:
			if self.dict['similarity'] <= 0.75:
				self.dict['ageya'] = None
			else:
	#				d = -3.0/4*log(1-4.0/3*(1-self.dict['similarity']/100))
	#				self.dict['ageya'] = d / (2*mu)
				self.dict['ageya'] = None
		except:
			self.dict['ageya'] = None
		for key, value in self.dict.items():
			try: exec 'self.{} = value'.format(key)
			except SyntaxError: pass
		self.INT_str, self.INT_end = self.lLTR_end+1, self.rLTR_str-1
	def write(self, fout):
		self.line = [self.start, self.end, self.len, self.lLTR_str, self.lLTR_end, self.lLTR_len, 
				self.rLTR_str, self.rLTR_end, self.rLTR_len, self.similarity, 
				self.seqid, self.chr, self.direction, self.TSD, self.lTSD, self.rTSD, 
				self.motif, self.superfamily, self.family, self.ageya]
		self.line = ['' if value is None else str(value) for value in self.line]
		print >> fout, '\t'.join(self.line)
class Retriever():
	def __init__(self, genome):
		self.genome = genome
		if glob.glob(genome+'.mod.*'):
			self.genome = self.genome + '.mod'
		self.pass_list = self.genome + '.pass.list'
		self.pass_gff3 = self.genome + '.pass.gff3'
		self.nmtf_pass_list = self.genome + '.nmtf.pass.list'
		self.pass_lists = [self.pass_list, self.nmtf_pass_list]
		self.retriever_all_scn = self.genome + '.retriever.all.scn'
		self.ltrlib = self.genome + '.LTRlib.fa'
		self.retriever_all_scn2 = self.retriever_all_scn + '2'
		if not (os.path.exists(self.retriever_all_scn2) \
		   and os.path.getsize(self.retriever_all_scn2) > 1000):
			print >>sys.stderr, 're-organize', self.retriever_all_scn
			self.re_scn()
		self.retriever_all_scn = self.retriever_all_scn2
	def get_full_seqs(self, fout=sys.stdout):
		d_seqs = seq2dict(self.genome)
		for rc in self.intact_list():
			ltr_seq = d_seqs[rc.chr].seq[rc.start-1:rc.end]
			ltr_cls = '{}/{}'.format(rc.TE_type, rc.SuperFamily)
			print >> fout, '>{}#{}\n{}'.format(rc.LTR_loc, ltr_cls, ltr_seq)
	def re_scn(self): # remove redundant
		idmap = self.seqIdmap
		lrt_set = set([])
		f = open(self.retriever_all_scn2, 'w')
		i,j,k = 0,0,0
		for line in open(self.retriever_all_scn):
			temp = line.strip().split()
			if line.startswith('#'):
				f.write(line)
			if line.startswith('#start'):
				temp = line.strip().strip('#()').replace('(', '').split()
				title = temp
				continue
			elif line.startswith('#'):
				continue
			rc = CandidateRecord(title, temp)
			i += 1
			if rc.chr is None or rc.chr == 'NA':
				j += 1
				rc.chr = idmap[rc.seqid]
			key = (rc.chr, rc.start, rc.end, rc.lLTR_str, rc.lLTR_end, rc.rLTR_str, rc.rLTR_end)
			if key in lrt_set:
				k += 1
				continue
			lrt_set.add(key)
			rc.write(f)
		f.close()
		print >>sys.stderr, '{} total {}, {} without chr, {} discarded, {} retained'.format(self.retriever_all_scn, i, j, k, i-k)
	@property
	def seqIdmap(self):
		i = 0
		d = {}
		for rc in SeqIO.parse(self.genome, 'fasta'):
			d[i] = rc.id
			i += 1
		return d
	def intact_list(self):
		for pass_list in self.pass_lists:
			for line in open(pass_list):
				temp = line.strip().split('\t')
				if line.startswith('#'):
					temp = line.strip().split()
					temp[0] = temp[0].strip('#')
					title = temp
					continue
				yield IntactRecord(title, temp)
	def all_scn(self):
		for line in open(self.retriever_all_scn):
			temp = line.strip().split()
			if line.startswith('#start'):
				temp = line.strip().strip('#()').replace('(', '').split()
				title = temp
				continue
			elif line.startswith('#'):
				continue
			yield CandidateRecord(title, temp)
	def re_classify(self, seqtype='dna', db='rexdb'):
		ltrlib = self.ltrlib
		gff, geneSeq, aaSeq = LTRlibAnn(ltrlib, seqtype=seqtype, db=db)
		gff = ltrlib + '.' + db + '.gff3'
		annout = '{}.anno'.format(gff)
		newlib = '{}.reclassified'.format(ltrlib)
		fann = open(annout, 'w')
		flib = open(newlib, 'w')
		Classifier(gff, fout=fann).replace_annotation(ltrlib, fout=flib, idmap=self.ltr_map)
		fann.close()
		flib.close()
	@property
	def ltr_map(self):
		d = {}
		for rc in self.all_scn():
#			print rc.chr, rc.lLTR_str, rc.lLTR_end
			LTR1 = '{}:{}..{}_LTR'.format(rc.chr, rc.lLTR_str, rc.lLTR_end)
			LTR2 = '{}:{}..{}_LTR'.format(rc.chr, rc.rLTR_str, rc.rLTR_end)
			INT  = '{}:{}..{}_INT'.format(rc.chr, rc.INT_str, rc.INT_end)
			for id in [LTR1, LTR2, INT]:
				d[id] = INT
		return d

def InsertionTimePlot(genome, type, mu=1.3e-8):
	if type == 'intact':
		outfig = genome + '.Intact.Insertion_Time.pdf'
		records = Retriever(genome).intact_list()
		typeStr = '{TE_type}/{SuperFamily}'
		timeStr = 'Insertion_Time'
	elif type == 'candidate':
		outfig = genome + '.Candidate.Insertion_Time.pdf'
		records = Retriever(genome).all_scn()
		typeStr = '{superfamily}/{family}'
		timeStr = 'ageya'
	else:
		raise ValueError('Unknown type: {}'.format(type))
	tmpfile = genome + '.pass.insert_time'	
	f = open(tmpfile, 'w')
	line = ['TE_Type', 'Insertion_Time']
	print >>f, '\t'.join(line)
	for rc in records:
		Type = typeStr.format(**rc.dict)
		if rc.dict[timeStr] is None:
			continue
		Insertion_Time = rc.dict[timeStr] / 1e6 * (1.3e-8 / mu) # Mya
		line = [Type, Insertion_Time]
		line = map(str, line)
		print >>f, '\t'.join(line)
	f.close()
	
	# plot
	r_src = '''
data <- read.table('{}', head=T)
library(ggplot2)
p <- ggplot(data, aes(x=Insertion_Time, color=TE_Type)) + geom_line(stat="density")
ggsave('{}', p)
'''.format(tmpfile, outfig)
	r_file = tmpfile + '.r'
	with open(r_file, 'w') as f:
		print >>f, r_src
	cmd = 'Rscript {}'.format(r_file)
	os.system(cmd)
class Classifier():
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
		line = ['#TE', 'Superfamily', 'Family', 'Clade', 'Code', 'Strand', 'hmmmatchs']
		print >> self.fout, '\t'.join(line)
		for rc in self.parse():
			rc_flt = rc #[line for line in rc if line.gene in self.markers]
			#if len(rc_flt) == 0:
			#	rc_flt = [line for line in rc]
			strands = [line.strand for line in rc_flt]
			if len(set(strands)) >1 :
				strand = '?'
			else:
				strand = strands[0]
			if strand == '-':
				rc_flt.reverse()
				rc.reverse()
			lid = rc_flt[0].ltrid
			domains = ' '.join(['{}|{}'.format(line.gene, line.clade)  for line in rc])
			genes  = [line.gene  for line in rc_flt]
			clades = [line.clade for line in rc_flt]
			names = [line.name for line in rc_flt]
			if self.db.startswith('rexdb'):
				order, superfamily, max_clade, coding = self.identify_rexdb(genes, names)
			elif self.db == 'gydb':
				order, superfamily, max_clade, coding = self.identify(genes, clades)
			line = [lid, order, superfamily, max_clade, coding, strand, domains]
			print >> self.fout, '\t'.join(line)
			self.ltrid, self.order, self.superfamily, self.clade, self.code, self.strand, self.domains = line
			yield self
	def identify_rexdb(self, genes, clades):
		perfect_structure = {
            ('LTR', 'Copia'): ['Ty1-GAG', 'Ty1-PROT', 'Ty1-INT', 'Ty1-RT', 'Ty1-RH'],
            ('LTR', 'Gypsy'): ['Ty3-GAG', 'Ty3-PROT', 'Ty3-RT', 'Ty3-RH', 'Ty3-INT'],
			}
		clade_count = Counter(clades)
		max_clade = max(clade_count, key=lambda x: clade_count[x])
		order, superfamily = self._parse_rexdb(max_clade)
		if len(clade_count) == 1:
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
				coding = 'cmpl' # completed gene structure
			else:
				coding = 'lost'
		except KeyError:
			coding = 'unknown'
		return order, superfamily, max_clade, coding
	def _parse_rexdb(self, clade): # full clade name
		if clade.startswith('Class_I/LTR/Ty1_copia'):
			order, superfamily = 'LTR', 'Copia'
		elif clade.startswith('Class_I/LTR/Ty3_gypsy'):
			order, superfamily = 'LTR', 'Gypsy'
		elif clade.startswith('Class_I/'): # LINE, pararetrovirus, Penelope, DIRS
			order, superfamily = clade.split('/')[1], 'unknown'
		elif clade.startswith('Class_II/'):
			try: order, superfamily = clade.split('/')[2:4]
			except ValueError: order, superfamily = clade.split('/')[2], 'unknown'
		return order, superfamily
	def identify(self, genes, clades):
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
			print >>sys.stderr, 'unknown clade: {}'.format(max_clade)
		try:
			ordered_genes = perfect_structure[(order, superfamily)]
			my_genes = [gene for gene in genes if gene in set(ordered_genes)]
			if ordered_genes == my_genes:
				coding = 'cmpl' # completed gene structure
			else:
				coding = 'lost'
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
					print >>sys.stderr, '[Warning] skipped', e
			if intid in d_class:
				neword, newfam = d_class[intid]
				re_org = self.re_orgnize(rc.id, neword, newfam)
				if re_org:
					i += 1
					rc.id = re_org
#					print >>sys.stderr, rc.description, len(rc.seq)
			SeqIO.write(rc, fout, 'fasta')
		print >> sys.stderr, i, 'sequences re-classified'
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
			if self.clade == '17.6': # exception
				self.clade = '17_6'
			self.superfamily = self.dict['Family'].split('/')[-1]
			if self.superfamily == 'Retroviridae':	# deltaretroviridae gammaretroviridae
				self.clade = self.dict['Cluster_or_genus'].replace('virus', 'viridae')
			if self.superfamily == 'Retrovirus':	# an exception
				self.superfamily = 'Retroviridae'
			self.order = 'LTR' if self.dict['System'] in {'LTR_retroelements', 'LTR_Retroelements', 'LTR_retroid_elements'} else self.dict['System']
			self.clade = self.clade.replace('-', '_') # A-clade V-clade C-clade
			
			yield self
			self.clade = self.clade.lower()
			yield self
		self.order, self.superfamily, self.clade, self.dict = ['LTR', 'Copia', 'ty1/copia', {}]  # AP_ty1/copia
		yield self
		for clade, order in zip(['retroelement', 'shadow', 'all'], ['LTR', 'Unknown', 'Unknown']): # CHR
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

class HmmScan():
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
class HmmDomRecord():
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
	from Bio import SeqIO
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

def hmm2best(inSeq, inHmmouts, prefix=None, db='rexdb', seqtype='dna', mincov=20, maxeval=1e-3):
	if prefix is None:
		prefix = inSeq
	d_besthit = {}
	for inHmmout in inHmmouts:
		for rc in HmmScan(inHmmout):
			suffix = rc.qname.split('|')[-1]
			if suffix.startswith('aa') or suffix.startswith('rev_aa'):
				qid = '|'.join(rc.qname.split('|')[:-1])
			else:
				qid = rc.qname
			domain,clade = parse_hmmname(rc.tname, db=db)
			if db.startswith('rexdb'):
				cdomain = domain.split('-')[1]
				if cdomain == 'aRH':
					cdomain = 'RH'
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
			else:
				key = (qid, domain)
				if key in d_besthit:
					if rc.score > d_besthit[key].score:
						d_besthit[key] = rc
				else:
					d_besthit[key] = rc
#	print d_besthit
	d_seqs = seq2dict(inSeq)
	lines = []
	for (qid, domain), rc in d_besthit.items():
		if rc.hmmcov < mincov or rc.evalue > maxeval:
			continue
#		gid = '{}|{}|{}'.format(qid, domain, rc.tname)
		rawid = qid
#		clade = '_'.join(rc.tname.split('_')[1:])
		gene,clade = parse_hmmname(rc.tname, db=db)
		if db.startswith('rexdb'):
			domain = gene
		gid = '{}|{}'.format(qid, rc.tname)
		gseq = d_seqs[rc.qname].seq[rc.envstart-1:rc.envend]
		if seqtype == 'dna':
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
		gffline = [qid, 'ltrapl', 'CDS', nuc_start, nuc_end, rc.score, strand, frame, attr, rc.evalue, rc.hmmcov, rc.acc, rawid, gid, gseq]
		lines.append(gffline)
	gff, seq, tsv = '{}.gff3'.format(prefix), '{}.faa'.format(prefix), '{}.tsv'.format(prefix)
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
	prog = 'perl {}/bin/Six-frame_translate.pl'.format(bindir)
	outSeq = prefix + '.aa'
	cmd = '{} {} > {}'.format(prog, inSeq, outSeq)
	os.system(cmd)
	return outSeq
def hmmscan(inSeq, hmmdb='rexdb', prefix=None):
	if prefix is None:
		prefix = inSeq
	outDomtbl = prefix + '.domtbl'
	cmd = 'hmmscan --notextw -E 0.01 --domE 0.01 --noali --cpu 4 --domtblout {} {} {} > /dev/null'.format(outDomtbl, hmmdb, inSeq)
	os.system(cmd)
	return outDomtbl
def LTRlibAnn(ltrlib, seqtype='dna', hmmdb='rexdb', prefix=None):
#	ltrlib = Retriever(genome).ltrlib
	if prefix is None:
		prefix = '{}.{}'.format(ltrlib, hmmdb)
	if seqtype == 'dna':
		print >>sys.stderr, 'translating {} in six frames'.format(ltrlib)
		aaSeq = translate(ltrlib)
#		aaSeq = ltrlib + '.aa'
	elif seqtype == 'prot':
		aaSeq = ltrlib
	print >>sys.stderr, 'HMM scanning against {}'.format(DB[hmmdb])
	domtbl = hmmscan(aaSeq, hmmdb=DB[hmmdb], prefix=prefix)
#	domtbl = prefix + '.domtbl'
	print >>sys.stderr, 'generating gene anntations'
	gff, geneSeq = hmm2best(aaSeq, [domtbl], db=hmmdb, prefix=prefix, seqtype=seqtype)
	return gff, geneSeq, aaSeq
def replaceCls(ltrlib, seqtype='dna', db='rexdb'):
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
					
def main():
	subcmd = sys.argv[1]
	if subcmd == 'InsertionTimePlot':
		genome = sys.argv[2]
		try: type = sys.argv[3]
		except IndexError: type = 'intact'
		try: mu = float(sys.argv[4])
		except IndexError: mu=1.3e-8
		InsertionTimePlot(genome, type, mu=mu)
	elif subcmd == 'LTRlibAnn': # hmmscan + HmmBest
		ltrlib = sys.argv[2]	# input is LTR library (fasta)
		try: 
			hmmdb = sys.argv[3] # rexdb, gydb, pfam, etc.
			try: seqtype = sys.argv[4]
			except IndexError: seqtype = 'dna'
			LTRlibAnn(ltrlib, hmmdb=hmmdb, seqtype=seqtype)
		except IndexError:
			LTRlibAnn(ltrlib)
	elif subcmd == 'HmmBest':
		inSeq = sys.argv[2] # aa seq # input: LTR library (translated protein)
		prefix = inSeq
		inHmmouts = sys.argv[3]     # input: hmmscan output (inSeq search against hmmdb)
		try: db = sys.argv[4]
		except IndexError: db = 'rexdb'
		try: seqtype = sys.argv[5]
		except IndexError: seqtype = 'dna'
		hmm2best(inSeq, [inHmmouts], prefix, db=db, seqtype=seqtype)
	elif subcmd == 'Classifier':
		gff = sys.argv[2]	# input: gff3 output by LTRlibAnn or HmmBest
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
	elif subcmd == 'get_full_seqs':
		genome = sys.argv[2]
		Retriever(genome).get_full_seqs()
	else:
		raise ValueError('Unknown command: {}'.format(subcmd))
if __name__ == '__main__':
	main()
