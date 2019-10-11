import sys
import re
class RMOutRecord():
	def __init__(self, line):
		d_class = dict([
            ('LTR', 'ClassI'),
            ('LINE', 'ClassI'),
            ('SINE', 'ClassI'),
            ('DNA', 'ClassII'),
            ('RC', 'ClassII'),
            ('Unknown', 'Others'),
            ('Satellite', 'Others'),
            ('Simple_repeat', 'Others'),
            ('Low_complexity', 'Others'),
            ('rRNA', 'Others'),
            ('snRNA', 'Others'),
            ])

		temp = line.strip().split()
		title = ['score', 'perc_div', 'perc_del', 'perc_ins',
				 'query_id', 'query_begin', 'query_end', 'query_left', 'strand',
				 'repeat_family', 'super_class', 'repeat_begin', 'repeat_end', 'repeat_left', 'repeat_id',
				]
		convert = [int, float, float, float,
				   str, int, int, str, str,
				   str, str, str, int, str, str ]
		assert len(temp) == len(title) or (len(temp)-1 == len(title) and temp[-1] == "*") or (len(temp) == len(title)-1)
#		try: assert len(temp) == len(title) or (len(temp)-1 == len(title) and temp[-1] == "*")
#		except AssertionError: print >> sys.stderr, temp, '\n', title
		self.__dict__ = {key: func(value) for key,value,func in zip(title, temp, convert)}
		self.query_left = int(self.query_left.strip('()'))
		self.repeat_begin = int(self.repeat_begin.strip('()'))
		self.repeat_left = int(self.repeat_left.strip('()'))
		if self.strand == 'C':
			self.strand = '-'
		if temp[-1] == "*": # there is a higher-scoring match whose domain partly (<80%) includes the domain of this match
			self.overlap = True
		else:
			self.overlap = False
		self.order = self.super_class.split('/')[0]
		try: self.superfamily = self.super_class.split('/')[1]
		except IndexError: self.superfamily = ''
		try: self.Class = d_class[self.order]
		except KeyError: self.Class = 'Others'
		self.query_length = self.query_end + self.query_left
		self.query_match_length = self.query_begin - self.query_end + 1
		if self.strand == '-':
			self.repeat_begin, self.repeat_end, self.repeat_left = self.repeat_left, self.repeat_end, self.repeat_begin 
		self.target = '{} {} {}'.format(self.repeat_family, self.repeat_begin, self.repeat_end)
	def write(self, fout=sys.stdout):
		if self.strand == '+':
			line = [self.score, self.perc_div, self.perc_del, self.perc_ins, \
				self.query_id, self.query_begin, self.query_end, '({})'.format(self.query_left), self.strand, \
				self.repeat_family, self.super_class, \
				self.repeat_begin, self.repeat_end, '({})'.format(self.repeat_left), self.repeat_id]
		elif self.strand == '-':
			line = [self.score, self.perc_div, self.perc_del, self.perc_ins, \
                self.query_id, self.query_begin, self.query_end, '({})'.format(self.query_left), self.strand, \
                self.repeat_family, self.super_class, \
                '({})'.format(self.repeat_left), self.repeat_end, self.repeat_begin, self.repeat_id]
		if self.overlap:
			line += ['*']
		line = map(str, line)
		print >> fout, '\t'.join(line)
	def write_gff3(self, fout=sys.stdout):
		attr = 'Target={target}'.format(self.__dict__)
		line = [self.query_id, 'RepeatMasker', 'dispersed_repeat', self.query_begin, self.query_end, \
				self.score, self.strand, '.', attr]
		line = map(str, line)
		print >> fout, '\t'.join(line)
	def get_seq(self, seqRecord):
		id = '{}:{}..{}|{}#{}'.format(self.query_id, self.query_begin, self.query_end, self.repeat_family, self.super_class)
		teRecord = seqRecord[self.query_begin-1:self.query_end]
		teRecord.id = id
		teRecord.description = id
		return teRecord

class RMOutParser():
	def __init__(self, inRMout):
		self.inRMout = inRMout
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.inRMout):
			if not re.compile(r'.*\d').match(line): # title
				continue
			yield RMOutRecord(line)
	def get_seqs(self, genome, fout=sys.stdout):
		from Bio import SeqIO
		d_seqs = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
		for rc in self._parse():
			seqRecord = d_seqs[rc.query_id]
			SeqIO.write(rc.get_seq(seqRecord), fout, 'fasta')

def RMOut2mnd(inRMout, outMnd, binsize0=300, scale=1, mapq=20, minlen=1000):
	from itertools import combinations
	d_family = {}
	for rc in RMOutParser(inRMout):
		if rc.query_match_length < minlen:
			continue
		try: d_family[rc.repeat_family] += [rc]
		except KeyError: d_family[rc.repeat_family] = [rc]
	i = 0
	for repeat_family, records in d_family.iteritems():
		for rc1, rc2 in combinations(records, 2):
			binsize = int(binsize0 * scale)
			bins = max(1, (rc1.query_match_length+rc2.query_match_length) / 2 / binsize)
			interval1 = (rc1.query_end - rc1.query_begin) / bins
			interval2 = (rc2.query_end - rc2.query_begin) / bins
			bins1 = range(rc1.query_begin, rc1.query_end, interval1)
			bins2 = range(rc2.query_begin, rc2.query_end, interval2)
			if rc1.strand == '-':
				bins1.reverse()
			if rc2.strand == '-':
				bins2.reverse()
			str1, str2 = 0,0
			chr1, chr2 = rc1.query_id, rc2.query_id
			frag1, frag2 = 0, 1
			mapq1, mapq2 = mapq, mapq
			cigar1,sequence1,cigar2,sequence2,readname1,readname2 = ['-'] * 6
			for bin1, bin2 in zip(bins1,bins2):
				pos1, pos2 = bin1, bin2
				line = [str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, mapq1,cigar1,sequence1,mapq2,cigar2,sequence2,readname1,readname2]
				line = map(str, line)
				print >> outMnd, ' '.join(line)
				i += 1
	print >>sys.stderr, '{} links'.format(i)

def main():
	subcmd = sys.argv[1]
	if subcmd == 'out2mnd':
		inRMout = sys.argv[2]
		outMnd = sys.stdout
		RMOut2mnd(inRMout, outMnd)
	elif subcmd == 'out2seqs':
		inRMout = sys.argv[2]
		genome = sys.argv[3]
		RMOutParser(inRMout).get_seqs(genome)
	else:
		raise ValueError('Unknown command: {}'.format(subcmd))

if __name__ =='__main__':
    main()

