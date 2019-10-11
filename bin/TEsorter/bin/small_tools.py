import sys
import os

import gzip
import re
import time
import shutil
import subprocess
from Bio import SeqIO
ISOTIMEFORMAT='%Y-%m-%d %X'

def open_file(infile, mode='r'):
	suffix = os.path.splitext(infile)[-1]
	if suffix == '.gz':
		import gzip
		return gzip.open(infile, mode)
	elif suffix == '.bz2':
		import bz2
		return bz2.BZ2File(infile, mode)
	else:
		return open(infile, mode)

def get_hex_colors(n):
	import matplotlib.pyplot as plt
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	values = range(n)
	jet = cm = plt.get_cmap('jet')
	cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	colorVal = [scalarMap.to_rgba(v) for v in values]
	return [colors.to_hex(v) for v in colorVal]
def getHtml(url):
        import urllib2,socket
        socket.setdefaulttimeout(5)
        socket.setdefaulttimeout(5)
        try:
                page = urllib2.urlopen(url)
                html = page.read()
                return html
        except socket.timeout:
                print 'time out, trying again!'
                time.sleep(5)
                return getHtml(url)
        except urllib2.HTTPError, e:
                if e.code == 404:
                        return e.code
                print 'HTTPError code: ', e.code, ', trying again!'
                time.sleep(5)
                return getHtml(url)
        except urllib2.URLError, e:
                print 'URLError reason',e.reason,', trying again!'
                time.sleep(5)
                return getHtml(url)
        except:
                'UnknownError, trying again!'
                time.sleep(5)
                return getHtml(url)

def remove_short_seqs(inSeq, outSeq, minLen=200, format='fasta'):
	'''Remove sequences shorter than the cutoff.'''
	f = open(outSeq, 'w')
	for record in SeqIO.parse(inSeq, format):
		if len(record.seq) >= minLen:
			SeqIO.write(record, f, format)
		else:
			continue
	f.close()

def fq_switch(inFq):
	if os.path.exists(inFq):
		return inFq
	elif os.path.exists(inFq+'.gz'):
		return inFq+'.gz'
	else:
		raise IOError('File %s is NOT exists' % inFq)

def mkdirs(DIR):
	'''compatible while multiprocessing with pp module'''
	if os.path.exists(DIR):
		pass
	else:
		try: os.makedirs(DIR)
		except OSError: pass

def is_complete_cds(Seq, translate_table=1):
	'''A Bio.Seq object (CDS) is complete (both start and end codons present)?'''
	try:
		pro_seq = Seq.translate(table=translate_table, to_stop=True, cds=True)
		return True
	except CodonTable.TranslationError:
		return False

def time_convert(number):
	'''Simply convert format of date.
usage:
import time
ISOTIMEFORMAT='%Y-%m-%d %X'
start = time.time()
end = time.time()
time_convert(end-start)'''
	hour = int(number/3600)
	min = int((number-3600*hour)/60)
	sec = number - 3600*hour - 60*min
	return '%sh%sm%0.2fs' %(hour,min,sec)

def is_gz(input_file):
	'''Is a gzip file?'''
	if os.path.splitext(input_file)[-1] == '.gz':
		return True
	else:
		return False

def table2xls(table_file, xls_file = None):
	import xlwt
	# init
	if not xls_file:
		xls_file = table_file + '.xls'

	wb = xlwt.Workbook()
	ws1= wb.add_sheet('Sheet1')

	i = 0
	for line in open(table_file,'r'):
		temp = line.rstrip().split('\t')
		j = 0
		for value in temp:
			# try: value = value.strip().decode('gbk')
			# except UnicodeDecodeError: value = 'NULL'
			ws1.write(i,j,value)
			j += 1
		i += 1
	wb.save(xls_file)

def combine_tabs_2xls(table_files, xls_file = None, sheets = None):
	import xlwt
	# init
	if not xls_file:
		xls_file = table_file[0] + '.xls'
	if not sheets:
		sheets = table_file

	wb = xlwt.Workbook()
	for table_file, sheet in zip(table_files, sheets):
		ws1= wb.add_sheet(sheet)
		i = 0
		for line in open(table_file,'r'):
			temp = line.rstrip().split('\t')
			j = 0
			for value in temp:
				ws1.write(i,j,value)
				j += 1
			i += 1
	wb.save(xls_file)


def run_time(func):
	def _run_time():
		start = time.time()
		logFile = 'main.logfile'
		f = open(logFile,'a')
		print >>f, 'Start at %s ...' % time.strftime(ISOTIMEFORMAT, time.localtime())
		print >>f, 'Runing commands:'

		print >>f, '  %s' % func()

		print >>f, 'End   at %s .' % time.strftime(ISOTIMEFORMAT, time.localtime())
		end = time.time()
		print >>f, 'Total Time Used: %s\n' % time_convert(end-start)
		f.close()
	return _run_time

def bk_not_overwrite(input_file_bk):
	if os.path.exists(input_file_bk):
		return bk_not_overwrite(input_file_bk + '.1')
	else:
		return input_file_bk

def backup_file(input_file):
	input_file_bk = input_file + '.bk'
	input_file_bk = bk_not_overwrite(input_file_bk)
	shutil.move(input_file, input_file_bk)
	return (input_file_bk, input_file)

def count_record(inFile, format):
	if is_gz(inFile):
		cat = 'zcat'
	else:
		cat = 'cat'
	if format in set(['fastq', 'fq']):
		cmd = '%s "%s"| grep -c "^+$"' % (cat, inFile,)
	elif format in set(['fasta', 'fa', 'fna', 'faa', 'fas']):
		cmd = '%s "%s"| grep -c "^>"' % (cat, inFile,)
	elif format in set(['text', 'txt']):
		cmd = '%s "%s"| wc -l' % (cat, inFile,)
	job = subprocess.Popen(cmd,stdout=subprocess.PIPE,\
					stderr=subprocess.PIPE,shell=True)
	output = job.communicate()
	try: return int(output[0])
	except ValueError: print output

def flattern(nested):
	try:
		for sublist in nested:
			for element in flattern(sublist):
				yield element
	except TypeError:
		yield nested
def flattern2(nested):
	for sublist in nested:
		for element in sublist:
			yield element
