#!/bin/env python
#coding utf-8
'''RUN system CoMmanDS in Multi-Processing'''
import sys
import os, stat
import shutil
#import psutil
import subprocess
from optparse import OptionParser
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s -%(levelname)s- %(message)s')
logger = LOGGER = logging.getLogger(__name__)

try:
	import pp
except ImportError as e:
	logger.warn('{}\nparallel computing is not available'.format(e))
try: 
	import drmaa    # for grid
	from tempfile import NamedTemporaryFile
except ImportError as e:
	logger.warn('{}\ngrid computing is not available'.format(e))

__version__ = '1.0'

class Grid(object):
	def __init__(self, cmd_list=None, 
			work_dir=None, out_path=None, 
			err_path=None, grid_opts='',
			cpu=1, mem='1g', template=None, 
			tc_tasks = None, script=None,
			join_files=True):
		self.cmd_list = cmd_list
		self.grid = self.which_grid()
		self.grid_opts = grid_opts
		self.template = template
#		if tc_tasks is not None:
#			if self.grid == 'sge':
#				self.grid_opts += ' -tc {}'.format(tc_tasks)
#		if self.grid == 'sge':
		self.grid_opts = self.grid_opts.format(cpu=cpu, mem=mem, tc=tc_tasks)
		if self.cmd_list is not None:
			self.script = script
			if self.script is None:
				fp = NamedTemporaryFile('w+t', delete=False)
				self.script = fp.name
			else:
				fp = open(self.script, 'w')
			self.make_script(fp)
			fp.close()
			os.chmod(self.script, stat.S_IRWXU)
			self.script = os.path.realpath(self.script)
			self.work_dir = work_dir
			if self.work_dir is None:
				self.work_dir = os.getcwd()
			self.out_path = out_path	# The path to a file representing job's stdout.
			self.err_path = err_path	# The path to a file representing job's stderr
			if self.out_path is None:
				self.out_path = ':' + self.script + '.out'
			elif not self.out_path.startswith(':'):
				self.out_path = ':' + self.out_path
			if self.err_path is None:
				self.err_path = ':' + self.script + '.err'
			elif not self.err_path.startswith(':'):
				self.err_path = ':' + self.err_path
			self.join_files = join_files	# True if stdin and stdout should be merged, False otherwise.
			
	def make_script(self, fout=sys.stdout):
		print >> fout, '#!/bin/bash'
		if self.template is None and self.grid == 'sge':
#			print >> fout, '#$ {}'.format(self.grid_opts)
			self.template = 'if [ $SGE_TASK_ID -eq {id} ]; then\n{cmd}\nfi'
		for i, cmd in enumerate(self.cmd_list):
			grid_cmd = self.template.format(id=i+1, cmd=cmd)
			print >> fout, grid_cmd
	def submit(self):
		job_status = []
		s = drmaa.Session()
		s.initialize()
		jt = s.createJobTemplate()
		jt.jobEnvironment = os.environ.copy()
		jt.nativeSpecification = self.grid_opts
		jt.remoteCommand = self.script
		jt.workingDirectory = self.work_dir
		jt.outputPath = self.out_path
#		jt.errorPath = self.err_path
		logger.info('submiting {}'.format([self.script, self.grid_opts]))
		jt.joinFiles = self.join_files
		joblist = s.runBulkJobs(jt, 1, len(self.cmd_list), 1)
		jobid = joblist[0].split('.')[0]
		_qsub_log(jobid, self.work_dir, self.script)
		#jid = s.runJob(jt)
		s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
		logger.info('waiting for {} tasks'.format(len(joblist)))
		for curjob in joblist:
			try:
				retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
			except Exception as e:
				logger.warn('Job {}: {}'.format(curjob, e))
				job_status += [(None, 1)]
				continue
			logger.info('Task: {0} finished with status {1}, exit code {2}'.format(
					retval.jobId, retval.hasExited, retval.exitStatus))
			task_id, status = retval.jobId, (not retval.hasExited) + retval.exitStatus
			job_status += [(task_id, status)]
		s.deleteJobTemplate(jt)
		s.exit()
		return job_status
	def which_grid(self):
		with drmaa.Session() as s:
			name = s.drmaaImplementation
		grid, version = name.split()
		if grid == 'OGS/GE':
			return 'sge'
def run_tasks(cmd_list, tc_tasks=None, mode='grid', grid_opts='', cpu=1, mem='1g', cont=1,
			retry=1, script=None, out_path=None, completed=None, cmd_sep='\n', template=None, **kargs):
	if not cmd_list:
		logger.info('cmd_list with 0 command. exit with 0')
		return 0

	variables = vars()
	del variables['cmd_list']
	logger.info( 'VARS: {}'.format(variables))

	close_cmp = False
	# file name
	xmod = 'a' if cont else 'w'
	if completed is not None and not isinstance(completed, file):
		completed = open(completed, xmod)
		close_cmp = True
	ntry = 0
	out_path0 = out_path
	tc_tasks0 = tc_tasks
	while True:
		ntry += 1
		logger.info('running {} commands: try {}'.format(len(cmd_list), ntry))
#		if out_path is not None:
#			out_path = '{}.{}'.format(out_path0, ntry)
		if mode == 'grid':
			avail_tasks = [tc_tasks0, len(cmd_list)]
			tc_tasks = min(avail_tasks)
			logger.info('reset tc_tasks to {} by {}'.format(tc_tasks, avail_tasks))
			job = Grid(cmd_list=cmd_list, tc_tasks=tc_tasks, grid_opts=grid_opts, 
						script=script, out_path=out_path, cpu=cpu, mem=mem, **kargs)
			logger.info('submiting jobs with {}'.format(job.grid))
			job_status = job.submit()
			exit_codes = [status for (task_id, status) in job_status]
		elif mode == 'local':
			avail_tasks = [tc_tasks0, avail_cpu(cpu), avail_mem(mem), len(cmd_list)]
			tc_tasks = min(avail_tasks)
			logger.info('reset tc_tasks to {} by {}'.format(tc_tasks, avail_tasks))
			job_status = pp_run(cmd_list, processors=tc_tasks)
			exit_codes = []
			fout = open(out_path, xmod) if out_path is not None else None
			for (stdout, stderr, status) in job_status:
				if fout is not None:
					print >>fout, '>>STATUS:\t{}\n>>STDOUT:\n{}\n>>STDERR:\n{}'.format(status, stdout, stderr)
				exit_codes += [status]
			if fout is not None:
				fout.close()
		uncmp = []
		for cmd, status in zip(cmd_list, exit_codes):
			if status > 0:
				uncmp += [cmd]
			elif completed is not None:
				completed.write(cmd + cmd_sep)
		if ntry >= retry or not uncmp:
			logger.info('finished with {} commands uncompleted'.format(len(uncmp)))
			break
		cmd_list = uncmp
	# close
	if close_cmp:
		completed.close()
	# if completed?
	return len(uncmp)
def avail_cpu(cpu):
	import psutil
	cpu_count = psutil.cpu_count()
	return max(1, int(1.0*cpu_count/cpu))
def avail_mem(mem):
	import psutil
	memory = psutil.virtual_memory()
	mem_free = memory.available
	mem = mem2float(mem)
	return max(1, int(1.0*mem_free/mem))
def mem2float(mem):
	import re
	d_mem = {'':1e1, 'k':1e3, 'm':1e6, 'g':1e9, 't':1e12}
	try:
		num, unit = re.compile(r'(\d+\.?\d*)([kmgt]?)', re.I).match(mem).groups()
		return float(num) * d_mem[unit.lower()]
	except AttributeError:
		raise AttributeError('Illegal MEMORY string `{}` (legal: 2g, 100m, 0.3t).'.format(mem))
def _qsub_log(jid, pwd, cmd):
    wlog = '''LOGFILE=/share/sge/default/common/working_dirs
JID={}
PWD={}
CMD={}
DATE=`date +"%Y-%m-%d-%H-%M-%S"`
echo "$JID:$PWD:$CMD:$(whoami):$DATE" >> $LOGFILE
'''.format(jid, pwd, cmd)
    run_cmd(wlog)
		
def file2list(cmd_file, sep="\n"):
	if not '\n' in sep:
		sep += '\n'
	if not os.path.exists(cmd_file) or not os.path.getsize(cmd_file):
		cmd_list = []
	else:
		f = open(cmd_file, 'r')
		cmd_list = f.read().split(sep)
	return [cmd for cmd in cmd_list if cmd]

def run_cmd(cmd, logger=None, log=None):
	if log is True and logger is None:
		logger = LOGGER
	if logger is not None:
		logger.info('run CMD: `{}`'.format(cmd))
	job = subprocess.Popen(cmd,stdout=subprocess.PIPE,\
							stderr=subprocess.PIPE,shell=True)
	output = job.communicate()
	status = job.poll()
	if logger is not None and status > 0:
		 logger.warn("exit code {} for CMD `{}`: ".format(status, cmd))
		 logger.warn('\n\tSTDOUT:\n{0}\n\tSTDERR:\n{1}\n\n'.format(*output))
	return output + (status,)

def default_processors(actual=None):
	from multiprocessing import cpu_count
	available_cpus = cpu_count()
	if actual:
		if actual > available_cpus:
			return available_cpus
		else:
			return actual
	else:
		return available_cpus
def pp_run(cmd_list, processors='autodetect'):
	if processors is None:
		processors = 'autodetect'
	ppservers = ()
	job_server = pp.Server(processors, ppservers=ppservers)
	jobs = [job_server.submit(run_cmd, (cmd,), (), ('subprocess',)) for cmd in cmd_list]
	return [job() for job in jobs]
def get_cmd_list(cmd_file, cmd_cpd_file=None, cmd_sep="\n", cont=True):
	if not '\n' in cmd_sep:
		cmd_sep += '\n'
	if not os.path.exists(cmd_file):
		raise IOError('Commands file : %s does NOT exist.'% (cmd_file, ))
	if cmd_cpd_file is None:
		cmd_cpd_file = cmd_file + '.completed'
	cmd_list = file2list(cmd_file, cmd_sep)
	cmd_cpd_list = file2list(cmd_cpd_file, cmd_sep)
	logger.info('{} commands in {}, {} commands in {}'.format(
				len(cmd_list), cmd_file, len(cmd_cpd_list), cmd_cpd_file))
	if cont:
		cmd_uncpd_list = sorted(list(set(cmd_list)-set(cmd_cpd_list)), \
							key=lambda x: cmd_list.index(x))
		logger.info('continue to run {} commands'.format(len(cmd_uncpd_list)))
	else:
		cmd_uncpd_list = sorted(list(set(cmd_list)), \
							key=lambda x: cmd_list.index(x))
	return cmd_uncpd_list
def submit_pp(cmd_file, processors=None, cmd_sep="\n", cont=True):
	if not '\n' in cmd_sep:
		cmd_sep += '\n'
	if not os.path.exists(cmd_file):
		raise IOError('Commands file : %s does NOT exist.'% (cmd_file, ))
	cmd_cpd_file = cmd_file + '.completed'
	cmd_list = file2list(cmd_file, cmd_sep)
	cmd_cpd_list = file2list(cmd_cpd_file, cmd_sep)
	if cont:
		cmd_uncpd_list = sorted(list(set(cmd_list)-set(cmd_cpd_list)), \
							key=lambda x: cmd_list.index(x))
	else:
		cmd_uncpd_list = sorted(list(set(cmd_list)), \
								key=lambda x: cmd_list.index(x))
		# not continue, so NONE completed
		cmd_cpd_list = []
	# for Argument list too long
	for i,cmd in enumerate(cmd_uncpd_list):
		if len(cmd.split('\n')) > 100:
			son_cmd_file = '%s.%s.sh' % (cmd_file, i)
			with open(son_cmd_file, 'w') as f:
				print >>f, cmd
			cmd_uncpd_list[i] = 'sh %s' % (son_cmd_file,)

	print '''
	total commands:\t%s
	skipped commands:\t%s
	retained commands:\t%s
	''' % (len(set(cmd_list)), \
	len(set(cmd_list))-len(cmd_uncpd_list), \
	len(cmd_uncpd_list))

	if not processors:
		processors = default_processors(len(cmd_uncpd_list))
	# start pp
	ppservers = ()
	job_server = pp.Server(processors, ppservers=ppservers)
	jobs = [(cmd, job_server.submit(run_cmd, (cmd,), (), ('subprocess',))) \
			for cmd in cmd_uncpd_list]
	# recover stdout, stderr and exit status
	cmd_out_file, cmd_err_file, cmd_warn_file = \
	cmd_file + '.out', cmd_file+'.err', cmd_file+'.warning'
	if cont:
		i = len(set(cmd_list))-len(cmd_uncpd_list)
		mode = 'a'
	else:
		i = 0
		mode = 'w'
	f = open(cmd_cpd_file, mode)
	f_out = open(cmd_out_file, mode)
	f_err = open(cmd_err_file, mode)
	f_warn = open(cmd_warn_file, mode)
	for cmd, job in jobs:
		i += 1
		out, err, status = job()
		f_out.write('CMD_%s_STDOUT:\n' % i + out + cmd_sep)
		f_err.write('CMD_%s_STDERR:\n' % i + err + cmd_sep)
		if not status == 0:
			f_warn.write(cmd + cmd_sep)
		else:
			f.write(cmd + cmd_sep)
	f.close()
	f_out.close()
	f_err.close()
	f_warn.close()

	job_server.print_stats()

def main():
	usage = __doc__ + "\npython %prog [options] <commands.list>"
	parser = OptionParser(usage, version="%prog " + __version__)
	parser.add_option("-p","--processors", action="store",type="int",\
					dest="processors", default=None, \
					help="number of processors [default=all available]")
	parser.add_option("-s","--separation", action="store", type="string",\
					dest="separation", default='\n', \
					help='separation between two commands [default="\\n"]')
	parser.add_option("-c","--continue", action="store", type="int",\
					dest="to_be_continue", default=1, \
					help="continue [1] or not [0] [default=%default]")
	parser.add_option("-m","--mode", action="store", type="choice",\
					dest="mode", default='local', choices=['local', 'grid'], \
					help='run mode [default=%default]')
	parser.add_option("--retry", action="store", type="int",\
					dest="retry", default=1, \
					help="retry times [default=%default]")
	parser.add_option("--grid-opts", action="store", type="string",\
					dest="grid_opts", default='-tc {tc} -l h_vmem={mem} -pe mpi {cpu}', \
					help='grid options [default="%default"]')
	(options,args)=parser.parse_args()
	if not args:
		parser.print_help()
		sys.exit()
	cmd_file = args[0]
	processors = options.processors
	separation = options.separation
	mode = options.mode
	grid_opts = options.grid_opts
	retry = options.retry
	to_be_continue = options.to_be_continue
#	submit_pp(cmd_file, processors=processors, \
#				cmd_sep=separation, cont=to_be_continue)
	run_job(cmd_file, tc_tasks=processors, mode=mode, grid_opts=grid_opts,
                cont=to_be_continue, retry=retry, cmd_sep=separation)
def run_job(cmd_file, cmd_list=None, tc_tasks=8, mode='grid', grid_opts='', cont=1,
            ckpt=None, retry=1, out_path=None, cmd_sep='\n', **kargs):
	if cmd_list is not None:
		with open(cmd_file, 'w') as fp:
			for cmd in cmd_list:
				print >> fp, cmd
	if kargs.get('cpu') is None:
		kargs['cpu'] = 1
	if kargs.get('mem') is None:
		kargs['mem'] = '1g'
	ckpt = cmd_file + '.ok' if ckpt is None else ckpt
	if cont and os.path.exists(ckpt):
		logger.info('check point file `{}` exists, skipped'.format(ckpt))
		return 0
	script = cmd_file + '.sh' if mode == 'grid' else None
	out_path = cmd_file + '.out' if out_path is None else out_path
	cmd_cpd_file = cmd_file + '.completed'
	if not cont and os.path.exists(ckpt):
		os.remove(ckpt)
	if not cont and os.path.exists(cmd_cpd_file):
		os.remove(cmd_cpd_file)
	if not cont and os.path.exists(out_path):
		os.remove(out_path)
	cmd_list = get_cmd_list(cmd_file, cmd_cpd_file, cmd_sep=cmd_sep, cont=cont)
	exit = run_tasks(cmd_list, tc_tasks=tc_tasks, mode=mode, grid_opts=grid_opts, 
				retry=retry, script=script, out_path=out_path, cont=cont,
				completed=cmd_cpd_file, cmd_sep=cmd_sep, **kargs)
	if not exit == 0:
		raise ValueError('faild to run {}, see detail in {}'.format(cmd_file, out_path))
	else:
		os.mknod(ckpt)
	return exit

if __name__ == '__main__':
	main()
