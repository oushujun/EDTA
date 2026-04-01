import sys
import os
import pyfastx
import sqlite3
import bitarray
import re
import argparse

def options():
	parser = argparse.ArgumentParser(description='Your program description')

	# Required arguments
	parser.add_argument('--genome', required=True, help='Path to the genome file')
	parser.add_argument('--test', required=True, help='Path to the test RepeatMasker.out file')
	parser.add_argument('--reference', required=True, help='Path to the reference RepeatMasker.out file')

	# Optional flags
	parser.add_argument('--include_unknown', action='store_true', help='Include reference sequences which have labels outside of the normal TE categories.')
	parser.add_argument('--include_Ns', action='store_true', help='Include non-ATCG characters as part of the genome. These are inherently non-maskable and you probably should not use this.')
	parser.add_argument('--extended_report', action='store_true', help='Produce an extended report quantifying per-test sequence, per-TE category performance.')
	parser.add_argument('--min_entries', type=int, default=0, help='For the confusion matrix, exclude superfamilies with fewer than this many entries in the reference annotation. Default: 0 (no filtering)')

	args = parser.parse_args()
	
	return parser, args

class lib_test:
	def __init__(self, genome_file, ref, test, remove_Ns = True, include_unknown = False, extend_report = False):
		self.gf = genome_file
		self.ref = ref
		self.test = test
		
		self.remove_Ns = remove_Ns
		self.N_count = 0
		
		self.include_unknown = include_unknown
				
		self.gen_size = None
		self.chrom_sizes = None
		
		self.reference_cov = None
		self.test_cov = None
		
		self.extend = extend_report
		
		'''
		#Original script categories
		my %category;
		$category{'ltr'}="'RLG\\|RLC\\|RLB\\|RLR\\|RLE\\|\\\\s+LTR\\|RLX\\|Gypsy\\|Copia'";
		$category{'nonltr'}="'SINE\\|LINE\\|Penelope\\|RIT\\|RIL\\|RST\\|RIX\\|RSX\\|nonLTR\\|\\\\s+YR'";
		$category{'line'}="'LINE\\|RIL\\|RIT\\|RIX\\|Penelope'";
		$category{'sine'}="'SINE\\|RST\\|RSX'";
		$category{'tir'}="'TIR\\|MITE\\|hAT\\|hAT-Ac\\|MULE\\|MLE\\|MuDR\\|Tourist\\|CACT\\|PILE\\|POLE\\|Stowaway\\|TcMar-Stowaway\\|PIF\\|Harbinger\\|Tc1\\|En-Spm\\|EnSpm\\|CMC-EnSpm\\|PiggyBac\\|Mirage\\|P-element\\|Transib\\|DTA\\|DTH\\|DTT\\|DTM\\|DTC\\|DTA\\|TIR\\|DTX\\|DTR\\|DTE\\|Merlin\\|DTP\\|DTB\\|polinton'";
		$category{'mite'}="MITE";
		$category{'helitron'}="'Helitron\\|DHH\\|DHX\\|helitron'";
		$category{'total'}="[0-9]"; #grep any line with numbers
		$category{'classified'}="'Unknown\\|unknown\\/unknow\\|repeat_region\\|Unspecified'"; #unknown TEs of all kind
		'''
		
		script_dir = os.path.dirname(os.path.abspath(__file__))
		ontology_file = os.path.join(script_dir, 'bin', 'TE_Sequence_Ontology.txt')
		self.load_ontology(ontology_file)

		self.category_labels = {
								'ltr':re.compile('|'.join(['RLG', 'RLC', 'RLB', 'RLR', 'RLE', 'LTR', 'RLX', 'Gypsy', 'Copia']), re.IGNORECASE),
								'nonltr':re.compile('|'.join(['SINE', 'LINE', 'Penelope', 'RIT', 'RIL', 'RST', 'RIX', 'RSX', 'nonLTR', 'YR']), re.IGNORECASE),
								'line':re.compile('|'.join(['LINE', 'RIL', 'RIT', 'RIX', 'Penelope']), re.IGNORECASE),
								'sine':re.compile('|'.join(['SINE', 'RST', 'RSX']), re.IGNORECASE),
								'tir':re.compile('|'.join(['TIR', 'hAT', 'hAT-Ac', 'MULE', 'MLE', 'MuDR', 'Tourist', 'CACT', 'PILE', 'POLE', 'Stowaway', 'TcMar-Stowaway', 'PIF', 'Harbinger', 'Tc1', 'En-Spm', 'EnSpm', 'CMC-EnSpm', 'PiggyBac', 'Mirage', 'P-element', 'Transib', 'DTA', 'DTH', 'DTT', 'DTM', 'DTC', 'DTA', 'DTX', 'DTR', 'DTE', 'Merlin', 'DTP', 'DTB', 'polinton']), re.IGNORECASE),
								'mite':re.compile('|'.join(['MITE']), re.IGNORECASE),
								'helitron':re.compile('|'.join(['Helitron', 'DHH', 'DHX', 'helitron']), re.IGNORECASE),
								'others':re.compile('|'.join(['unknown','unknow','repeat_region','Unspecified']), re.IGNORECASE)
								}

		
	def load_genome(self):
		if not os.path.exists(f'{self.gf}.fxi'):
			print(f'Building self.gf index {self.gf}.fxi...')
			
		fa = pyfastx.Fasta(self.gf, full_index=True)
		self.gen_size = fa.size
		
		if 'N' in fa.composition:
			self.N_count += fa.composition['N']
		if 'n' in fa.composition:
			self.N_count += fa.composition['n']
			
		conn = sqlite3.connect(f'{self.gf}.fxi')
		curs = conn.cursor()
		self.chrom_sizes = dict(curs.execute('SELECT chrom, slen FROM seq').fetchall())
		'''
		#We actually CAN'T remove the Ns here - it has to be done at the very end. How does it contribute...?
		#Remove Ns from the sequence size for each chromosome/seq
		if self.remove_Ns:
			ok_atcg = set([65, 67, 71, 84])
			chrom_trans = dict(curs.execute('SELECT ID, chrom FROM seq').fetchall())
			for seqid, ascii_value, count in curs.execute('SELECT seqid, abc, num FROM comp').fetchall():
				if seqid in chrom_trans:
					chrom_value = chrom_trans[seqid]
					if ascii_value not in ok_atcg:
						self.chrom_sizes[chrom_value] -= count
		'''		
		curs.close()
		conn.close()
		
	def iterate_repeatmasker(self, file):
		header = 'placeholder'
		with open(file) as fh:
			#Read through the .out header until the first instance of a digit is encountered; that's the first real record
			while not header[0].isdigit():
				header = fh.readline().strip()
				if len(header) == 0:
					header = 'placeholder'
			
			segs = header.split()
			chrom = segs[4]
			origin_sequence = segs[9]
			label = segs[10]
			start, end = int(segs[5])-1, int(segs[6])
			
			yield chrom, origin_sequence, label, start, end
			
			for line in fh:
				segs = line.strip().split()
				chrom = segs[4]
				origin_sequence = segs[9]
				label = segs[10]
				start, end = int(segs[5])-1, int(segs[6])

				yield chrom, origin_sequence, label, start, end

	def create_report_card(self, conf, te_type):
		tpc = conf['tp']
		fpc = conf['fp']
		tnc = conf['tn']
		fnc = conf['fn']
		
		report_format = '''TE categoy:\t{cat}
		
TP:\t{tpc}
FN:\t{fnc}
TN:\t{tnc}
FP:\t{fpc}

\t\ttest_true\ttest_false
ref_true\t{tpc}\t{fnc}
ref_false\t{fpc}\t{tnc}

Sensitivity:\t{sens}
Specificity:\t{spec}
Accuracy:\t{acc}
Precision:\t{prec}
FDR:\t{fdr}
F1 measure:\t{f1}

#Metrics\tsens\tspec\taccu\tprec\tFDR\tF1\tTP\tTN\tFP\tFN
{tst}.{cat}.lib.report\t{sens}\t{spec}\t{acc}\t{prec}\t{fdr}\t{f1}\t{tpc}\t{tnc}\t{fpc}\t{fnc}
'''
		
		accuracy = (tpc+tnc) / (tpc + tnc + fpc + fnc)
		
		if tpc > 0 or fnc > 0:
			sensitivity = tpc / (tpc+fnc)
		else:
			sensitivity = 0
			
		if fpc > 0 or tnc > 0:
			specificity = tnc / (fpc+tnc)
		else:
			specificity = 0

		if tpc > 0 or fpc > 0:
			precision = tpc / (tpc+fpc)
			false_discovery = fpc / (tpc+fpc)
		else:
			precision = 0
			false_discovery = 1
			
		if tpc > 0 or fpc > 0 or fnc > 0:
			eff_one = (2*tpc) / (2*tpc + fpc + fnc)
		else:
			eff_one = 0

		report_format = report_format.format(genome_file=self.gf,
											std = self.ref,
											tst=self.test,
											tpc = tpc,
											fnc = fnc,
											tnc= tnc,
											fpc = fpc,
											sens = sensitivity,
											spec = specificity,
											acc = accuracy,
											prec = precision,
											fdr = false_discovery,
											f1 = eff_one,
											cat = te_type)
		
		return report_format
		
	def parent_teacher_meeting(self, extended_report):
		header = '\t'.join(['sequence_id', 
						'number_of_matches', 
						'tp_all_cats', 
						'fp_all_cats',
						'tp_line',
						'fp_line',
						'tp_sine',
						'fp_sine',
						'tp_ltr',
						'fp_ltr',
						'tp_nonltr',
						'fp_nonltr',
						'tp_tir',
						'fp_tir',
						'tp_mite',
						'fp_mite',
						'tp_helitron',
						'fp_helitron',
						'tp_unknown',
						'fp_unknown'
						])
		report = []
		report.append(header)
		for c in extended_report:
			next_row = '\t'.join([c, 
								str(extended_report[c]['any']['ct']), 
								str(extended_report[c]['any']['tp']),
								str(extended_report[c]['any']['fp']),
								str(extended_report[c]['line']['tp']),
								str(extended_report[c]['line']['fp']),
								str(extended_report[c]['sine']['tp']),
								str(extended_report[c]['sine']['fp']),
								str(extended_report[c]['ltr']['tp']),
								str(extended_report[c]['ltr']['fp']),
								str(extended_report[c]['nonltr']['tp']),
								str(extended_report[c]['nonltr']['fp']),
								str(extended_report[c]['tir']['tp']),
								str(extended_report[c]['tir']['fp']),
								str(extended_report[c]['mite']['tp']),
								str(extended_report[c]['mite']['fp']),
								str(extended_report[c]['helitron']['tp']),
								str(extended_report[c]['helitron']['fp']),
								str(extended_report[c]['others']['tp']),
								str(extended_report[c]['others']['fp'])
								])
			report.append(next_row)
		
		report = '\n'.join(report)
		
		return(report)
		
	def collect_out_file(self, outfile):
		groups = {}
		for chrom, origin_sequence, label, start, end in self.iterate_repeatmasker(outfile):
			if chrom not in groups:
				groups[chrom] = {}
			included = False
			if label != 'Unspecified':
				for cat in self.category_labels:
					if cat == 'others':
						continue
					if re.search(self.category_labels[cat], label):
						included = True
						if cat not in groups[chrom]:
							groups[chrom][cat] = []
						groups[chrom][cat].append((origin_sequence, start, end,))

						#Special case to always add mites to tirs
						if cat == 'mite':
							if 'tir' not in groups[chrom]:
								groups[chrom]['tir'] = []
							groups[chrom]['tir'].append((origin_sequence, start, end,))

			if not included:
				if 'others' not in groups[chrom]:
					groups[chrom]['others'] = []
				groups[chrom]['others'].append((origin_sequence, start, end,))

		return groups
		
	def load_ontology(self, ontology_file):
		self.alias_to_superfamily = {}
		self.superfamily_display_order = []
		seen = set()

		with open(ontology_file) as fh:
			for line in fh:
				line = line.strip()
				if not line or line.startswith('#'):
					continue
				parts = line.split('\t')
				if len(parts) < 3:
					continue
				so_name = parts[0].strip()
				aliases = [a.strip() for a in parts[2].split(',')]

				if so_name not in seen:
					self.superfamily_display_order.append(so_name)
					seen.add(so_name)

				for alias in aliases:
					if alias:
						self.alias_to_superfamily[alias] = so_name

	def superfamily_short_name(self, so_name):
		special = {
			'helitron': 'Helitron',
			'LINE_element': 'LINE/unknown',
			'SINE_element': 'SINE/unknown',
			'non_LTR_retrotransposon': 'nonLTR/unknown',
			'LTR_retrotransposon': 'LTR/unknown',
			'terminal_inverted_repeat_element': 'TIR/unknown',
			'DNA_transposon': 'DNA/unknown',
			'retrotransposon': 'Retro/unknown',
			'repeat_fragment': 'Unknown',
			'repeat_region': 'Unknown',
			'low_complexity': 'Low_complexity',
			'satellite_DNA': 'Satellite',
			'centromeric_repeat': 'Centromeric',
			'rRNA_gene': 'rRNA',
			'rDNA_intergenic_spacer_element': 'rDNA/IGS',
			'ERTBV_retrotransposon': 'LTR/ERTBV',
			'pararetrovirus': 'LTR/Pararetrovirus',
			'Penelope_retrotransposon': 'Penelope',
			'polinton': 'TIR/Polinton',
			'MITE': 'TIR/MITE',
			'TRIM': 'LTR/TRIM',
			'LARD': 'LTR/LARD',
		}
		if so_name in special:
			return special[so_name]
		prefix_map = {
			'_LTR_retrotransposon': 'LTR/',
			'_TIR_transposon': 'TIR/',
			'_LINE_retrotransposon': 'LINE/',
			'_SINE_retrotransposon': 'SINE/',
			'_YR_transposon': 'YR/',
			'_YR_retrotransposon': 'YR/',
		}
		for suffix, prefix in prefix_map.items():
			if so_name.endswith(suffix):
				return prefix + so_name[:-len(suffix)]
		return so_name

	def get_superfamily(self, label):
		if label in self.alias_to_superfamily:
			return self.alias_to_superfamily[label]
		return label

	def collect_by_superfamily(self, outfile):
		groups = {}
		for chrom, origin_sequence, label, start, end in self.iterate_repeatmasker(outfile):
			sf = self.get_superfamily(label)
			if chrom not in groups:
				groups[chrom] = {}
			if sf not in groups[chrom]:
				groups[chrom][sf] = []
			groups[chrom][sf].append((start, end))
		return groups

	def confusion_matrix(self, min_entries=0):
		ref_sf = self.collect_by_superfamily(self.ref)
		tst_sf = self.collect_by_superfamily(self.test)

		# Count reference entries per superfamily for filtering
		ref_entry_counts = {}
		for c in ref_sf:
			for sf in ref_sf[c]:
				ref_entry_counts[sf] = ref_entry_counts.get(sf, 0) + len(ref_sf[c][sf])

		# Collect all superfamily names
		all_ref_sfs = set(ref_entry_counts.keys())
		all_tst_sfs = set()
		for c in tst_sf:
			all_tst_sfs.update(tst_sf[c].keys())

		# Apply min_entries filter to reference superfamilies
		if min_entries > 0:
			all_ref_sfs = {sf for sf in all_ref_sfs if ref_entry_counts[sf] >= min_entries}

		# Order by ontology file order, then unexpected labels alphabetically
		def ordered(sfs):
			known = [s for s in self.superfamily_display_order if s in sfs]
			extra = sorted(sfs - set(self.superfamily_display_order))
			return known + extra

		ref_sfs = ordered(all_ref_sfs)
		tst_sfs = ordered(all_tst_sfs)

		# Short display names
		ref_names = [self.superfamily_short_name(s) for s in ref_sfs]
		tst_names = [self.superfamily_short_name(s) for s in tst_sfs]

		# Initialize matrix
		matrix = {r: {t: 0 for t in tst_sfs} for r in ref_sfs}
		missed = {r: 0 for r in ref_sfs}
		novel = {t: 0 for t in tst_sfs}

		for c in self.chrom_sizes:
			chrom_size = self.chrom_sizes[c]

			# Build per-superfamily ref bitarrays
			ref_arrays = {}
			ref_any = bitarray.bitarray(chrom_size)
			ref_any.setall(0)
			for sf in ref_sfs:
				arr = bitarray.bitarray(chrom_size)
				arr.setall(0)
				if c in ref_sf and sf in ref_sf[c]:
					for start, end in ref_sf[c][sf]:
						arr[start:end] = True
				ref_arrays[sf] = arr
				ref_any |= arr

			# Build per-superfamily test bitarrays
			tst_arrays = {}
			tst_any = bitarray.bitarray(chrom_size)
			tst_any.setall(0)
			for sf in tst_sfs:
				arr = bitarray.bitarray(chrom_size)
				arr.setall(0)
				if c in tst_sf and sf in tst_sf[c]:
					for start, end in tst_sf[c][sf]:
						arr[start:end] = True
				tst_arrays[sf] = arr
				tst_any |= arr

			# Compute pairwise overlaps
			for r_sf in ref_sfs:
				if ref_arrays[r_sf].count() == 0:
					continue
				for t_sf in tst_sfs:
					overlap = (ref_arrays[r_sf] & tst_arrays[t_sf]).count()
					if overlap > 0:
						matrix[r_sf][t_sf] += overlap
				missed[r_sf] += (ref_arrays[r_sf] & ~tst_any).count()

			for t_sf in tst_sfs:
				if tst_arrays[t_sf].count() == 0:
					continue
				novel[t_sf] += (tst_arrays[t_sf] & ~ref_any).count()

		# Build all rows as lists of strings for column-width calculation
		col_headers = ['Ref\\Test'] + tst_names + ['Missed', 'Total(bp)']

		rows = []
		rows.append(col_headers)

		for r_sf, r_name in zip(ref_sfs, ref_names):
			vals = [matrix[r_sf][t] for t in tst_sfs] + [missed[r_sf]]
			total = sum(vals)
			rows.append([r_name] + [str(v) for v in vals] + [str(total)])

		# Novel row
		novel_vals = [novel[t] for t in tst_sfs] + [0]
		novel_total = sum(novel_vals)
		rows.append(['Novel'] + [str(v) for v in novel_vals] + [str(novel_total)])

		# Column totals
		col_totals = []
		for t in tst_sfs:
			col_totals.append(sum(matrix[r][t] for r in ref_sfs) + novel[t])
		col_totals.append(sum(missed[r] for r in ref_sfs))
		grand_total = sum(col_totals)
		rows.append(['Total(bp)'] + [str(v) for v in col_totals] + [str(grand_total)])

		# Compute column widths and right-align numbers
		num_cols = len(col_headers)
		col_widths = [0] * num_cols
		for row in rows:
			for i, cell in enumerate(row):
				col_widths[i] = max(col_widths[i], len(cell))

		lines = []
		for row in rows:
			parts = []
			for i, cell in enumerate(row):
				if i == 0:
					parts.append(cell.ljust(col_widths[i]))
				else:
					parts.append(cell.rjust(col_widths[i]))
			lines.append('  '.join(parts))

		report = '\n'.join(lines)

		with open(f'{self.test}.superfamily_confusion_matrix.tsv', 'w') as out:
			print(report, file=out)

		return report

	def prepare_resources(self):
		self.load_genome()
		self.reference_cov = self.collect_out_file(self.ref)
		self.test_cov = self.collect_out_file(self.test)	
		
	def operate(self):
		shared_chromosomes = set(self.reference_cov.keys()).intersection(set(self.test_cov.keys()))
		ref_chroms_only = set(self.reference_cov.keys()) - shared_chromosomes
		tst_chroms_only = set(self.test_cov.keys()) - shared_chromosomes

		#Set up a repository of outputs
		record_keeper = {}
		if self.extend:
			extended_report = {}
		else:
			extended_report = None

		for cat in self.category_labels:
			record_keeper[cat] = {'tp':0,'tn':0,'fp':0,'fn':0}
		record_keeper['any'] = {'tp':0,'tn':0,'fp':0,'fn':0}
		if self.include_unknown:
			record_keeper['others'] = {'tp':0,'tn':0,'fp':0,'fn':0}

		for c in shared_chromosomes:
			#Build the all-categories test bitarray for "any" metric and extended reporting
			tst_all = bitarray.bitarray(self.chrom_sizes[c])
			tst_all.setall(0)
			for tst_cat in self.test_cov[c]:
				for origin, start, end in self.test_cov[c][tst_cat]:
					tst_all[start:end] = True
					if self.extend:
						if origin not in extended_report:
							extended_report[origin] = {}
							for ecat in self.category_labels:
								extended_report[origin][ecat] = {'ct':0, 'tp':0, 'fp':0}
							if self.include_unknown:
								extended_report[origin]['others'] = {'ct':0, 'tp':0, 'fp':0}
							extended_report[origin]['any'] = {'ct':0, 'tp':0, 'fp':0}

			#Set up the all-sequences reference array for this chromosome
			all_ref = bitarray.bitarray(self.chrom_sizes[c])
			all_ref.setall(0)

			#Reusable bitarrays
			ref = bitarray.bitarray(self.chrom_sizes[c])
			tst_cat_arr = bitarray.bitarray(self.chrom_sizes[c])

			#Calculate per-category values, update all-sequences array
			for cat in self.category_labels:
				#Reset the reference bitarray
				ref.setall(0)
				if cat in self.reference_cov[c]:
					for origin, start, end in self.reference_cov[c][cat]:
						if cat != 'others' or self.include_unknown:
							all_ref[start:end] = True
						ref[start:end] = True

				if self.include_unknown and cat != 'others':
					if 'others' in self.reference_cov[c]:
						for origin, start, end in self.reference_cov[c]['others']:
							ref[start:end] = True

				#Build category-specific test bitarray
				tst_cat_arr.setall(0)
				if cat in self.test_cov[c]:
					for origin, start, end in self.test_cov[c][cat]:
						tst_cat_arr[start:end] = True
				if self.include_unknown and cat != 'others':
					if 'others' in self.test_cov[c]:
						for origin, start, end in self.test_cov[c]['others']:
							tst_cat_arr[start:end] = True

				if self.extend:
					for tst_cat in self.test_cov[c]:
						for origin, start, end in self.test_cov[c][tst_cat]:
							extended_report[origin][cat]['ct'] += 1
							tp = ref[start:end].count()
							fp = end-start-tp
							extended_report[origin][cat]['tp'] += tp
							extended_report[origin][cat]['fp'] += fp

				tp = (ref & tst_cat_arr).count()
				fp = (~ref & tst_cat_arr).count()
				fn = (ref & ~tst_cat_arr).count()
				tn = self.chrom_sizes[c] - tp - fp - fn

				record_keeper[cat]['tp'] += tp
				record_keeper[cat]['fp'] += fp
				record_keeper[cat]['tn'] += tn
				record_keeper[cat]['fn'] += fn

			#Calculate all-sequences values
			tp = (all_ref & tst_all).count()
			fp = (~all_ref & tst_all).count()
			fn = (all_ref & ~tst_all).count()
			tn = self.chrom_sizes[c] - tp - fp - fn

			record_keeper['any']['tp'] += tp
			record_keeper['any']['fp'] += fp
			record_keeper['any']['tn'] += tn
			record_keeper['any']['fn'] += fn

			if self.extend:
				for tst_cat in self.test_cov[c]:
					for origin, start, end in self.test_cov[c][tst_cat]:
						extended_report[origin]['any']['ct'] += 1
						tp = all_ref[start:end].count()
						fp = end-start-tp
						extended_report[origin]['any']['tp'] += tp
						extended_report[origin]['any']['fp'] += fp

		#A chromosome only contained recoveries for the test sequence; by definition these are all false positives
		for c in tst_chroms_only:
			#For "any": all test annotations are FP
			tst_all = bitarray.bitarray(self.chrom_sizes[c])
			tst_all.setall(0)
			for tst_cat in self.test_cov[c]:
				for origin, start, end in self.test_cov[c][tst_cat]:
					tst_all[start:end] = True

			fp_all = tst_all.count()
			tn_all = self.chrom_sizes[c] - fp_all
			record_keeper['any']['fp'] += fp_all
			record_keeper['any']['tn'] += tn_all

			#Per-category: only count category-specific test annotations as FP
			tst_cat_arr = bitarray.bitarray(self.chrom_sizes[c])
			for cat in self.category_labels:
				tst_cat_arr.setall(0)
				if cat in self.test_cov[c]:
					for origin, start, end in self.test_cov[c][cat]:
						tst_cat_arr[start:end] = True
				if self.include_unknown and cat != 'others':
					if 'others' in self.test_cov[c]:
						for origin, start, end in self.test_cov[c]['others']:
							tst_cat_arr[start:end] = True

				fp = tst_cat_arr.count()
				tn = self.chrom_sizes[c] - fp
				record_keeper[cat]['fp'] += fp
				record_keeper[cat]['tn'] += tn

			if self.extend:
				for tst_cat in self.test_cov[c]:
					for origin, start, end in self.test_cov[c][tst_cat]:
						if origin not in extended_report:
							extended_report[origin] = {}
							for ecat in self.category_labels:
								extended_report[origin][ecat] = {'ct':0, 'tp':0, 'fp':0}
							if self.include_unknown:
								extended_report[origin]['others'] = {'ct':0, 'tp':0, 'fp':0}
							extended_report[origin]['any'] = {'ct':0, 'tp':0, 'fp':0}
						extended_report[origin]['any']['ct'] += 1
						extended_report[origin]['any']['fp'] += end - start

		#A chromosome contained only recoveries for the reference sequence; by definition these are all false negatives
		for c in ref_chroms_only:
			all_ref = bitarray.bitarray(self.chrom_sizes[c])
			all_ref.setall(0)

			ref = bitarray.bitarray(self.chrom_sizes[c])
			for cat in self.category_labels:
				ref.setall(0)
				if cat in self.reference_cov[c]:
					for origin, start, end in self.reference_cov[c][cat]:
						if cat != 'others' or self.include_unknown:
							all_ref[start:end] = True
						ref[start:end] = True

				fn = ref.count()
				tn = self.chrom_sizes[c] - fn

				record_keeper[cat]['fn'] += fn
				record_keeper[cat]['tn'] += tn

			fn_all = all_ref.count()
			tn_all = self.chrom_sizes[c] - fn_all
			record_keeper['any']['fn'] += fn_all
			record_keeper['any']['tn'] += tn_all
				
		
		if extended_report is not None:
			rep = self.parent_teacher_meeting(extended_report)
			with open(f'{self.test}.sequence_assessment.report', 'w') as out:
				print(rep, file = out)
			
		with open(f'{self.test}.repeatmasker.lib.report', 'w') as out:
			print(f'Genome:\t{self.gf}', file = out)
			print(f'Standard annotation:\t{self.ref}', file = out)
			print(f'Testing annotation:\t{self.test}', file = out)
			print('', file = out)
			print('sens=TP/(TP+FN)', file = out)
			print('spec=TN/(FP+TN)', file = out)
			print('accu=(TP+TN)/(TP+TN+FP+FN)', file = out)
			print('prec=TP/(TP+FP)', file = out)
			print('FDR=1-prec=FP/(TP+FP)', file = out)
			print('F1=2TP/(2TP+FP+FN)', file = out)
			print('', file = out)
			print(f'Ns removed:\t{self.remove_Ns}', file = out)
			print(f'Unknowns included:\t{self.include_unknown}', file = out)
			print('', file = out)

		
			for te in record_keeper:
				if self.remove_Ns:
					record_keeper[te]['tn'] -= self.N_count
				rep = self.create_report_card(record_keeper[te], te)
				print(rep, file = out)
			
			
		return None
				
def main():
	p, a = options()

	mn = lib_test(a.genome, a.reference, a.test, remove_Ns = not a.include_Ns, include_unknown = a.include_unknown, extend_report = a.extended_report)
	mn.prepare_resources()

	mn.operate()
	mn.confusion_matrix(min_entries=a.min_entries)
	
	
if __name__ == "__main__":
	main()
	