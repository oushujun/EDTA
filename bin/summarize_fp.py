#!/usr/bin/env python3
# Summarize per-family false positives of a TE annotation against a reference
# Shujun Ou (shujun.ou.1@gmail.com)
# v0.1: 03/13/2026
#
# Compares a "test" RepeatMasker .out (the annotation under evaluation) against a
# "reference" .out (trusted / gold-standard) at base-pair resolution and reports, per
# test family, how much of what it annotates is NOT backed by a reference annotation of
# the same TE subclass -- i.e. the family's false-positive footprint.
#
# Method: build a per-chromosome, per-subclass bitmask of the reference coverage
# (subclasses assigned from TE_Sequence_Ontology.txt), then walk the test .out and, per
# family, tally the bp overlapping each reference subclass, the bp with no reference at
# all (no_ref), and the bp that do not match the same subclass (fp_bp).
#
# Usage:
#   python3 summarize_fp.py --genome genome.fa --test test.out --reference ref.out
#   # --ontology PATH is optional; TE_Sequence_Ontology.txt is auto-detected by default.
# Output: an aligned per-family table sorted by fp_bp, largest first.
# Requires: pyfastx, bitarray (plus stdlib sqlite3, argparse).
import os
import sys
import pyfastx
import sqlite3
import bitarray
import argparse
from collections import defaultdict

SUBCLASSES = ['LTR', 'LINE', 'SINE', 'TIR', 'Helitron', 'nonTE', 'Unknown']

def options():
	parser = argparse.ArgumentParser(
		description='Per-family false-positive accounting of a test TE annotation against '
					'a reference annotation, by subclass, at base-pair resolution')
	parser.add_argument('--genome', required=True, help='Genome FASTA; used only for chromosome sizes (pyfastx indexes it)')
	parser.add_argument('--test', required=True, help='RepeatMasker .out being evaluated')
	parser.add_argument('--reference', required=True, help='Trusted / gold-standard RepeatMasker .out to score against')
	parser.add_argument('--ontology', default=None, help='TE_Sequence_Ontology.txt (default: auto-detect from the script dir)')
	return parser.parse_args()

def load_genome(genome_file):
	if not os.path.exists(f'{genome_file}.fxi'):
		print(f'Building index {genome_file}.fxi...', file=sys.stderr, flush=True)
	pyfastx.Fasta(genome_file, full_index=True)
	conn = sqlite3.connect(f'{genome_file}.fxi')
	curs = conn.cursor()
	chrom_sizes = dict(curs.execute('SELECT chrom, slen FROM seq').fetchall())
	curs.close()
	conn.close()
	return chrom_sizes

def iterate_repeatmasker(file):
	"""Yield (chrom, family, class, start, end) for each annotation line of a RepeatMasker .out.

	Skips the header block (every line until the first that begins with a digit). Columns are
	positional: seg[4]=chrom, seg[5]=query start (1-based), seg[6]=query end,
	seg[9]=repeat/family name, seg[10]=class/family label. start is returned 0-based.
	"""
	header = 'placeholder'
	with open(file) as fh:
		while not header[0].isdigit():
			header = fh.readline().strip()
			if len(header) == 0:
				header = 'placeholder'

		segs = header.split()
		yield segs[4], segs[9], segs[10], int(segs[5])-1, int(segs[6])

		for line in fh:
			segs = line.strip().split()
			yield segs[4], segs[9], segs[10], int(segs[5])-1, int(segs[6])

def parse_ontology(ontology_file):
	"""Parse TE_Sequence_Ontology.txt and return alias -> subclass mapping."""
	# Map SO names / section keywords to subclasses
	so_to_subclass = {}

	# LTR keywords
	for kw in ['LTR_retrotransposon', 'Copia', 'Gypsy', 'Bel_Pao', 'TRIM', 'LARD',
			   'Retrovirus', 'Endogenous_Retrovirus', 'pararetrovirus', 'ERTBV',
			   'long_terminal_repeat', 'RR_tract', 'primer_binding_site',
			   'Ngaro', 'DIRS', 'Viper', 'YR_retrotransposon']:
		so_to_subclass[kw] = 'LTR'

	# LINE keywords
	for kw in ['LINE_element', 'L1_LINE', 'R2_LINE', 'Jockey_LINE', 'I_LINE',
			   'RTE_LINE', 'CR1_LINE', 'CRE_LINE', 'Deceiver_LINE', 'Inkcap_LINE',
			   'Tad1_LINE', 'L2_LINE', 'R1_LINE', 'R4_LINE', 'Crack_LINE',
			   'Vingi_LINE', 'Tx1_LINE', 'Rex_LINE', 'Proto2_LINE',
			   'Penelope']:
		so_to_subclass[kw] = 'LINE'

	# SINE keywords
	for kw in ['SINE_element', 'tRNA_SINE', '5S_SINE', '7SL_SINE',
			   'Alu_SINE', 'B2_SINE', 'B4_SINE', 'ID_SINE', 'MIR_SINE']:
		so_to_subclass[kw] = 'SINE'

	# TIR keywords
	for kw in ['terminal_inverted_repeat_element', 'MITE', 'CACTA', 'hAT',
			   'Mutator', 'PIF_Harbinger', 'Tc1_Mariner', 'P_TIR', 'piggyBac',
			   'polinton', 'Transib', 'Merlin', 'terminal_inverted_repeat',
			   'PILE', 'POLE', 'Sola', 'Ginger', 'Kolobok', 'Dada', 'IS3EU',
			   'Zator', 'KDZ', 'DNA_transposon', 'Crypton']:
		so_to_subclass[kw] = 'TIR'

	# Helitron
	so_to_subclass['helitron'] = 'Helitron'

	# nonTE
	for kw in ['centromeric_repeat', 'knob', 'satellite_DNA', 'telomeric_repeat',
			   'subtelomere', 'low_complexity', 'chloroplast_DNA', 'mitochondrial_DNA',
			   'rRNA_gene', 'rDNA_intergenic_spacer_element',
			   'cytosolic_2S_rRNA', 'cytosolic_5S_rRNA', 'cytosolic_5_8S_rRNA',
			   'cytosolic_16S_rRNA', 'cytosolic_18S_rRNA', 'cytosolic_23S_rRNA',
			   'cytosolic_25S_rRNA', 'cytosolic_28S_rRNA',
			   'rRNA_5_external', 'rRNA_3_external',
			   'rRNA_internal_transcribed_spacer1', 'rRNA_internal_transcribed_spacer2',
			   'snRNA', 'scRNA']:
		so_to_subclass[kw] = 'nonTE'

	# Unknown
	for kw in ['repeat_region', 'repeat_fragment', 'retrotransposon', 'non_LTR_retrotransposon']:
		so_to_subclass[kw] = 'Unknown'

	# Parse the file
	alias_map = {}
	with open(ontology_file) as fh:
		for line in fh:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			parts = line.split('\t')
			if len(parts) < 3:
				continue

			so_name = parts[0]

			# Find subclass for this SO entry
			subclass = None
			for kw, sc in so_to_subclass.items():
				if kw in so_name:
					subclass = sc
					break

			if subclass is None:
				subclass = 'Unknown'

			# Parse aliases from column 3
			aliases = [a.strip() for a in parts[2].split(',')]
			for alias in aliases:
				if alias:
					alias_map[alias] = subclass

	return alias_map

def classify_label(label, alias_map):
	"""Classify a RepeatMasker class/family label into a subclass."""
	if label in alias_map:
		return alias_map[label]

	# Try prefix matching (e.g. "LTR/Gypsy" -> check "LTR/Gypsy")
	prefix = label.split('/')[0] if '/' in label else label
	if prefix in alias_map:
		return alias_map[prefix]

	# Fallback heuristics
	pl = prefix.lower()
	if pl in ('ltr',):
		return 'LTR'
	if pl in ('line',):
		return 'LINE'
	if pl in ('sine', 'sine?'):
		return 'SINE'
	if pl in ('dna', 'dnaauto', 'dnanona', 'tir', 'mite'):
		return 'TIR'
	if pl in ('helitron', 'rc'):
		return 'Helitron'
	if pl in ('simple_repeat', 'low_complexity', 'satellite', 'rrna', 'rdna', 'snrna', 'scrna'):
		return 'nonTE'

	return 'Unknown'

def main():
	args = options()

	# Auto-detect ontology path
	if args.ontology is None:
		script_dir = os.path.dirname(os.path.abspath(__file__))
		args.ontology = os.path.join(script_dir, 'TE_Sequence_Ontology.txt')

	alias_map = parse_ontology(args.ontology)
	chrom_sizes = load_genome(args.genome)

	# Build per-subclass reference bitmasks per chromosome
	ref_masks = {}  # {chrom: {subclass: bitarray}}
	for chrom, origin, label, start, end in iterate_repeatmasker(args.reference):
		if chrom not in ref_masks:
			ref_masks[chrom] = {}
			for sc in SUBCLASSES:
				ref_masks[chrom][sc] = bitarray.bitarray(chrom_sizes[chrom])
				ref_masks[chrom][sc].setall(0)

		subclass = classify_label(label, alias_map)
		ref_masks[chrom][subclass][start:end] = True

	# Iterate test annotations and compute per-family reference composition
	family_stats = defaultdict(lambda: {
		'class': '', 'test_subclass': '', 'annotations': 0, 'total_bp': 0,
		**{sc: 0 for sc in SUBCLASSES}, 'no_ref': 0, 'fp_bp': 0
	})

	for chrom, origin, label, start, end in iterate_repeatmasker(args.test):
		entry = family_stats[origin]
		entry['class'] = label
		test_sc = classify_label(label, alias_map)
		entry['test_subclass'] = test_sc
		entry['annotations'] += 1
		span = end - start
		entry['total_bp'] += span

		if chrom in ref_masks:
			for sc in SUBCLASSES:
				overlap = ref_masks[chrom][sc][start:end].count()
				entry[sc] += overlap
			# no_ref = positions not covered by any reference annotation
			all_ref_count = bitarray.bitarray(end - start)
			all_ref_count.setall(0)
			for sc in SUBCLASSES:
				all_ref_count |= ref_masks[chrom][sc][start:end]
			ref_total = all_ref_count.count()
			entry['no_ref'] += span - ref_total
			# fp_bp = positions not matching the SAME subclass in reference
			if test_sc in ref_masks[chrom]:
				same_sc_overlap = ref_masks[chrom][test_sc][start:end].count()
			else:
				same_sc_overlap = 0
			entry['fp_bp'] += span - same_sc_overlap
		else:
			entry['no_ref'] += span
			entry['fp_bp'] += span

	# Compute column widths for aligned output
	headers = ['family', 'test_class', 'test_subclass', 'annotations', 'total_bp', 'fp_bp'] + \
			  ['ref_' + sc for sc in SUBCLASSES] + ['no_ref']

	rows = []
	for family, stats in sorted(family_stats.items(), key=lambda x: x[1]['fp_bp'], reverse=True):
		row = [
			family,
			stats['class'],
			stats['test_subclass'],
			str(stats['annotations']),
			str(stats['total_bp']),
			str(stats['fp_bp']),
		]
		for sc in SUBCLASSES:
			row.append(str(stats[sc]))
		row.append(str(stats['no_ref']))
		rows.append(row)

	# Calculate max width per column
	widths = [len(h) for h in headers]
	for row in rows:
		for i, val in enumerate(row):
			widths[i] = max(widths[i], len(val))

	# Print header
	header_line = '  '.join(h.rjust(widths[i]) if i >= 3 else h.ljust(widths[i]) for i, h in enumerate(headers))
	print(header_line)

	# Print rows
	for row in rows:
		line = '  '.join(val.rjust(widths[i]) if i >= 3 else val.ljust(widths[i]) for i, val in enumerate(row))
		print(line)

if __name__ == "__main__":
	main()
