#!/usr/bin/python

# this program extract inverted repeats detected by Generic Repeat Finder (GRF) that are high-quality, which is defined as
#   ration between the spacer length and one arm (stem) length is <= 0.5. that is equivalent of that spacer is no more 
#   than 1/5th of total inverted repeat length

from __future__ import division
import sys

def getStem(tir):      # get the length of both stems
	num = 0
	d = 0
	for c in tir:  # c is each character in tir
		if c == 'm' or c == 'M':
			num += 2 * d
			d = 0
		elif c == 'I' or c == 'D':
			num += d
			d = 0
		else:
			d = d * 10 + int(c)
	return num

input = sys.argv[1]
ratio = float(sys.argv[2])
output = sys.argv[3]

flag = False
fin = open(input, 'r')
fout = open(output, 'w')
for line in fin:	
	if line[0] == '>':
                # print (line)
		lines = line.rstrip().split(':')
		l = int(lines[2]) - int(lines[1]) + 1
		spacer = l - getStem(lines[3])
                # l-spacer is of length of both stems
		if spacer / (l - spacer) * 2 <= ratio:
			fout.write(line)
			flag = True
		else:
			flag = False
	else:
		if flag:
			fout.write(line)
		
fin.close()
fout.close()
