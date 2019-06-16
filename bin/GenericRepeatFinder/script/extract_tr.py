import sys

def tir_len(tir):
	l1 = l2 = 0
	d = 0
	for c in tir:
		if c == 'm' or c == 'M':
			l1 += d
			l2 += d
			d = 0
		elif c == 'I':
			l1 += d
			d = 0
		elif c == 'D':
			l2 += d
			d = 0
		else:
			d = d * 10 + int(c)
	return (l1, l2)
	
if len(sys.argv) < 3:
	print("Extract only the TRs (without spacer) of TIR/TDR output of GRF.")
	print('Usage: python extract_tr.py <input> <output>\n')
	sys.exit(1)	

input = sys.argv[1]
output = sys.argv[2]

fin = open(input, 'r')
fout = open(output, 'w')

id = ''
for line in fin:
	if line[0] == '>':
		id = line.rstrip()		
	else:
		lines = id.split(':')
		l1, l2 = tir_len(lines[3])
		seq = line.rstrip()
		fout.write(id + ':l\n')
		fout.write(seq[0:l1] + '\n')
		fout.write(id + ':r\n')
		fout.write(seq[(-l2):] + '\n')
fin.close()
fout.close()
