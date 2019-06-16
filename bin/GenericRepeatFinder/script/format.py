import sys

if len(sys.argv) < 3:
	print("Format the filtered interspersed repeat FASTA file to original GRF format.")
	print('Usage: python format.py <input> <output>\n')
	sys.exit(1)

input = sys.argv[1]
output = sys.argv[2]

prev = ''

fin = open(input, 'r')
fout = open(output, 'w')
for line in fin:
	if line[0] == '>':		
		lines = line.rstrip().split(':')
		id = ':'.join(lines[0:4])
		if lines[4] != prev:
			fout.write('--------------------------------------------------\n')		
			prev = lines[4]
		fout.write(id + '\n')
	else:
		fout.write(line)
fin.close()
fout.close()
