import sys

if len(sys.argv) < 3:
	print("Covert the interspersed repeat output of GRF (i.e., \"interspersed_repeat.out\") to FASTA format.")
	print('Usage: python convert.py <input> <output>\n')
	sys.exit(1)

input = sys.argv[1]
output = sys.argv[2]

family = 0

fin = open(input, 'r')
fout = open(output, 'w')
for line in fin:
	if line[0] == '-':
		family += 1
	elif line[0] == '>':
		fout.write(line.rstrip() + ':' + str(family) + '\n')
	else:
		fout.write(line)
fin.close()
fout.close()
