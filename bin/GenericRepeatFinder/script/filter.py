import sys

if len(sys.argv) < 4:
	print("Filter out TIRs/TDRs containing tandem repeats.")
	print('Usage: python filter.py <input> <phobos_output> <output>\n')
	sys.exit(1)	

input1 = sys.argv[1]
input2 = sys.argv[2]
output = sys.argv[3]

fin1 = open(input1, 'r')
fin2 = open(input2, 'r')
fout = open(output, 'w')

id = ''
prev = ''
count = 0

for line in fin2:
	if line[0] == '>':
		prev = line.rstrip()[:-2]
		break
		
for line in fin2:
	if line[0] == '#':
		continue
	elif line[0] == '>':
		id = line.rstrip()[:-2]
		if id != prev:
			if count > 4:
				fin1.readline()
				fin1.readline()
			else:
				fout.write(fin1.readline())
				fout.write(fin1.readline())
			count = 0
			prev = id
	else:
		count += 1
if count > 4:
	fin1.readline()
	fin1.readline()
else:
	fout.write(fin1.readline())
	fout.write(fin1.readline())
fin1.close()
fin2.close()
fout.close()
