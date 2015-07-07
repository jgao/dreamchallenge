import random

infile = open('staph.frag_1.fastq')
outfile = open('staph.frag_1.fastq.rand', 'w')

tag = ['A', 'T', 'C', 'G']

i=0;

for line in infile:
	if i%2500 == 0:
		print i
	l = list(line)
	l = map(lambda x: tag[random.randint(0,3)] if x == 'N' else x, l)
	outfile.write("".join(l))
	i += 1

infile.close()
outfile.close()

print "done";
'''

infile = open('/run/media/janet/2CDAC4ABDAC4731E/Jshlorv/IS1data/SRR1718501.fastq')
outfile = open('/run/media/janet/2CDAC4ABDAC4731E/Jshlorv/IS1data/SRR1718501.fastq.rand', 'w')

tag = ['A', 'T', 'C', 'G']

i=0;

for line in infile:
	if i%2500 == 0:
		print i
	l = list(line)
	l = map(lambda x: tag[random.randint(0,3)] if x == 'N' else x, l)
	outfile.write("".join(l))
	i += 1

infile.close()
outfile.close()

print "done";'''