import sys

f = open(sys.argv[1],'r')
sums = {}
for line in f:
	num = int(line.strip().split()[2])
	sums[num] = sums.setdefault(num, 0) + 1
f.close()
outf = open(sys.argv[2],'w')
for item in sums.iteritems():
	outf.write('%d\t%d\n' % item)
outf.close()