import sys

f = open(sys.argv[1])
n = int(sys.argv[2])
summ = 0
i = 0
for i, line in enumerate(f):
	num = int(line.strip().split()[n])
	summ += num
print '%d\t%d' % (i + 1, summ)
f.close()
