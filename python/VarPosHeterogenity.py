import sys
import os

def get_diversity(freqs):
	sum_squares = sum(x ** 2 for x in freqs)
	cover = sum(freqs)
	res = (cover ** 2 - sum_squares) / float(cover * (cover - 1))
	return(res)

TRESHOLDS = (0, 2, 5, 10)
sum_pos = [0 for tresh in TRESHOLDS]
sum_diversity = [0.0 for tresh in TRESHOLDS]

f = open(sys.argv[1], 'r')
directory = os.path.dirname(sys.argv[1])
outf = open('%s/diversity.tsv' % directory, 'a')
for line in f:
	splitted = line.strip().split()
	freqs = [int(freq) for freq in splitted[2:6]]
	summ = int(splitted[6])
	if summ <= 1:
		diversity = 0.0
	else:
		diversity = get_diversity(freqs)
	for i, tresh in enumerate(TRESHOLDS):
		if summ >= tresh:
			sum_diversity[i] += diversity
			sum_pos[i] += 1

sum_diversity = '\t'.join('%.3f' % el for el in sum_diversity)
sum_pos = '\t'.join(str(el) for el in sum_pos)
outf.write('cumulative\t%s\t%s\n' % (sum_pos, sum_diversity))
outf.close()
