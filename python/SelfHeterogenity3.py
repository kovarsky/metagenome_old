import sys
import re
import os

def get_diversity(freqs):
	sum_squares = sum(x ** 2 for x in freqs)
	cover = sum(freqs)
	res = (cover ** 2 - sum_squares) / float(cover * (cover - 1))
	return(res)

LETTERS = ('A', 'T', 'G', 'C')
TRESHOLDS = (10,20,25,30)

def pileup_process(pileupname):	

	sum_pos = [0 for tresh in TRESHOLDS]
	sum_diversity = [0.0 for tresh in TRESHOLDS]
	pileup = open(pileupname, 'r')

	for line in pileup:
		seq, coord, ref, summ, reads, quals = line.split('\t')
		summ = int(summ)
		reads = re.sub('[\^].|\$', '', reads)
		if (len(reads) != summ):
			print 'Unequal lengths'
			print len(reads), summ
			print line
			print reads
			sys.exit()
		freqs = {letter: reads.count(letter) for letter in LETTERS}
		alt_len = sum(freqs.values())
		#I can add here summ_treshold
		if any((el >= 2 and (el / float(summ)) > 0.1)for el in freqs.values()):
			freqs[ref] = summ - alt_len
			diversity = get_diversity(freqs.values())
		else:
			diversity = 0.0
			#freqs[ref] = summ - alt_len
			#diversity = get_diversity(freqs.values())
		#print freqs
		#print line
		#raw_input("Press Enter to continue...")
		for i, tresh in enumerate(TRESHOLDS):
			if summ >= tresh:
				sum_diversity[i] += diversity
				sum_pos[i] += 1
	pileup.close()
	sum_diversity = '\t'.join('%.3f' % el for el in sum_diversity)
	sum_pos = '\t'.join(str(el) for el in sum_pos)
	return(sum_pos, sum_diversity)

pileup_list = open(sys.argv[1], 'r')
flag = sys.argv[2]
for i, pileupname in enumerate(pileup_list):
	pileupname = pileupname.strip()
	sum_pos, sum_diversity = pileup_process(pileupname)
	if i == 0:
		directory = os.path.dirname(pileupname)
		outf = open('%s/diversity_%s.tsv' % (directory, flag), 'w')
	pileupname_short = re.sub('\.pileup', '', os.path.basename(pileupname))
	outf.write('%s\t%s\t%s\n' % (pileupname_short, sum_pos, sum_diversity))
outf.close()
pileup_list.close()