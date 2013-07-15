import sys
import os
import numpy as np
import re
import itertools
#import CrossHeterogenity

LETTERS = ('A', 'T', 'G', 'C')
class PileupObj():
	"""docstring for PileupObj"""
	def __init__(self, pileupname):
		self.data = open(pileupname, 'r')
		self.counts = np.zeros(20)
		self.readline()
		self.not_eof = True
		# self.readflag = True
		# self.contigflag = False

	def readline(self):
		line = self.data.readline()
		if line != '':
			self.seq, self.coord, ref, self.summ, reads, quals = line.split('\t')
			self.coord = int(self.coord)
			reads = re.sub('[\^].|\$', '', reads)
			self.freqs = np.array([reads.count(letter) for letter in LETTERS])
			ref_index = LETTERS.index(ref)
			self.summ = int(self.summ)
			self.freqs[ref_index] = int(self.summ) - sum(self.freqs)
			#self.summ = float(self.summ)
		else:
			self.not_eof = False

	def get_probs(self):
		sorted_freqs = sorted(self.freqs)
		frac = sorted_freqs[2] / float(self.summ)
		ind = int(frac // 0.05)
		return(sorted_freqs[2], ind)



	def close(self):
		self.data.close()


def main():
	plistfile = open(sys.argv[1], 'r')
	pileuplist = [line.strip() for line in plistfile]
	directory = os.path.dirname(pileuplist[0])
	res = np.zeros((len(pileuplist), 20))
	for i, pileup in enumerate(pileuplist):
		data = PileupObj(pileup)
		while data.not_eof:
			if data.summ >= 10:
				max_alt, ind = data.get_probs()
				res[i][19] += 1
				if max_alt >= 2:
					res[i][ind] += 1
			data.readline()
	np.savetxt('%s/alt_spectra.out' % directory, res, delimiter='\t',fmt='%d')

if __name__ == '__main__':
	main()



		
