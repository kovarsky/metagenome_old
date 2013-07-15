import sys
import random
import os
import numpy as np
import re
import itertools
import CrossHeterogenity


plistname = sys.argv[1]

LETTERS = ('A', 'T', 'G', 'C')

def getmax(vector):
    """Shows 'letter' with max coverage"""
    maxcov = max(vector)
    ind_of_max = [i for i in xrange(4) if vector[i] == maxcov]
    ind = random.choice(ind_of_max)
    return(ind)

class PileupObj():
	"""docstring for PileupObj"""
	def __init__(self, pileupname):
		self.data = open(pileupname, 'r')
		self.readline()
		self.not_eof = True
		# self.readflag = True
		# self.contigflag = False

	def readline(self):
		line = self.data.readline()
		if line != '':
			self.line = line
			self.seq, self.coord, ref, self.summ, reads, quals = line.split('\t')
			self.coord = int(self.coord)
			reads = re.sub('[\^].|\$', '', reads)
			self.freqs = np.array([reads.upper().count(letter) for letter in LETTERS])
			ref_index = LETTERS.index(ref)
			self.ref_index = ref_index
			self.freqs[ref_index] = int(self.summ) - sum(self.freqs)
			self.summ = int(self.summ)
			self.res = (getmax(self.freqs) != self.ref_index)
			self.line = line.strip()
		else:
			self.seq = None
			self.not_eof = False

	def close(self):
		self.data.close()



def main():
	pnames = (line.strip() for line in open(plistname, 'r'))
	for pname in pnames:
		pilefile = PileupObj(pname)
		directory = os.path.dirname(pname)
		outname = (os.path.splitext(pname)[0] + '.refsnphist')
		counts = {}
		pilefile.readline()
		while pilefile.not_eof:
			if pilefile.res:
				#print(pilefile.freqs)
				#print (pilefile.line)
				#print('-----------------')
				counts[pilefile.summ] = counts.setdefault(pilefile.summ, 0) + 1
			pilefile.readline()
		outf = open(outname, 'w')
		for item in counts.items():
			outf.write('%d\t%d\n' % item)
		outf.close()
		pilefile.close()




if __name__ == '__main__':
	main()
