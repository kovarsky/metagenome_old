import sys
import os
import numpy as np
import re
import itertools

def get_diversity(freqs):
	sum_squares = np.dot(freqs, freqs)
	cover = sum(freqs)
	res = (cover ** 2 - sum_squares) / float(cover * (cover - 1))
	return res

LETTERS = ('A', 'T', 'G', 'C')
TRESHOLD = 5

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
			self.seq, self.coord, ref, self.summ, reads, quals = line.split('\t')
			self.coord = int(self.coord)
			reads = re.sub('[\^].|\$', '', reads)
			self.freqs = np.array([reads.upper().count(letter) for letter in LETTERS])
			ref_index = LETTERS.index(ref)
			self.ref_index = ref_index
			self.freqs[ref_index] = int(self.summ) - sum(self.freqs)
			self.summ = int(self.summ)
			self.line = line.strip()
		else:
			self.seq = None
			self.not_eof = False

	def close(self):
		self.data.close()


class BunchOfObjects(object):
	"""docstring for BunchOfObjects"""
	def __init__(self, pileupnames):
		self.bunch = [PileupObj(pileupname) for pileupname in pileupnames]
		length = len(pileupnames)
		self.cross_cov = np.zeros((length, length))
		self.cross_hetero = np.zeros((length, length))

	def get_current_cont_indices(self, current_cont):
		self.indices = [i for i, obj in enumerate(self.bunch) if (obj.seq == current_cont and obj.not_eof)]
		if len(self.indices) < 2:
			self.indices = range(len(self.bunch))
			return False
		else:
			return True
		
	def comparison(self):
		positions = [self.bunch[i].coord for i in self.indices]
		min_pos = min(positions)
		self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
		if len(self.indices) >= 2:
			for i, j in itertools.combinations(self.indices, 2):
				s1 = float(sum(self.bunch[i].freqs))
				s2 = float(sum(self.bunch[j].freqs))
				summ = s1 + s2
				freqs = (self.bunch[i].freqs / s1 + self.bunch[j].freqs / s2) * summ * 0.5
				if summ >= TRESHOLD:
					self.cross_cov[i][j] += 1
					self.cross_hetero[i][j] += get_diversity(freqs)
		for i in self.indices:
			self.bunch[i].readline()

	def close(self):
		for obj in self.bunch:
			obj.close()


def main():
	plistname = sys.argv[1]
	contlistname = sys.argv[2]
	plist = open(plistname, 'r')
	contlist = open(contlistname, 'r')
	pileupnames = [line.strip() for line in plist]
	directory = os.path.dirname(pileupnames[0])
	conts = [line.strip() for line in contlist]
	pileup_bunch = BunchOfObjects(pileupnames)
	i = 1
	for cont in conts:
		flag = pileup_bunch.get_current_cont_indices(cont)
		while flag:
			pileup_bunch.comparison()
			flag = pileup_bunch.get_current_cont_indices(cont)
	np.savetxt('%s/cross_cov_%d.out' % (directory, TRESHOLD), pileup_bunch.cross_cov,delimiter='\t',fmt='%.1f')
	np.savetxt('%s/cross_hetero_%d.out' % (directory, TRESHOLD), pileup_bunch.cross_hetero,delimiter='\t',fmt='%.3f')
	pileup_bunch.close()


if __name__ == '__main__':
	main()



		
