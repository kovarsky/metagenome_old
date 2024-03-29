import sys
import random
import os
import numpy as np
import re
import itertools
import CrossHeterogenity


def getmax(vector):
    """Shows 'letter' with max coverage"""
    maxcov = max(vector)
    ind_of_max = [i for i in xrange(4) if vector[i] == maxcov]
    ind = random.choice(ind_of_max)
    return(ind)

def check_sample_cond(vector):
	summ_cond = (sum(vector) >= 2)
	sorted_vector = sorted(vector)
	if summ_cond == 2:
		alt_cond = (sorted_vector[2] != 1)
	else:
		alt_cond = True
	return(alt_cond and summ_cond)

suffix = sys.argv[3]
treshold = int(sys.argv[4])


class BunchOfObjects(CrossHeterogenity.BunchOfObjects):
	"""docstring for BunchOfObjects"""
	def __init__(self, pileupnames):
		self.bunch = [CrossHeterogenity.PileupObj(pileupname) for pileupname in pileupnames]
		length = len(pileupnames)
		self.cross_cov = np.zeros((length + 1, length + 1))
		self.cross_distance = np.zeros((length + 1, length + 1))
		self.cov = np.zeros((length))
		self.refnum = length

	def get_current_cont_indices(self, current_cont):
		self.indices = [i for i, obj in enumerate(self.bunch) if (obj.seq == current_cont and obj.not_eof)]
		if len(self.indices) < 1:
			self.indices = range(len(self.bunch))
			return False
		else:
			return True

	def distance(self, meansdata):
		positions = [self.bunch[i].coord for i in self.indices]
		min_pos = min(positions)
		self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
		maxes = {i: getmax(self.bunch[i].freqs) for i in self.indices}
		procindices = [i for i in self.indices if (self.bunch[i].freqs[maxes[i]] >= treshold 
										and self.bunch[i].freqs[maxes[i]] >= meansdata[i])]
		if len(procindices) >= 2:
			for i, j in itertools.combinations(procindices, 2):
				self.cross_cov[i][j] += 1
				self.cross_distance[i][j] += (maxes[i] != maxes[j])
		for i in procindices:
			self.cov[i] += 1
			self.cross_cov[i][self.refnum] += 1
			self.cross_distance[i][self.refnum] += (maxes[i] != self.bunch[i].ref_index)
		for i in self.indices:
			self.bunch[i].readline()

	def save(self,directory):
		np.savetxt('%s/cross_diffs_%s.out' % (directory, suffix), self.cross_distance, delimiter='\t',fmt='%.1f')
		np.savetxt('%s/cross_cov_%s.out' % (directory, suffix), self.cross_cov, delimiter='\t',fmt='%.1f')
		np.savetxt('%s/cov_%s.out' % (directory, suffix), self.cov, delimiter='\t',fmt='%.1f')


def main():
	plistname = sys.argv[1]
	contlistname = sys.argv[2]
	plist = open(plistname, 'r')
	contlist = open(contlistname, 'r')
	pileupnames = [line.strip() for line in plist]
	directory = os.path.dirname(pileupnames[0])
	means = [float(line.strip()) for line in open('%s/means.list' % directory)]
	conts = [line.strip() for line in contlist]
	pileup_bunch = BunchOfObjects(pileupnames)
	for cont in conts:
		flag = pileup_bunch.get_current_cont_indices(cont)
		while flag:
			pileup_bunch.distance(means)
			flag = pileup_bunch.get_current_cont_indices(cont)
	pileup_bunch.save(directory)
	pileup_bunch.close()

if __name__ == '__main__':
	main()