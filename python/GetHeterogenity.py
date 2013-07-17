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
occtreshold = int(sys.argv[5])

class PileupObj(CrossHeterogenity.PileupObj):



class BunchOfObjects(CrossHeterogenity.BunchOfObjects):
	"""docstring for BunchOfObjects"""
	def __init__(self, pileupnames):
		self.bunch = [CrossHeterogenity.PileupObj(pileupname) for pileupname in pileupnames]
		length = len(pileupnames)
		self.cross_cov = np.zeros((length + 1, length + 1))
		self.cross_distance = np.zeros((length + 1, length + 1))
		self.cov = np.zeros((length))
		self.refnum = length
		self.index = 0

	def get_current_cont_indices(self, current_cont):
		self.indices = [i for i, obj in enumerate(self.bunch) if (obj.seq == current_cont and obj.not_eof)]
		if len(self.indices) < 1:
			self.indices = range(len(self.bunch))
			return False
		else:
			return True

	def get_alleles(self, meansdata):
		positions = [self.bunch[i].coord for i in self.indices]
		min_pos = min(positions)
		self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
		maxes = {i: getmax(self.bunch[i].freqs) for i in self.indices}
		procindices = [i for i in self.indices if (self.bunch[i].freqs[maxes[i]] >= treshold)] 
		#								and self.bunch[i].freqs[maxes[i]] >= meansdata[i])]
		bins = np.bincount(maxes.values(), minlength = 4)
		counts = sorted(bins)
		major = getmax(bins)
		condition = (counts[2] >= occtreshold and len(procindices) >= 2)
		if condition:
			self.lines = {}
			for i in procindices:
				eq_reads = str(self.bunch[i].freqs[major])
				ind = getmax(self.bunch[i].freqs)
				max_reads = str(self.bunch[i].freqs[ind])
				sum_reads = str(self.bunch[i].summ)
				self.lines[i] = '%d\t%s\n' % (self.index, '\t'.join([str(inds), seq_reads, max_reads, sum_reads]))
			self.variant_freqs = '%d\t%s\n' % (self.index, '\t'.join([str(el) for el in bins]))
			self.index += 1
		for i in self.indices:
			self.bunch[i].readline()

		return condition

	# def save(self,directory):
	# 	np.savetxt('%s/cross_diffs_%s.out' % (directory, suffix), self.cross_distance, delimiter='\t',fmt='%.1f')
	# 	np.savetxt('%s/cross_cov_%s.out' % (directory, suffix), self.cross_cov, delimiter='\t',fmt='%.1f')
	# 	np.savetxt('%s/cov_%s.out' % (directory, suffix), self.cov, delimiter='\t',fmt='%.1f')


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
	outfs = [open((os.path.splitext(pname)[0] + '.het'), 'w') for pname in pileupnames]
	variant_freqs_file = open('%s/variant_freqs_%s.tsv' % (directory, suffix), 'w')
	for cont in conts:
		flag = pileup_bunch.get_current_cont_indices(cont)
		while flag:
			if pileup_bunch.get_alleles(means):
				for num, line in pileup_bunch.lines.items():
					outfs[num].write(line)
				variant_freqs_file.write(pileup_bunch.variant_freqs)
			flag = pileup_bunch.get_current_cont_indices(cont)
	pileup_bunch.close()

if __name__ == '__main__':
	main()