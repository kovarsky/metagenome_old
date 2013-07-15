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
LETTERS = ['A', 'T', 'G', 'C']
#ALT_TRESH = 'new'
#SUM_TRESH = 2
#FRAC_TRESH = 0
#TRESHOLDS = '_'.join([str(ALT_TRESH), str(SUM_TRESH), str(FRAC_TRESH)])

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
		# print 'positions: %s' % positions
		min_pos = min(positions)
		# print 'min position: %s' % min_pos
		self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
		# print 'indices: %s' % self.indices
		# print '-------------'
		# for j in self.indices:
		# 	print 'line%d: %s' % (j, self.bunch[j].line)
		# print '-------------'
		maxes = {i: getmax(self.bunch[i].freqs) for i in self.indices}
		# print 'maxes: %s' % maxes
		procindices = [i for i in self.indices if (self.bunch[i].freqs[maxes[i]] >= treshold 
										and self.bunch[i].freqs[maxes[i]] >= meansdata[i])]
		# print 'procindices: %s' % procindices
		if len(procindices) >= 2:
			# all_freqs = [self.bunch[i].freqs for i in self.indices]
			# summary = sum(all_freqs)
			# sum_reads = sum(summary)
			# summary[getmax(summary)] = 0
			# alt_reads = sum(summary)
			# #frac = float(alt_reads) / sum_reads
			#print alt_reads, sum_reads, frac
			# condition = (alt_reads >= ALT_TRESH and sum_reads >= SUM_TRESH)
			# if (condition):
			self.line = ''.join([ LETTERS[maxes[i]] if i in procindices else '-' for i in range(self.refnum)]) + '\n'
			# for i, j in itertools.combinations(procindices, 2):
			# 	###Another treshold!!!
			# 	#if (check_sample_cond(self.bunch[i].freqs) and check_sample_cond(self.bunch[j].freqs)):
			# 	self.cross_cov[i][j] += 1
			# 	self.cross_distance[i][j] += (maxes[i] != maxes[j])
		# for i in procindices:
		# 	self.cov[i] += 1
		# 	self.cross_cov[i][self.refnum] += 1
		# 	self.cross_distance[i][self.refnum] += (maxes[i] != self.bunch[i].ref_index)
		for i in self.indices:
			#print 'freqs%d: %s' % (i, self.bunch[i].freqs)
			self.bunch[i].readline()
		return(len(procindices) >= 2)

	#def save(self,directory):

		#np.savetxt('%s/cross_diffs_%s.out' % (directory, suffix), self.cross_distance, delimiter='\t',fmt='%.1f')
		#np.savetxt('%s/cross_cov_%s.out' % (directory, suffix), self.cross_cov, delimiter='\t',fmt='%.1f')
		#np.savetxt('%s/cov_%s.out' % (directory, suffix), self.cov, delimiter='\t',fmt='%.1f')


def main():
	plistname = sys.argv[1]
	contlistname = sys.argv[2]
	plist = open(plistname, 'r')
	contlist = open(contlistname, 'r')
	pileupnames = [line.strip() for line in plist]
	directory = os.path.dirname(pileupnames[0])
	means = [float(line.strip()) for line in open('%s/means.list' % directory)]
	conts = [line.strip() for line in contlist]
	outf = open('%s/snp_chars_%s' % (directory, suffix),'w')
	pileup_bunch = BunchOfObjects(pileupnames)
	for cont in conts:
		flag = pileup_bunch.get_current_cont_indices(cont)
		while flag:
			# print '=================='
			if pileup_bunch.distance(means):
				outf.write(pileup_bunch.line)
			flag = pileup_bunch.get_current_cont_indices(cont)
	#pileup_bunch.save(directory)
	pileup_bunch.close()
	outf.close()

if __name__ == '__main__':
	main()