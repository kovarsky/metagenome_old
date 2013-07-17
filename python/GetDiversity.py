import sys
import random
import os
import numpy as np
import re
import itertools
import CrossHeterogenity


def get_diversity(freqs):
	sum_squares = sum(x ** 2 for x in freqs)
	cover = sum(freqs)
	res = (cover ** 2 - sum_squares) / float(cover * (cover - 1))
	if res < 0:
		print res
		print cover
		print freqs
		sys.exit()
	return(res)

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
		self.diversities = [{} for pileupname in pileupnames]
		self.diversity = {}
		self.norm_diversity = {}
		self.covers = {}
		self.hetero_pos = [{} for pileupname in pileupnames]
		#self.lines = {}

	def distance(self, meansdata):
		positions = [self.bunch[i].coord for i in self.indices]
		# print 'positions: %s' % positions
		min_pos = min(positions)
		# print 'min position: %s' % min_pos
		self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
		# print 'indices: %s' % self.indices
		# print '-------------'
		# for j in self.indices:
			# print 'line%d: %s' % (j, self.bunch[j].line)
		# print '-------------'
		maxes = {i: getmax(self.bunch[i].freqs) for i in self.indices}
		# print 'maxes: %s' % maxes
		procindices = [i for i in self.indices if (self.bunch[i].summ >= treshold)]
		# print procindices
		if len(procindices) > 0:
			sum_freqs = sum(self.bunch[i].freqs for i in procindices)
			# print 'sum_freqs: ', sum_freqs
			sum_diversity = get_diversity(sum_freqs)
			# print 'sum_diversity: ', sum_diversity
			# norm_sum_freqs = sum(self.bunch[i].freqs / meansdata[i] for i in procindices)
			# # print 'norm_sum_freqs: ', norm_sum_freqs
			# norm_sum_diversity = get_diversity(norm_sum_freqs)
			# # print 'norm_sum_diversity: ', norm_sum_diversity
			bins = np.bincount(maxes.values(), minlength = 4)
			counts = sorted(bins)
			major = getmax(bins)
			sum_cover = sum(sum_freqs)
			self.diversity[sum_cover] = self.diversity.setdefault(sum_cover, 0) + sum_diversity
			self.covers[sum_cover] = self.covers.setdefault(sum_cover, 0) + 1
			# self.norm_diversity[sum_cover] = self.norm_diversity.setdefault(sum_cover, 0) + norm_sum_diversity
			#norm_sum_cover = sum(norm_sum_cover)
			for i in procindices:
				sorted_freqs = sorted(self.bunch[i].freqs)
				cover = self.bunch[i].summ
				alt_frac = float(sorted_freqs[2]) / cover
				# #major_frac = float(sorted_freqs[3]) / cover
				# if (sorted_freqs >= 2 and alt_frac > 0.25):
				# 	self.hetero_pos[i][cover] = self.hetero_pos[i][cover].setdefault(cover, 0) + 1
				if (sorted_freqs >= 2 and alt_frac > 0.1):
					diversity = get_diversity(self.bunch[i].freqs)
					self.hetero_pos[i][cover] = self.hetero_pos[i].setdefault(cover, 0) + 1
					self.diversities[i][cover] = self.diversities[i].setdefault(cover, 0) + diversity
			if (counts[2] >= occtreshold):
				self.lines = {}
				for i in procindices:
					eq_reads = str(self.bunch[i].freqs[major])
					ind = getmax(self.bunch[i].freqs)
					max_reads = str(self.bunch[i].freqs[ind])
					sum_reads = str(self.bunch[i].summ)
					self.lines[i] = '\t'.join([eq_reads, max_reads, sum_reads]) + '\n'


		#if maxes.values()
		#print 'meansdata: ', meansdata
		#print 'procindices: %s' % procindices
		#if len(procindices) >= 2:
			# all_freqs = [self.bunch[i].freqs for i in self.indices]
			# summary = sum(all_freqs)
			# sum_reads = sum(summary)
			# summary[getmax(summary)] = 0
			# alt_reads = sum(summary)
			# #frac = float(alt_reads) / sum_reads
			#print alt_reads, sum_reads, frac
			# condition = (alt_reads >= ALT_TRESH and sum_reads >= SUM_TRESH)
			# if (condition):
			#for i, j in itertools.combinations(procindices, 2):
				###Another treshold!!!
				#if (check_sample_cond(self.bunch[i].freqs) and check_sample_cond(self.bunch[j].freqs)):
			#	self.cross_cov[i][j] += 1
			#	self.cross_distance[i][j] += (maxes[i] != maxes[j])
		#for i in procindices:
		#	self.cov[i] += 1
		#	self.cross_cov[i][self.refnum] += 1
		#	self.cross_distance[i][self.refnum] += (maxes[i] != self.bunch[i].ref_index)
		for i in self.indices:
			self.bunch[i].readline()

		return(len(procindices) > 0 and counts[2] >= occtreshold)

	def save(self,directory):
		f = open('%s/sumheterogenity%s' % (directory, suffix), 'w')
		for cover in self.diversity.keys():
			f.write('%d\t%f\t%d\n' % (cover, self.diversity[cover], self.covers[cover]))
		f.close()

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
	outfs_hist = [open((os.path.splitext(pname)[0] + '.hethist2'), 'w') for pname in pileupnames]
	for cont in conts:
		flag = pileup_bunch.get_current_cont_indices(cont)
		while flag:
			# print '=================='
			if pileup_bunch.distance(means):
				for num, line in pileup_bunch.lines.items():
					outfs[num].write(line)
			flag = pileup_bunch.get_current_cont_indices(cont)
	pileup_bunch.save(directory)
	for i, fs in enumerate(outfs_hist):
		for cover, div in pileup_bunch.diversities[i].items():
			fs.write('%d\t%f\t%d\n' % (cover, div, pileup_bunch.hetero_pos[i][cover]))
	pileup_bunch.close()

	directory = os.path.dirname(pname)
	outname = (os.path.splitext(pname)[0] + '.refsnphist')

if __name__ == '__main__':
	main()
