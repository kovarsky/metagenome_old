import numpy as np
import sys
import os


def main():
	filename = sys.argv[1]
	window_size = int(sys.argv[2])
	directory = os.path.dirname(filename)
	f = open(filename, 'r')
	means_f = open('%s/good_means.list' % directory, 'r')
	means = np.array([float(line.strip()) for line in means_f])
	outf = open('%s/cov_prof_windows_%d.csv' % (directory, window_size), 'w')
	outf_norm = open('%s/normalized_cov_prof_windows_%d.csv' % (directory, window_size), 'w')
	for i, line in enumerate(f):
		vector = np.array([int(el) for el in line.strip().split(',')])
		if not(i % window_size ):
			if i:
				norm_vector = summ_vector / float(window_size)
				# print '------------------'
				# print i
				# print summ_vector
				# print norm_vector
				weight_norm_vector = norm_vector / means
				outf.write(','.join('%.2f' % num for num in norm_vector) + '\n')
				outf_norm.write(','.join('%.3f' % num for num in weight_norm_vector) + '\n')
			summ_vector = vector
		else:
			# print 'vector: ', vector
			summ_vector += vector
	if (i % window_size):
		norm_vector = summ_vector / float(i % window_size + 1)
		weight_norm_vector = norm_vector / means
		outf.write(','.join('%.2f' % num for num in norm_vector) + '\n')
		outf_norm.write(','.join('%.3f' % num for num in weight_norm_vector) + '\n')
	outf.close()
	outf_norm.close()

if __name__ == '__main__':
	main()
