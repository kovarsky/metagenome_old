import sys
import os
import numpy as np
import random
import itertools


path = sys.argv[1]
times = int(sys.argv[2])
frac = float(sys.argv[3])



def main():
    numfile = open('%s/nums.list' % path, 'r')
    nums = [(int(line.strip()) - 1) for line in numfile]
    # print nums
    # pileupsfile = open('%s/pileups.list' % path, 'r')
    # namesdict = {i: os.path.basename(os.path.splitext(line.strip())[0])
    #              for i, line in enumerate(pileupsfile) if (i in nums)}
    # pileupsfile = open('%s/pileups.list' % path, 'r')
    # names = [os.path.basename(os.path.splitext(line.strip())[0])
    #              for i, line in enumerate(pileupsfile) if (i in nums)]
    length = len(nums)
    for time in range(times):
        datafile = open('%s/snp_chars_filtered_4' % path, 'r')
        cross_cov = np.zeros((length, length))
        cross_distance = np.zeros((length, length))
        for line in datafile:
            if random.random() <= frac:
                goodnums = [i for i in xrange(length) if line[i] != '-']
                for i, j in itertools.combinations(goodnums, 2):
                    cross_cov[i][j] += 1
                    cross_distance[i][j] += (line[i] != line[j])
        datafile.close()
        np.savetxt('%s/sampled_cross_diffs_%.1f_%d.out' % (path, frac, time), cross_distance, delimiter='\t',fmt='%.1f')
        np.savetxt('%s/sampled_cross_cov_%.1f_%d.out' % (path, frac, time), cross_cov, delimiter='\t',fmt='%.1f')


    numfile.close()

if __name__ == '__main__':
    main()