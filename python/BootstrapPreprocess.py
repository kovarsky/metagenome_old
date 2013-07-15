import sys
import os
#import numpy as np
#import random
import itertools


path = sys.argv[1]
#times = int(sys.argv[2])
#frac = float(sys.argv[3])



def main():
    numfile = open('%s/nums.list' % path, 'r')
    nums = [(int(line.strip()) - 1) for line in numfile]
    datafile = open('%s/snp_chars_alt_4' % path, 'r')
    outf = open('%s/snp_chars_filtered_4' % path, 'w')
    # pileupsfile = open('%s/pileups.list' % path, 'r')
    # names = [os.path.basename(os.path.splitext(line.strip())[0])
    #              for i, line in enumerate(pileupsfile) if (i in nums)]
    length = len(nums)
    treshold = length - 2
    for line in datafile:
        newline = ''.join([line[num] for num in nums])
        #goodnums = [i for i, num in enumerate(nums) if line[num] != '-']
        if newline.count('-') < treshold:
            outf.write(newline + '\n')
    datafile.close()
    

if __name__ == '__main__':
    main()