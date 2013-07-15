import sys
import os
import json
#import numpy as np
#import random
import itertools


path = sys.argv[1]
#times = int(sys.argv[2])
#frac = float(sys.argv[3])

class Split():
    """docstring for Split"""
    def __init__(self, lists, namesdict):
        self.gr_A = set([namesdict[name] for name in lists[0]])
        self.gr_B = set([namesdict[name] for name in lists[1]])
    def consistance(self, gr_a, gr_b):
        gr_A = self.gr_A
        gr_B = self.gr_B
        aA = len(gr_a & gr_A)
        aB = len(gr_a & gr_B)
        bA = len(gr_b & gr_A)
        bB = len(gr_b & gr_B)
        inter1 = aA + bB
        inter2 = aB + bA
        if inter1 >= inter2:
            self.cons = 0.5 * abs((aA - aB) / len(gr_a) + (bB - bA) / len(gr_b))
        else:
            self.cons = 0.5 * abs((aB - aA) / len(gr_a) + (bA - bB) / len(gr_b))
        self.cond = ((len(gr_a & gr_A) + len(gr_a & gr_B) > 0)
                    & (len(gr_b & gr_A) + len(gr_b & gr_B) > 0))

class Splits():
    """docstring for Splits"""
    def __init__(self, tree):
        self.data = [Split(lists) for lists in tree]
        splitnum = len(tree)
        self.summary = [0 for i in splitnum]
        self.consistances = [0.0 for i in splitnum]

    def get_consistance(self, gr_a, gr_b):
        for i, split in enumerate(self.data):
            if split.cond:
                self.summary[i] += 1
                self.consistances[i] += split.cons
            else:
                pass

    def close(self, path):
        f = open('%s/consistances.list', 'w')
        for i, split in enumerate(self.data):
            f.write('%.2f\n' % (self.consistances[i] / self.summary[i]))
        f.close()

def main():
    numfile = open('%s/nums.list' % path, 'r')
    nums = [(int(line.strip()) - 1) for line in numfile]
    datafile = open('%s/snp_chars_alt_4' % path, 'r')
    pileupsfile = open('%s/pileups.list' % path, 'r')
    namesdict = {os.path.basename(os.path.splitext(line.strip())[0]): i
                  for i, line in enumerate(pileupsfile) if (i in nums)}
    length = len(nums)
    treef = open('%s/tree.json' % path, 'r')
    tree = json.load(treef.readline().strip())
    splits = Splits(tree, namesdict)
    for line in datafile:
        letters = [line[num] for num in nums]
        bases = set([letter for letter in letters if letters != '-'])
        if (newline.count('-') < treshold) and (len(bases) == 2):
            gr_a = set([num for num in nums if line[num]=bases[0]])
            gr_b = set([num for num in nums if line[num]=bases[1]])
            splits.get_consistance(gr_a, gr_b)
    splits.close()
    datafile.close()
    pileupsfile.close()
    treef.close()
    numfile.close()



if __name__ == '__main__':
    main()