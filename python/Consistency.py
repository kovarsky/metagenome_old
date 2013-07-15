import sys
import os
import json
#import numpy as np
#import random
import itertools


path = sys.argv[1]
#times = int(sys.argv[2])
#frac = float(sys.argv[3])

def lister(elements):
    if isinstance(elements, basestring):
        return([elements])
    else:
        return(elements)

class Split():
    """docstring for Split"""
    def __init__(self, lists, namesdict):
        self.gr_A = set([namesdict[name] for name in lister(lists[0])])
        self.gr_B = set([namesdict[name] for name in lister(lists[1])])
    def consistance(self, gr_a, gr_b):
        gr_A = self.gr_A
        gr_B = self.gr_B
        aA = len(gr_a & gr_A)
        #print('aA ', aA)
        aB = len(gr_a & gr_B)
        #print('aB ', aB)
        bA = len(gr_b & gr_A)
        #print('bA ', bA)
        bB = len(gr_b & gr_B)
        #print('bB ', bB)
        inter1 = aA + bB
        inter2 = aB + bA
        self.cond = ((len(gr_a & gr_A) + len(gr_a & gr_B) > 0)
                    & (len(gr_b & gr_A) + len(gr_b & gr_B) > 0))
        if self.cond:
            self.cons = 0.5 * abs((aA - aB) / float(aA + aB) + (bB - bA) / float(bB + bA))
            # if inter1 >= inter2:
            #     self.cons = 0.5 * abs((aA - aB) / float(aA + aB) + (bB - bA) / float(bB + bA))
            # else:
            #     self.cons = 0.5 * abs((aB - aA) / float(aA + aB) + (bA - bB) / float(bB + bA))
        else:
            self.cons = 0

class Splits():
    """docstring for Splits"""
    def __init__(self, tree, namesdict):
        self.data = [Split(lists, namesdict) for lists in tree]
        splitnum = len(tree)
        self.summary = [0 for i in range(splitnum)]
        self.consistances = [0.0 for i in range(splitnum)]

    def get_consistance(self, gr_a, gr_b):
        for i, split in enumerate(self.data):
            # raw_input("--------")
            split.consistance(gr_a, gr_b)
            #print ("Polymorphism splits: ")
            #print ("a ", gr_a, len(gr_a))
            #print ("b ", gr_b, len(gr_b))
            #print ("Tree splits:")
            #print ("A", split.gr_A, len(split.gr_A))
            #print ("B", split.gr_B, len(split.gr_B))
            #print ("consistence: ", split.cons)
            #print ("condition: ", split.cond)
            if split.cond:
                self.summary[i] += 1
                self.consistances[i] += split.cons
            else:
                pass

    def close(self, path):
        f = open('%s/consistancies.list' % path, 'w')
        for i, split in enumerate(self.data):
            try:
                f.write('%.2f\n' % (self.consistances[i] / self.summary[i]))
            except ZeroDivisionError:
                f.write('0.00\n')
        f.close()

def main():
    numfile = open('%s/nums.list' % path, 'r')
    nums = [(int(line.strip()) - 1) for line in numfile]
    datafile = open('%s/snp_chars_alt_4' % path, 'r')
    pileupsfile = open('%s/pileups.list' % path, 'r')
    namesdict = {os.path.basename(os.path.splitext(line.strip())[0]): i
                  for i, line in enumerate(pileupsfile) if (i in nums)}
    #print namesdict
    #print nums
    length = len(nums)
    treshold = length - 2
    treef = open('%s/newtree.json' % path, 'r')
    tree = json.loads(treef.readline().strip())
    splits = Splits(tree, namesdict)
    for line in datafile:
        letters = [line[num] for num in nums]
        #print ('=======')
        #print(letters)
        bases = set([letter for letter in letters if letter != '-'])
        #print bases
        if (letters.count('-') < treshold) and (len(bases) == 2):
            gr_a = set([num for num in nums if line[num] == list(bases)[0]])
            gr_b = set([num for num in nums if line[num] == list(bases)[1]])
            splits.get_consistance(gr_a, gr_b)
    splits.close(path)
    datafile.close()
    pileupsfile.close()
    treef.close()
    numfile.close()



if __name__ == '__main__':
    main()
