import sys
import random
import re
#from time import time
#from sets import Set
#import json
import numpy as np
import os
import itertools



POSITIONS = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
TRESHOLD = 3

def getmax(vector):
    """Shows 'letter' with max coverage"""
    maxcov = max(vector)
    ind_of_max = [i for i in xrange(4) if vector[i] == maxcov]
    ind = random.choice(ind_of_max)
    return(ind)

        
class Reference():
    """docstring for Reference"""
    def __init__(self, filename):
        #super(Reference, self).__init__()
        #self.arg = arg
        self.data = open(filename, 'r')
        self.next()

    def next(self):
        letter = self.data.readline().strip()
        try:
            self.maxnum = POSITIONS[letter]
        except KeyError:
            pass


class Tab():
    """docstring for Tab"""
    def __init__(self, tabname):
        self.data = open(tabname, 'r')
        self.readline()
        self.not_eof = True
        # self.readflag = True
        # self.contigflag = False

    def readline(self):
        line = self.data.readline().strip()
        if line != '':
            splitted_line = line.split('\t')
            self.seq = splitted_line[0]
            self.pos = int(splitted_line[1])
            self.freqs = [int(splitted_line[i]) for i in xrange(2,6)]
        else:
            self.seq = None
            self.pos = None
            self.not_eof = False

    def close(self):
        self.data.close()        

class BunchOfTabs(object):
    """docstring for BunchOfObjects"""
    def __init__(self, tabnames, refname):
        self.bunch = [Tab(tabname) for tabname in tabnames]
        self.ref = Reference(refname)
        length = len(tabnames)
        self.refnum = length
        self.cross_cov = np.zeros((length + 1, length + 1))
        self.distance = np.zeros((length + 1, length + 1))

    def get_indices(self, current_cont, current_pos):
        self.indices = [i for i, obj in enumerate(self.bunch) if 
                        (obj.seq == current_cont and obj.pos == current_pos and obj.not_eof)]
        
    def comparison(self):
        maxes = {i: getmax(self.bunch[i].freqs) for i in self.indices}
        sums = {i: float(sum(self.bunch[i].freqs)) for i in self.indices}
        def condfunc(ind):
            max_cond = (maxes[ind] >= 2)
            sums_cond = (sums[ind] >= 2)
            return (max_cond and sums_cond)
        self.compindices = [i for i in self.indices if condfunc(i)]
        for i in self.compindices:
            self.cross_cov[i][self.refnum] += 1
            self.distance[i][self.refnum] += (maxes[i] != self.ref.maxnum)
        if len(self.compindices) >= 2:
            for i, j in itertools.combinations(self.compindices, 2):
                self.cross_cov[i][j] += 1
                self.distance[i][j] += (maxes[i] != maxes[j])
    def next(self):
        for i in self.indices:
            self.bunch[i].readline()
        self.ref.next()
    
    def close(self):
        for obj in self.bunch:
            obj.close()




def main(index_filename, tab_list_filename, prefix):
    #plistname = sys.argv[1]
    #contlistname = sys.argv[2]
    index_data = open(index_filename, 'r')
    tab_list = open(tab_list_filename, 'r')
    tabnames = [line.strip() for line in tab_list]
    short_names = [os.path.basename(os.path.splitext(el)[0]) for el in tabnames]
    directory = os.path.dirname(tabnames[0])
    refname = '%s/reference.list' % directory
    tabs = BunchOfTabs(tabnames, refname)
    for entry in index_data:
        entry_splitted = entry.strip().split()
        cont = entry_splitted[0]
        pos = int(entry_splitted[1])
        tabs.get_indices(cont, pos)
        #TODO:
        condition = True
        #
        if condition:
            # print entry
            # print tabs.ref.maxnum
            # print '----'
            tabs.comparison()
        tabs.next()

    head = '\t'.join(short_names)
    np.savetxt('%s/dist_cross_cov_%s.out' % (directory, prefix), tabs.cross_cov,delimiter='\t',fmt='%.1f',header= head)
    np.savetxt('%s/dist_cross%s.out' % (directory, prefix), tabs.distance,delimiter='\t',fmt='%.1f', header = head)





if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tBigTab3.py <in: sorted index file> <in: sorted *.tab filename list> <in: prefix>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
