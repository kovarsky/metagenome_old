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
#ALT_TRESH = 'new'
#SUM_TRESH = 2
#FRAC_TRESH = 0
#TRESHOLDS = '_'.join([str(ALT_TRESH), str(SUM_TRESH), str(FRAC_TRESH)])
LETTERS = ['A', 'T', 'G', 'C']


class BunchOfObjects(CrossHeterogenity.BunchOfObjects):
    """docstring for BunchOfObjects"""
    def __init__(self, pileupnames):
        self.bunch = [CrossHeterogenity.PileupObj(pileupname) for pileupname in pileupnames]
        length = len(pileupnames)
        # self.cross_cov = np.zeros((length + 1, length + 1))
        # self.cross_distance = np.zeros((length + 1, length + 1))
        # self.cov = np.zeros((length))
        self.refnum = length
        self.position_count = {}
        self.read_count = {}
        #self.position_count = np.zeros

    def get_current_cont_indices(self, current_cont):
        self.indices = [i for i, obj in enumerate(self.bunch) if (obj.seq == current_cont and obj.not_eof)]
        #print self.bunch[17].seq, current_cont, self.bunch[17].coord, self.bunch[17].not_eof, len(self.indices)
        if len(self.indices) < 1:
            self.indices = range(len(self.bunch))
            return False
        else:
            return True

    def getline(self):
        positions = [self.bunch[i].coord for i in self.indices]
        # print 'positions: %s' % positions
        min_pos = min(positions)
        # print 'min position: %s' % min_pos
        self.min_pos = min_pos
        self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
        # print 'indices: %s' % self.indices
        # print '-------------'
        # for j in self.indices:
        #   print 'line%d: %s' % (j, self.bunch[j].line)
        # print '-------------'
        maxes = {i: getmax(self.bunch[i].freqs) for i in self.indices}
        # print 'maxes: %s' % maxes
        procindices = [i for i in self.indices if self.bunch[i].freqs[maxes[i]] >= treshold ]
        self.position_count[len(procindices)] = self.position_count.setdefault(len(procindices), 0) + 1
        new_maxes = [maxes[i] for i in procindices]
        #counts = sorted(np.bincount(new_maxes, minlength = 4))
        #print counts
        #occur = counts[2]
        # print 'procindices: %s' % procindices
        #print ()
        self.procindices_num = len(procindices)
        #self.procsumm_num = sum(self.bunch[i].summ for i in procindices)
        self.summs = [self.bunch[i].summ for i in procindices]
        self.procindices = procindices
        #if (occur >= occtreshold):
        if len(procindices) >= 1:
                cond_consistancy = True
                #self.freqs = '\t'.join(str(el) for el in sum(self.bunch[i].freqs for i in procindices))
                self.line = ','.join([ str(self.bunch[i].summ) if i in procindices else '0' for i in range(self.refnum)]) + '\n'
        else:
            cond_consistancy = False

            #   # all_freqs = [self.bunch[i].freqs for i in self.indices]
            #   # summary = sum(all_freqs)
            #   # sum_reads = sum(summary)
            #   # summary[getmax(summary)] = 0
            #   # alt_reads = sum(summary)
            #   # #frac = float(alt_reads) / sum_reads
            #   #print alt_reads, sum_reads, frac
            #   # condition = (alt_reads >= ALT_TRESH and sum_reads >= SUM_TRESH)
            #   # if (condition):
            #   for i, j in itertools.combinations(procindices, 2):
            #       ###Another treshold!!!
            #       #if (check_sample_cond(self.bunch[i].freqs) and check_sample_cond(self.bunch[j].freqs)):
            #       self.cross_cov[i][j] += 1
            #       self.cross_distance[i][j] += (maxes[i] != maxes[j])
            # for i in procindices:
            #   self.cov[i] += 1
            #   self.cross_cov[i][self.refnum] += 1
            #   self.cross_distance[i][self.refnum] += (maxes[i] != self.bunch[i].ref_index)
        for i in self.indices:
            #print 'freqs%d: %s' % (i, self.bunch[i].freqs)
            self.bunch[i].readline()
        return(cond_consistancy)

    def get_normalized(self, cover_data):
        new_summs = sum(self.summs[i] / cover_data[j] if cover_data[j] != 0 else 0 for i, j in enumerate(self.procindices))
        return (new_summs)


    # def save(self,directory):
    #   np.savetxt('%s/cross_diffs_%s.out' % (directory, suffix), self.cross_distance, delimiter='\t',fmt='%.1f')
    #   np.savetxt('%s/cross_cov_%s.out' % (directory, suffix), self.cross_cov, delimiter='\t',fmt='%.1f')
    #   np.savetxt('%s/cov_%s.out' % (directory, suffix), self.cov, delimiter='\t',fmt='%.1f')


def main():
    plistname = sys.argv[1]
    contlistname = sys.argv[2]
    plist = open(plistname, 'r')
    contlist = open(contlistname, 'r')
    pileupnames = [line.strip() for line in plist]
    directory = os.path.dirname(pileupnames[0])
    conts = [line.strip() for line in contlist]
    pileup_bunch = BunchOfObjects(pileupnames)
    total_cov_f = open('%s/means.list' % directory, 'r')
    #total_cov_gen = (line.strip().split() for line in total_cov_f)
    total_cov_data = [float(line.strip()) for line in total_cov_f]
    #print total_cov_data
    logf = open('%s/logfile_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    outf = open('%s/cov_corr_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    # outf_poscount = open('%s/poscount_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    # outf_profile = open('%s/mutual_profile_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    # outf_normprofile = open('%s/normalized_mutual_profile_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    # outf_readprofile = open('%s/mutual_readprofile_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    # outf_freqs = open('%s/mutual_freqs_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
    linenum = 0
    for cont in conts:
        #print '___________%s__' % cont
        flag = pileup_bunch.get_current_cont_indices(cont)
        while flag:
            # print '=================='
            # pileup_bunch.distance()
            if pileup_bunch.getline():
                #print pileup_bunch.min_pos
                # outf_freqs.write('%d\t%s\n' % (pileup_bunch.min_pos, pileup_bunch.freqs))
                outf.write('%s' % (pileup_bunch.line))
                #print pileup_bunch.min_pos, pileup_bunch.line
            flag = pileup_bunch.get_current_cont_indices(cont)
            if not(linenum % 100000):
                logf.write( '%d\n' % linenum)
            linenum += 1
            # outf_profile.write('%d\t%d\n' % (pileup_bunch.min_pos, pileup_bunch.procindices_num))
            # outf_readprofile.write('%d\t%d\n' % (pileup_bunch.min_pos, pileup_bunch.procsumm_num))
            # outf_normprofile.write('%d\t%f\n' % (pileup_bunch.min_pos, pileup_bunch.get_normalized(total_cov_data)))
    # for key in sorted(pileup_bunch.position_count.keys()):
        # outf_poscount.write('%d\t%d\n' % (key, pileup_bunch.position_count[key]))
    logf.close()
    # outf_poscount.close()
    # pileup_bunch.save(directory)
    pileup_bunch.close()

if __name__ == '__main__':
    main()
