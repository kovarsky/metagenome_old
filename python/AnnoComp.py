import sys
import os
import numpy as np
import re
import itertools
import random


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


#suffix = sys.argv[3]
treshold = 4
LETTERS = ('A', 'T', 'G', 'C')


class PileupObj():
    """docstring for PileupObj"""
    def __init__(self, pileupname):
        self.data = open(pileupname, 'r')
        self.readline()
        self.not_eof = True
        # self.readflag = True
        # self.contigflag = False

    def readline(self):
        line = self.data.readline()
        if line != '':
            self.seq, self.coord, ref, self.summ, reads, quals = line.split('\t')
            self.coord = int(self.coord)
            reads = re.sub('[\^].|\$', '', reads)
            ref_index = LETTERS.index(ref)
            freqs = np.array([reads.upper().count(letter) for letter in LETTERS])
            freqs[ref_index] = int(self.summ) - sum(freqs)
            self.summ = float(self.summ)
            self.maxx = getmax(freqs)
            self.snp = (self.maxx != ref_index)
            self.support = freqs[self.maxx]
        else:
            self.seq = None
            self.not_eof = False

    def close(self):
        self.data.close()


class BunchOfObjects():
    """docstring for BunchOfObjects"""
    def __init__(self, pileupnames):
        self.bunch = [PileupObj(pileupname) for pileupname in pileupnames]
        self.len = len(pileupnames)

    def get_current_cont_indices(self, current_cont):
        self.indices = [i for i, obj in enumerate(self.bunch) if (obj.seq == current_cont and obj.not_eof)]
        if len(self.indices) < 1:
            self.indices = range(len(self.bunch))
            return False
        else:
            return True

    def getline(self):
        positions = [self.bunch[i].coord for i in self.indices]
        min_pos = min(positions)
        self.indices = [i for i in self.indices if self.bunch[i].coord == min_pos]
        covers = [self.bunch[i].summ if i in self.indices else 0 for i in range(self.len)]
        procindices = [i for i in self.indices if (self.bunch[i].support >= treshold)]
        total = len(procindices)
        variants = [self.bunch[i].maxx for i in procindices if self.bunch[i].snp]
        uniq_vars = set(variants)
        snps = [(str(min_pos), LETTERS[var], str(variants.count(var)), str(total)) for var in uniq_vars]
        self.pos = min_pos
        self.info = [covers, snps]
        for i in self.indices:
            self.bunch[i].readline()

    def close(self):
        for pileup in self.bunch:
            pileup.close()

            #variants = [[pos, letter, occurence, total]]


class Genes(object):
    """docstring for Genes"""
    def __init__(self, annofilename, coverfilename, snpfilename, samplenums):
        print coverfilename
        print snpfilename
        self.coverfile = open(coverfilename, 'w')
        self.snpfile = open(snpfilename, 'w')
        self.annofile = open(annofilename, 'r')
        self.data = [Gene(line, samplenums) for line in self.annofile]

    def get_genes(self, pos, info):
        iterator = itertools.takewhile(lambda x: x.start <= pos, self.data)
        cur_gene_inds = (i for i, gene in enumerate(iterator))
        for i in cur_gene_inds:
            if self.data[i].end < pos:
                self.data[i].write(self.coverfile, self.snpfile)
                del self.data[i]
            else:
                self.data[i].add(info)

    def close_files(self):
        self.coverfile.close()
        self.snpfile.close()
        self.annofile.close()


class Gene(object):
    """docstring for Gene"""

    def __init__(self, line, samplenums):
        elems = line.strip().split()
        self.id = elems[0]
        self.start = int(elems[1])
        self.end = int(elems[2])
        self.snps = []
        self.covers = np.zeros(samplenums)
        self.counts = np.zeros(samplenums)

    def add(self, info):
        self.covers += np.array(info[0])
        self.counts += np.array([cover != 0 for cover in info[0]])
        self.snps += info[1]
        #self.snps.append(info[1])

    def write(self, coverfile, snpfile):
        # print self.covers
        # print self.counts
        # print self.snps
        frac = self.covers / self.counts
        frac[np.where(np.isnan(frac))] = 0
        res = '\t'.join(['%.2f' % el for el in frac])
        coverfile.write('%s\t%s\n' % (self.id, res))
        for snp in self.snps:
            snpfile.write('%s\t%s\n' % (self.id, '\t'.join(snp)))


def main():
    plistname = sys.argv[1]
    contlistname = sys.argv[2]
    annofilename = sys.argv[3]
    prefix = sys.argv[4]

    plist = open(plistname, 'r')
    contlist = open(contlistname, 'r')
    pileupnames = [line.strip() for line in plist]
    samplenums = len(pileupnames)
    pileup_bunch = BunchOfObjects(pileupnames)

    directory = os.path.dirname(pileupnames[0])
    
    coverfilename = '%s/%s_gene_covers.tsv' % (directory, prefix)
    snpfilename = '%s/%s_gene_snps.tsv' % (directory, prefix)
    genes = Genes(annofilename, coverfilename, snpfilename, samplenums)

    conts = [line.strip() for line in contlist]
    for cont in conts:
        flag = pileup_bunch.get_current_cont_indices(cont)
        while flag:
            pileup_bunch.getline()
            genes.get_genes(pileup_bunch.pos, pileup_bunch.info)
            flag = pileup_bunch.get_current_cont_indices(cont)

    genes.close_files()
    pileup_bunch.close()

if __name__ == '__main__':
    main()