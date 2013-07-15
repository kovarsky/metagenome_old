import sys
import random
import os

def getmax(vector):
    """Shows 'letter' with max coverage"""
    maxcov = max(vector)
    ind_of_max = [i for i in xrange(len(vector)) if vector[i] == maxcov]
    ind = random.choice(ind_of_max)
    return(ind)

LETTERS = ("A", "T", "G", "C")
TRESHOLDS = [0.005, 0.01, 0.02, 0.05]

class TabStat():
    """docstring for ClassName"""
    def __init__(self, file_handler, strain_info):
        self.file = open(file_handler, 'r')
        dirname = os.path.dirname(file_handler)
        sample_name = os.path.basename(os.path.splitext(file_handler)[0])
        self.outf = open('%s/%s.fasta' % (dirname, sample_name), 'w')
        self.outf.write('>%s: %s\n' % (strain_info, sample_name))
        self.read_line_or_not = True
        self.symbnum = 0

    def printed(self, index_chrom, index_pos):
        print index_chrom, index_pos

    def readline(self, index_chrom, index_pos):
        try:
            if self.read_line_or_not:
                self.line = self.file.readline().strip().split() 
            chrom, pos = self.line[:2]
            self.read_line_or_not = ((chrom, pos) == (index_chrom, index_pos))
        except ValueError:
            self.read_line_or_not = False

    def getstat(self, index_max_variant, tresh):
        if self.read_line_or_not:
            variants = map(int, self.line[2:6])
            max_variant = getmax(variants)
            self.output = LETTERS[max_variant]
        else:
            self.output = 'N'

    def fastaWrite(self):
        self.outf.write(self.output)
        self.symbnum += 1
        if (self.symbnum % 70 == 0):
            self.outf.write('\n')

    def close(self):
        self.file.close()
        self.outf.close()

def main(index_handler, tab_handlers_list):
    #try:
    index_file = open(index_handler, 'r')
    tab_handlers = [line.strip() for line in open(tab_handlers_list, 'r')]
    dirname = os.path.dirname(index_handler)
    strain_info = open(dirname + '/strain_info','r').readline().strip().split('\t')[3]
    tab_stats = [TabStat(tab_handler, strain_info) for tab_handler in tab_handlers]
    #strain_info = open(dirname + '/strain_info','r').readline().strip().split('\t')[3]
    tresh = 0.01
    for i, index_line in enumerate(index_file):
        stripped_line = index_line.strip().split()
        index_chrom, index_pos = stripped_line[0:2]
        alt_num = float(stripped_line[7])
        sum_num = float(stripped_line[6])
        frac = alt_num/sum_num
        cond = ((alt_num > 5) and (sum_num > 300))
        index_variants = map(int, stripped_line[2:6])
        index_max_variant = getmax(index_variants)
        none_cond = []
        macrocond = (frac >= tresh and cond)
        for tab_stat in tab_stats:
            tab_stat.readline(index_chrom, index_pos)
            if (macrocond):
                tab_stat.getstat(index_max_variant, tresh)
                none_cond.append(tab_stat.output != 'N')

        if (macrocond and any(none_cond)):
            for tab_stat in tab_stats:
                tab_stat.fastaWrite()

        if (i % 10000 == 0):
            print '%dth step' % i

    for tab_stat in tab_stats:
        tab_stat.close()    
    index_file.close()

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tMakeMetaFasta2.py <in: tabs> <tab handler file>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2])
