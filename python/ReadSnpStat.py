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
    def __init__(self, file_handler):
        self.file = open(file_handler, 'r')
        self.read_line_or_not = True
        self.total_cover = {tresh: 0 for tresh in TRESHOLDS}
        self.snp_num = {tresh: 0 for tresh in TRESHOLDS}
        self.sample_name = os.path.basename(os.path.splitext(file_handler)[0])

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
            self.total_cover[tresh] += 1
            variants = map(int, self.line[2:6])
            max_variant = getmax(variants)
            self.snp_num[tresh] += (max_variant != index_max_variant)
    def write(self, fileObj):
        for tresh in TRESHOLDS:
            fileObj.write('%s\t%s\t%s\n' % (self.sample_name, self.total_cover[tresh], self.snp_num[tresh]))

    def close(self):
        self.file.close()

def main(index_handler, tab_handlers_list):
    #try:
    index_file = open(index_handler, 'r')
    tab_handlers = [line.strip() for line in open(tab_handlers_list, 'r')]
    tab_stats = [TabStat(tab_handler) for tab_handler in tab_handlers]
    dirname = os.path.dirname(index_handler)
    #strain_info = open(dirname + '/strain_info','r').readline().strip().split('\t')[3]

    for i, index_line in enumerate(index_file):
        stripped_line = index_line.strip().split()
        index_chrom, index_pos = stripped_line[0:2]
        alt_num = float(stripped_line[7])
        sum_num = float(stripped_line[6])
        try:
            frac = alt_num/sum_num
        except ZeroDivisionError:
            frac = 0
        cond = ((alt_num > 5) and (sum_num > 300))
        index_variants = map(int, stripped_line[2:6])
        index_max_variant = getmax(index_variants)
        for tab_stat in tab_stats:
            tab_stat.readline(index_chrom, index_pos)
            for tresh in TRESHOLDS:
                if (frac >= tresh and cond):
                    tab_stat.getstat(index_max_variant, tresh)
        if (i % 10000 == 0):
            print '%dth step' % i

    outf = open('%s/new_snp_stat2.tsv' % dirname, 'w')
    for tab_stat in tab_stats:
        tab_stat.write(outf)
        tab_stat.close()    
    outf.close()
    index_file.close()

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tReadSnpStat.py <in: sorted index file> <tab handler file>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2])
