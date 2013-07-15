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


def main(index_handler, tab_handler):
    #try:
    index_file = open(index_handler, 'r')
    tab_file = open(tab_handler, 'r')
    dirname = os.path.dirname(index_handler)

    strain_info = open(dirname + '/strain_info','r').readline().strip().split('\t')[3]
    sample_name = os.path.basename(os.path.splitext(tab_handler)[0])
    outf = open('%s/%s.fasta' % (dirname, sample_name), 'w')
    outf.write('>%s: %s\n' % (strain_info, sample_name))
    read_line_or_not = True
    j = 0
    for i, index_line in enumerate(index_file):
        stripped_line = index_line.strip().split()
        ind_chrom, ind_pos = stripped_line[0:2]
        alt_read_num = int(stripped_line[7])
        try:
            if read_line_or_not:
                line = tab_file.readline().strip().split()
                chrom, pos = line[:2]
                variants = map(int, line[2:])
            if (chrom, pos) == (ind_chrom, ind_pos):
                read_line_or_not = True
                value = LETTERS[getmax(variants)]
            else:
                read_line_or_not = False;
                value = 'N'
        except ValueError:
            read_line_or_not = True;
            value = 'N'
        if alt_read_num > 9:
            j += 1 
            outf.write(value)
            if ((j + 1) % 70 == 0):
                outf.write('\n')
        else:
            pass
        if (i % 10000 == 0):
            print '%dth step' % i
    outf.close()
    tab_file.close()
    index_file.close()

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tMakeMetaFasta.py <in: sorted index file> <*.tab file>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2])
