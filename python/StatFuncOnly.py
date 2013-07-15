import sys
import random
import re
import sets
import json
from sets import Set



def heterogen(vector):
    """Heteronity measure function. It should be noted
    that I substitute sum(vector) from sum_squared.
    Some mistakes can arise due to this purpose"""
    vector = [el for el in vector if el > 0]
    sum_squared = (sum(vector)) ** 2
    squares = map(lambda x: x * (x - 1), vector)
    delta = sum_squared - sum(vector)
    if delta > 0:
        prob = (delta - sum(squares)) / float(delta)
    else:
        prob = 0
    return prob

def main(handler, chromlistfile, tresholds):
    tresholds = [int(tresh) for tresh in tresholds]

    chroms = Set([line.strip() for line in open(chromlistfile, 'r')])
    fnamelist = [line for line in open(handler, 'r')]
    filenum = len(fnamelist)

    total_cover_data = [{chrom: [0] * filenum for chrom in chroms} for i in tresholds]
    heterogen_data = {chrom: [0] * filenum for chrom in chroms}

    def intconv(arg):
        """Aux func for conversion second element of line
        into int"""
        return([arg[0], int(arg[1])] + arg[2:])

    def strconv(arg):
        """Aux func for conversion second element of line
        into str"""
        return([arg[0], str(arg[1])] + arg[2:])

    def statfunc(line, num_of_sample):
        vector = line[2:6]
        vector = map(lambda x: int(x), vector)
        for j, treshold in enumerate(tresholds):
            if sum(vector) >= treshold:
                total_cover_data[j][line[0]][num_of_sample] += 1
        heterogen_data[line[0]][num_of_sample] += heterogen(vector)

    f_tab_list = open(handler, 'r')
    for num_of_sample, fname in enumerate(f_tab_list):
        print ('working with ' + fname)
        fname = fname.strip()
        f_tab = open(fname, 'r')
        for line in f_tab:
            el = intconv(line.strip().split())
            statfunc(el, num_of_sample)
        f_tab.close()

    path = re.match(r'.+(?=/\w)', handler).group()

    outf = open('%s/json.statistics.txt' % path, 'w')
    outf.write(json.dumps(total_cover_data) + '\n')
    outf.write(json.dumps(heterogen_data) + '\n')

    outf.close()


if __name__ == '__main__':
    if (len(sys.argv) < 3):
        print 'usage:\tVarTabCreate.py <in: snptab list> <in: contigs list> <in: tresholds>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3:])
