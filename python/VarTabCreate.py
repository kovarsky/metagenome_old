import sys
import random
import re
import sets
import json
from sets import Set


def heterogen(vector):
    """Heteronity meadure function. It should be noted
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


def pileup_dict_create(f_pileup):
    """From pileup file creates dictionary
    in the form {chrom:{pos:[letter,deep]}}
    deep means coverage in this position"""
    pileup_dict = {}
    for i, line in enumerate(f_pileup):
        try:
            words = line.strip().split()
            #print(words)
            chrom = words[0]
            pos = int(words[1])
            letter = words[2]
            deep = words[3]
            if chrom not in pileup_dict.keys():
                pileup_dict[chrom] = {}
            pileup_dict[chrom][pos] = [letter, deep]
        except IndexError:
            print('IndexError:')
            print line
            print i
    return(pileup_dict)

def joining(varpos_set, snptab_set, pileup_dict, f_vartab):
    """First, this function computes difference between varpos
    file and *.snptab file. Second positions that are contained
    in variable position list but not in *.snptab files are added
    into *.vartab file"""
    letterdict = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

    diff = varpos_set - snptab_set
    diff_list = []
    for el in diff:
        el = el.strip().split()
        chrom = el[0]
        pos = int(el[1])
        diff_list.append([chrom, pos])

    print 'Joining'

    for el in diff_list:
        try:
            letters = ['0', '0', '0', '0']
            buff = pileup_dict[el[0]][el[1]]
            letters[letterdict[buff[0]]] = buff[1]
            f_vartab.write('%s\t%d\t%s\n' % (el[0], el[1], '\t'.join(letters)))
        except KeyError:
            pass

def main(handler, varposfilename, chromlistfile, tresholds):
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

    f_snptab_list = open(handler, 'r')
    f_varpos = open(varposfilename, 'r')
    varpos = {line.strip() for line in f_varpos}
    varpos_set = sets.Set(varpos)
    for num_of_sample, fname in enumerate(f_snptab_list):
        print ('working with ' + fname)
        fname = fname.strip()
        f_snptab = open(fname, 'r')
        snptab = {'\t'.join(line.split()[0:2]) for line in f_snptab}
        snptab_set = sets.Set(snptab)

        filename_pileup = re.sub('snptab$', 'pileup', fname)
        f_pileup = open(filename_pileup, 'r')
        filename_vartab = re.sub('snptab$', 'vartab', fname)
        f_vartab = open(filename_vartab, 'w')

        print 'Creating pileup_dict'
        pileup_dict = pileup_dict_create(f_pileup)

        joining(varpos_set, snptab_set, pileup_dict, f_vartab)

        f_vartab.close()
        f_pileup.close()
        f_snptab.close()

        f_snptab = open(fname, 'r')
        f_vartab = open(filename_vartab, 'r')
        vartab = [intconv(line.strip().split()) for line in f_vartab]
        snptab = [intconv(line.strip().split()) for line in f_snptab]
        tab = vartab + snptab
        tab = sorted(tab)
        f = open(fname[:-6] + 'tab', 'w')
        for el in tab:
            statfunc(el, num_of_sample)
            f.write('\t'.join(strconv(el)) + '\n')
        f.close()
        f_vartab.close()
        f_snptab.close()

    path = re.match(r'.+(?=/\w)', handler).group()

    outf = open('%s/json.statistics.txt' % path, 'w')
    outf.write(json.dumps(total_cover_data) + '\n')
    outf.write(json.dumps(heterogen_data) + '\n')

    outf.close()


if __name__ == '__main__':
    if (len(sys.argv) < 3):
        print 'usage:\tVarTabCreate.py <in: snptab list> <in: varposfilename> <in: contigs list> <in: tresholds>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
