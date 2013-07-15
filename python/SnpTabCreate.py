#!/bin/env python
import sys
import random
import re


SUBTABLE = (('A\tR',  'A\tG'),
            ('G\tR',  'G\tA'),
            ('G\tS',  'G\tC'),
            ('C\tS',  'C\tG'),
            ('C\tY',  'C\tT'),
            ('T\tY',  'T\tC'),
            ('A\tW',  'A\tT'),
            ('T\tW',  'T\tA'),
            ('G\tK',  'G\tT'),
            ('T\tK',  'T\tG'),
            ('A\tM',  'A\tC'),
            ('C\tM',  'C\tA'))


def translateiupac(line):
    for sub in SUBTABLE:
        line = re.sub(sub[0], sub[1], line)
    return(line)


def main(inputsnps):

    letterdict = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    outname = re.sub('\.\w+$', '.snptab', inputsnps)
    outf = open(outname, 'w')
    with open(inputsnps) as f:
        header = f.readline()
        letters = ['0', '0', '0', '0']
        for i, line in enumerate(f):
            line = translateiupac(line)
            words = line.strip().split()
            pos = words[1]
            chrom = words[0]
            if i == 0:
                oldpos = pos
                oldchrom = chrom
            if oldpos != pos:
                outf.write('%s\t%s\t%s\n' %
                    (oldchrom,  oldpos,  '\t'.join(letters)))
                letters = ['0', '0', '0', '0']
            ref = words[2]
            alt = words[3]
            refcov = words[4]
            altcov = words[5]
            if ref != 'N':
                letters[letterdict[ref]] = refcov
            letters[letterdict[alt]] = altcov
            oldpos = pos
            oldchrom = chrom
        try:
            outf.write('%s\t%s\t%s\n' %
                    (oldchrom,  oldpos,  '\t'.join(letters)))
        except UnboundLocalError:
            print outname
    outf.close()
        #patternlist = [line.strip('\n\r') for line in f]


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tFastaProcessing.py <in: *.snps file>'
        sys.exit(0)
    else:
        main(sys.argv[1])
