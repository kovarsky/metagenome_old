import sys
import random
import re
from time import time
from sets import Set
import json


def getmax(vector):
    """Shows 'letter' with max coverage"""
    maxcov = max(vector)
    ind_of_max = [i for i in xrange(len(vector)) if vector[i] == maxcov]
    ind = random.choice(ind_of_max)
    return(ind)

def check_variants(procdata):
    variants = [0, 0, 0, 0]
    for vector in procdata:
        # maxcov = max(vector)
        # max_indices = [i for i in xrange(0, 4) if vector[i] == maxcov]
        # for i in max_indices:
            # variants[i] += 1
        ind = getmax(vector)
        variants[ind] += 1
    summa = sum(variants)
    diff = summa - max(variants)
    sorted_variants = sorted(variants, reverse = True)
    second_letter = sorted_variants[1]
    return([map(str,variants + [summa, diff, second_letter]), getmax(variants)])

def dict_ini(chrom_set, filenum):
    """Auxilary func for initializing dictionary of triangle matrices"""
    sampledict = {}
    for chrom in chrom_set:
        sampledict[chrom] = []
        for i in xrange(filenum):
            sampledict[chrom].append([])
            for j in xrange(filenum):
                sampledict[chrom][i].append(0)
    return(sampledict)

def read_all_files(index_handler, filename_list, outf, chrom_set):
    """Actually, it is the main function. Here we read all files line by line
    and determine whether variable position is covered in particular
    sample or not. Then only covered (!) samples are passed into
    dist_matrix_calculate() function that gives us desired result"""
    filenum = len(filename_list)
    #crossh_dict = dict_ini(chrom_set, filenum)
    #mutualcov_dict = [dict_ini(chrom_set, filenum) for i in tresholds]
    #diffcov_dict = [dict_ini(chrom_set, filenum) for i in tresholds]

    prev_time = time()
    index_file = open(index_handler, 'r')
    files = [open(filename) for filename in filename_list]
    # outfiles = [open(filename + '.freq', 'w') for filename in filename_list]
    # outdict = [dict() for filename in filename_list]
    print 'Num of files: %d' % len(files)
    read_line_or_not = [True] * filenum
    data = ['0'] * 4 * filenum
    lines_of_files = [''] * filenum

    loop_time = time()
    for i, line in enumerate(index_file):
        filenums_for_proc = []
        if (i % 10000) == 0:
            print ('%dth step. Time is %f\n' % (i, (time() - loop_time)))
        chrom, position = line.strip().split()
        for j, file_ in enumerate(files):
            if read_line_or_not[j]:
                lines_of_files[j] = file_.readline().strip().split()
            if len(lines_of_files[j]) > 1:
                position_in_file = lines_of_files[j][1]
                if position == position_in_file:
                    read_line_or_not[j] = True
                    data[j] = lines_of_files[j][2:]
                    filenums_for_proc.append(j)
                else:
                    read_line_or_not[j] = False
        #print lines_of_files
        #print read_line_or_not
        # if position in ['10219','15505','16819','16822','16823']:
        #     print '\n____%s\n' % position
        #     for i, el in enumerate(lines_of_files):
        #         print '%s\t%s' % ('\t'.join(el), read_line_or_not[i])
            
        procdata = [map(int, data[i]) for i in filenums_for_proc]
        num_of_proc_files = len(procdata)
        variants, max_ind = check_variants(procdata)
        # for j, file_ in enumerate(outfiles):
        #     if read_line_or_not[j] and getmax(procdata[j]) != max_ind:
        #         outdict[j][variants[-3:-1]] = outdict[j].setdefault(variants[-3:-1], 0) + 1
        #         file_.write('\t'.join(out_info) + '\n')

        outlist = [chrom, position] + variants
        outf.write('\t'.join(outlist) + '\n')

    for file_ in files:
        file_.close()
    index_file.close()


def main(index_handler, filename_list_handler, chrom_list_handler, prefix):
    #try:
    path = re.match(r'.+(?=/\w)', chrom_list_handler).group()
    #get list of *.tab filenames
    filename_list_file = open(filename_list_handler, 'r')
    filename_list = [line.strip() for line in filename_list_file]
    filename_list_shorter = [re.search('(?<=\/)\w+(?=\.)',el).group()
                             for el in filename_list]
    filename_list_file.close()
    filenum = len(filename_list)

    outf = open('%s/%svarpos_variants.txt' % (path, prefix), 'w')

    #get list of contignames
    chrom_list_file = open(chrom_list_handler, 'r')
    chrom_set = Set([line.strip() for line in chrom_list_file])
    read_all_files(index_handler, filename_list, outf, chrom_set)

    chrom_list_file.close()
    outf.close()

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tBigTab2.py <in: sorted index file> <in: sorted *.tab filename list> <in: contig list>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
