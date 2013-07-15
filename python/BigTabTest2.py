import sys
import random
import re
from time import time
from sets import Set
import json


def getmax(vector):
    """Shows 'letter' with max coverage"""
    maxcov = max(vector)
    ind_of_max = [i for i in xrange(0, 4) if vector[i] == maxcov]
    ind = random.sample(ind_of_max, 1)
    return(ind[0])


def check_variants(procdata):
    variants = [0, 0, 0, 0]
    for vector in procdata:
        maxcov = max(vector)
        max_indices = [i for i in xrange(0, 4) if vector[i] == maxcov]
        for i in max_indices:
            variants[i] += 1
    variants = sorted(variants)
    max_alt = variants[2]
    sum_of_alt = sum(variants[:3])
    return([max_alt, sum_of_alt])    

def heterogen(vector):
    """Heterogenity meadure function. It should be noted
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


def crossheterogen(vector1, vector2):
    """Computes crossheterogenity distance between two vectors"""
    try:
        summ = map(lambda x, y: x + y, vector1, vector2)
    except TypeError:
        print vector1
        print vector2
        sys.exit()
    prob = heterogen(summ)
    return prob


def distance(vector1, vector2):
    """Function returns whether letters at same position with maximum
    coverage equal or not"""
    diffcov = (getmax(vector1) != getmax(vector2))
    return(diffcov)


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





def read_all_files(index_handler, filename_list, logf, chrom_set, path,tresholds, max_alt_treshold,sum_alt_treshold):
    """Actually, it is the main function. Here we read all files line by line
    and determine whether variable position is covered in particular
    sample or not. Then only covered (!) samples are passed into
    dist_matrix_calculate() function that gives us desired result"""
    filenum = len(filename_list)
    crossh_dict = dict_ini(chrom_set, filenum)
    mutualcov_dict = [dict_ini(chrom_set, filenum) for i in tresholds]
    diffcov_dict = [dict_ini(chrom_set, filenum) for i in tresholds]

    def dist_matrix_calculate(filenums_for_proc, data, chrom):
        procdata = [map(int, data[i]) for i in filenums_for_proc]
        num_of_proc_files = len(procdata)
        max_alt,sum_of_alt = check_variants(procdata)

        max_alt_cond = (max_alt >= int(max_alt_treshold))
        sum_alt_cond = (sum_of_alt >= int(sum_alt_treshold))
        if ((num_of_proc_files > 1) & max_alt_cond & sum_alt_cond):
            for i in xrange(num_of_proc_files - 1):
                for j in xrange(i + 1, num_of_proc_files):
                    ii = filenums_for_proc[i]
                    jj = filenums_for_proc[j]
                    crossh_dict[chrom][ii][jj] += crossheterogen(procdata[i],
                                                                 procdata[j])
                    for k, treshold in enumerate(tresholds):
                        if ((sum(procdata[i]) >= treshold) & 
                            (sum(procdata[j]) >= treshold)):
                            mutualcov_dict[k][chrom][ii][jj] += 1
                            diffcov_dict[k][chrom][ii][jj] += distance(procdata[i],
                                                            procdata[j])
    prev_time = time()
    index_file = open(index_handler, 'r')
    logf.write('Creating stat tables: \n')
    files = [open(filename) for filename in filename_list]
    read_line_or_not = [True] * filenum
    data = ['0'] * 4 * filenum
    lines_of_files = [''] * filenum

    loop_time = time()
    for i, line in enumerate(index_file):
        filenums_for_proc = []
        if (i % 10000) == 0:
            print ('%dth step. Time is %f\n' % (i, (time() - loop_time)))
            logf.write('%dth step. Time is %f\n' % (i, (time() - loop_time)))
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
        dist_matrix_calculate(filenums_for_proc, data, chrom)

    for file_ in files:
        file_.close()
    index_file.close()
    logf.write('Done! Total time is %.1f' % (time() - prev_time))

    outf = open('%s/json.distance_matrices_2.txt' % path, 'a')
    outf.write(json.dumps(crossh_dict) + '\n')
    outf.write(json.dumps(mutualcov_dict) + '\n')
    outf.write(json.dumps(diffcov_dict) + '\n')
    outf.close()


def main(index_handler, filename_list_handler, chrom_list_handler, max_alt_treshold, sum_alt_treshold, tresholds):
    logf = open('logfile.log', 'w')
    #try:
    path = re.match(r'.+(?=/\w)', chrom_list_handler).group()
    #get list of *.tab filenames
    filename_list_file = open(filename_list_handler, 'r')
    filename_list = [line.strip() for line in filename_list_file]
    filename_list_shorter = [re.search('(?<=\/)\w+(?=\.)',el).group()
                             for el in filename_list]
    filename_list_file.close()
    filenum = len(filename_list)
    logf.write('Number of files: %d\n' % filenum)

    outf = open('%s/json.distance_matrices_2.txt' % path, 'w')
    message = """
    List contains names of samples specified by
    tablistfile: %s
    -------------------------------------------------------------
    First dictionary means crossheterogenity measure.
    This number is defined for the positions that are covered in
    both samples in pair. That's wrong.
    Structure of this dict:
    {chrom:[[] * num_of_samples] * num_of_samples]}
    In this case num_of_samples = %d
    --------------------------------------------------------------
    Second dictionary means number of positions in pair of samples
    that have coverage > treshold in both samples.
    Structure of this dict:
    [{chrom:[[] * num_of_samples] * num_of_samples]}]*num_of_tresholds
    Tresholds are %s
    --------------------------------------------------------------
    Third dictionary means number of positions with coverage > treshold
    where major letter are not similar.
    --------------------------------------------------------------
    """ % (filename_list_handler, filenum, str(tresholds))

    outf.write(message + '\n')
    outf.write(json.dumps(filename_list_shorter) + '\n')
    outf.close()

    #get list of contignames
    chrom_list_file = open(chrom_list_handler, 'r')
    chrom_set = Set([line.strip() for line in chrom_list_file])
    tresholds = [int(el) for el in tresholds]

    read_all_files(index_handler, filename_list, logf, chrom_set, path,tresholds, max_alt_treshold,sum_alt_treshold)

    chrom_list_file.close()
    logf.close()

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tBigTab2.py <in: sorted index file> <in: sorted *.tab filename list> <in: contig list>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:])
