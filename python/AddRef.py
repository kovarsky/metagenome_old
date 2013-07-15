import re
import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq
import os

#LETTERS


def main():
    inputfastaname = '/data3/bio/metagenome/reference/HMP_2012_02/HMP_2012_02.fasta'
    input_f = open(inputfastaname, "r")
    input_seq_iterator = SeqIO.parse(input_f, "fasta")

    cont_dict = {}
    contig_filelist = open(sys.argv[1], 'r')
    for line in contig_filelist:
        fname = line.strip()
        contigs = open(fname,'r')
        directory = os.path.dirname(fname)
        f = open('%s/reference.fasta' % directory,'w')
        for contig in contigs:
            cont_dict[contig.strip()] = f

    for seq in input_seq_iterator:
        if seq.id in cont_dict.keys():
            SeqIO.write([seq], cont_dict[seq.id], "fasta")
    input_f.close()            


main()    

