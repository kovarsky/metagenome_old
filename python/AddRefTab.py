import re
import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq
import os

#LETTERS


def main():
    varpos_filename = sys.argv[1]
    directory = os.path.dirname(varpos_filename)
    varpos_data = open(varpos_filename, 'r') 
    input_fasta = open('%s/reference.fasta' % directory, "r")
    input_seq_iterator = SeqIO.parse(input_fasta, "fasta")

    seqs = {seq.id: seq for seq in input_seq_iterator}
    outf = open('%s/reference.list' % directory, 'w')
    for line in varpos_data:
        chrom, pos = line.strip().split()
        outf.write(seqs[chrom][int(pos) - 1] + '\n')

    outf.close()
    varpos_data.close()
    input_fasta.close()


main()    

