import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq


LETTERS = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}


def get_max_variant(variants):
    """ This function takes a part of
    line in *.snptab file with number of reads
    supporting particular letter and gets letter
    with max support"""
    maxcov = max(variants)
    inds_of_max = [i for i in xrange(0, 4) if variants[i] == maxcov]
    ind = random.choice(inds_of_max)
    return(LETTERS[ind])


class SNP(object):
    """docstring for SNP"""
    def __init__(self, fields):
        id_ = fields[0]
        pos = int(fields[1]) - 1
        variants = [int(field) for field in fields[2:]]
        self.id = id_
        self.snps = {}
        self.snps[pos] = get_max_variant(variants)

    def add(self, fields):
        """Checks whether ids match to each other. If True
        add position and supported letter to dict"""
        id_ = fields[0]
        if id_ == self.id:
            pos = int(fields[1]) - 1
            variants = [int(field) for field in fields[2:]]
            self.snps[pos] = get_max_variant(variants)
            return(True)
        else:
            return(False)


def insert_snp(record, snp_record):
    """Takes Seq object and change some
    nucleotides according to SNP record object"""
    sequence = list(record.seq)
    # newseq = [snp_record.snps[i] if i in snp_record.snps.keys()
    #           else sequence[i] for i in xrange(len(sequence))]
    newseq = [snp_record.snps.get(i, sequence[i])
                        for i in xrange(len(sequence))]
    record.seq = Seq(''.join(newseq))
    return(record)


def get_snp_iter(snptabfilename):
    """ Input: *.snptab filename
        Output: generator yielding SNP-objects
        consisting of contig id and dictionary
        {pos1: letter1, pos2: letter2}"""
    f = open(snptabfilename, 'r')
    for i, line in enumerate(f):
        fields = line.strip().split()
        if i == 0:
            snp_record = SNP(fields)
        else:
            if snp_record.add(fields):
                pass
            else:
                yield(snp_record)
                snp_record = SNP(fields)
    yield(snp_record)


def get_output_fasta_iter(input_seq_iterator, snp_iter):
    """ Input: fasta file iterator and snps iterator
        output: iterator yielding fasta records with
        changes in particular positions """
    for i, record in enumerate(input_seq_iterator):
        if i == 0:
            snp_record = next(snp_iter, 'END')
        if snp_record != 'END':
            if record.id == snp_record.id:
                yield(insert_snp(record, snp_record))
                snp_record = next(snp_iter, 'END')
            else:
                yield(record)
        else:
            yield(record)


def main(inputfastaname, snptabfilename, outputfastaname):
    """General function for editing fasta"""

    input_f = open(inputfastaname, "r")
    input_seq_iterator = SeqIO.parse(input_f, "fasta")
    snp_iter = get_snp_iter(snptabfilename)
    output_fasta_iter = get_output_fasta_iter(input_seq_iterator, snp_iter)

    output_f = open(outputfastaname, "w")
    SeqIO.write(output_fasta_iter, output_f, "fasta")
    output_f.close()
    input_f.close()

if __name__ == '__main__':
    if (len(sys.argv) < 4):
        print 'usage:\tFastaSnpInsert.py <input fasta file> <*.snptab file> <output fasta file>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
