import sys, os

path = sys.argv[1]

class FastaBunch(object):
    """docstring for Fastas"""
    def __init__(self, names):
        self.fastas = [Fasta(name) for name in names]

    def add(self, letters):
        if all(base!='-' for base in letters):
            if len(set(letters)) > 1:
                for i, fasta in enumerate(self.fastas):
                    fasta.add(letters[i])

    def close(self):
        for fasta in self.fastas:
            fasta.close()


class Fasta(object):
    """docstring for Fasta"""
    def __init__(self, name):
        self.file = open(name, 'w')
        self.file.write('>%s\n' % name)
        self.line = ''

    def add(self, base):
        self.line = self.line + base
        if len(self.line) == 70:
            self.file.write('%s\n' % self.line)
            self.line = ''

    def close(self):
        self.file.write('%s\n' % self.line)
        self.file.close()
        

def main():
    numfile = open('%s/nums.list' % path, 'r')
    nums = [(int(line.strip()) - 1) for line in numfile]
    datafile = open('%s/snp_chars_alt_4' % path, 'r')
    pileupsfile = open('%s/pileups.list' % path, 'r')
    names = ['%s.fasta' % os.path.splitext(line.strip())[0] 
                for i, line in enumerate(pileupsfile) if (i in nums)]
    print names
    fastas = FastaBunch(names)
    for line in datafile:
        letters = [line[num] for num in nums]
        fastas.add(letters)
    fastas.close()


if __name__ == '__main__':
    main()