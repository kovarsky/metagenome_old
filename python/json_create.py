import sys
#import re

def main(handler, filename):
    data_structures = []
    jsonfile = open(handler, 'r')
    data_structures = [line.strip() for line in jsonfile if ((line[0] != ' ') & (len(line) > 1))]
    #print line[:10]
    jsonfile.close()
    newjsonfile = open(filename,'w')
    newjsonfile.write('[%s]\n' % ','.join(data_structures))
    newjsonfile.close()


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print 'usage:\tBigTab2.py <in: raw jsonfile> <output filename>'
        sys.exit(0)
    else:
        main(sys.argv[1], sys.argv[2])
