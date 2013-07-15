import sys

# def get_chrom_and_pos(filename):
# 	f = open(filename,'r')
# 	for line in f:
# 		words = line.strip().split()
# 		yield([])
data = {}
filelist = (line.strip() for line in open(sys.argv[1],'r'))
for filename in filelist:
	chroms_and_pos = (line.strip().split()[:2] for line in open(filename))	
	for chrom, pos in chroms_and_pos:
		data.setdefault(chrom,{})
		data[chrom][pos] = data.setdefault(chrom, {}).setdefault(pos, 0) + 1

outf = open(sys.argv[2],'w')
for chrom in data.keys():
	for pos in data[chrom].keys():
		outf.write('%s\t%s\t%d\n' % (chrom, pos, data[chrom][pos]))