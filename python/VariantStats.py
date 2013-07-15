import sys
import re


LETTERS = ['A', 'T', 'G', 'C']
directory = sys.argv[1]
suffix = sys.argv[2]
treshold = int(sys.argv[3])
occtreshold = int(sys.argv[4])
f = open('%s/var_sites_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'r')
f2 = open('%s/mutual_readprofile_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'r')
outf = open('%s/snp_counts_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')
outf2 = open('%s/read_counts_%s_%d_%d' % (directory, suffix, treshold, occtreshold), 'w')

for line in f:
	pos, variants = line.strip().split()
	var_counts = [variants.count(letter) for letter in LETTERS]
	str_var_counts = '\t'.join(str(v) for v in var_counts)
	summ = sum(var_counts)
	snpnum = summ - max(var_counts)
	outf.write('%s\t%s\t%s\t%s\n' % (pos, str_var_counts, summ, snpnum))

occ_dict = {}
for line in f2:
	cover = int(line.strip().split()[1])
	occ_dict[cover] = occ_dict.setdefault(cover,0) + 1

for key, val in occ_dict.items():
	outf2.write('%d\t%d\n' % (key, val))

f.close()
outf.close()

f2.close()
outf2.close()
