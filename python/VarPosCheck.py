import sys
import random
import re
import sets


def main(handler,varposfilename):
	def intconv(x):
		try:
			return([x[0],int(x[1])]+x[2:])
		except IndexError:
			pass
	
		

	f_snptab_list = open(handler,'r')
	f_varpos = open(varposfilename,'w')
	union = sets.Set({})
	for fname in f_snptab_list:
		f_snptab = open(fname.strip(),'r')
		snptab = {'\t'.join(line.split()[0:2]) for line in f_snptab}
		snptab_set = sets.Set(snptab)
		print(len(snptab_set))
		union = union | snptab_set
		print(len(union))

	varpos = [ intconv(line.strip().split()) for line in union]
	varpos = sorted(varpos)
	l = len(varpos)
	for i,line in enumerate(varpos):
		if i!=(l-1):
			f_varpos.write(line[0]+'\t'+str(line[1])+'\n')
		else:
			f_varpos.write(line[0]+'\t'+str(line[1]))


if __name__ == '__main__':
	flag = 'none'
	if (len(sys.argv) < 2):
		print 'usage:\tVarPosCheck.py <in: snptab list> <out: varposfilename>'
		sys.exit(0)
	else:
		main(sys.argv[1],sys.argv[2])
