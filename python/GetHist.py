import sys
import os

plistname = sys.argv[1]

def main():
	pnames = (line.strip() for line in open(plistname, 'r'))
	for pname in pnames:
		f = open(pname, 'r')
		directory = os.path.dirname(pname)
		outname = (os.path.splitext(pname)[0] + '.hist')
		counts = {}
		for line in f:
			try:
				if line != '':
					count = int(line.split('\t')[3])
					counts[count] = counts.setdefault(count, 0) + 1
				else:
					pass
			except IndexError:
				print(line.split('\t'))
		outf = open(outname, 'w')
		for item in counts.items():
			outf.write('%d\t%d\n' % item)
		outf.close()
		f.close()




if __name__ == '__main__':
	main()
