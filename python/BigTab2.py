import sys
import random
import re
import sets
from time import time
from sets import Set
from collections import defaultdict
import json

def getmax(vector):
	#vector = map(lambda x: int(x), vector)
	maxcov = max(vector)
	ind_of_max = [i for i in xrange(0,4) if vector[i]==maxcov]
	ind = random.sample(ind_of_max,1)
	return(ind[0])



def heterogen(vector):
	sum_squared = (sum(vector))**2
	squares = map(lambda x: x**2, vector)
	if sum_squared!=0:
		prob = (sum_squared - sum(squares))/float(sum_squared)
	else:
		prob = 0
	return prob

def crossheterogen(vector1,vector2):
	try:
		summ = map(lambda x,y: x+y,vector1,vector2)
	except TypeError:
		print vector1
		print vector2
		sys.exit()
	prob = heterogen(summ)
	return prob

def distance(vector1,vector2,treshold):
	sumcov1 = sum(vector1)

	sumcov2 = sum(vector2)
	mutualcov = ((sumcov1+sumcov2) > 0)
	commoncov = (sumcov1 >= treshold) & (sumcov2 >= treshold)
	if commoncov:
		diffcov = (getmax(vector1)!=getmax(vector2))
	else:
		diffcov = False
	return([mutualcov,commoncov,diffcov])

#Initialize dictionary for distance matrices
def dict_ini(chroms,filenum):
	sampledict = {}
	for chrom in chroms:
		sampledict[chrom]=[]
		for i in xrange(filenum):
			sampledict[chrom].append([])
			for j in xrange(filenum):
				sampledict[chrom][i].append(0)
	return(sampledict)

#Cut line and convert to int
def vectorcut(l,n):
	v = l[4*n:4*n+4]
	v = map(lambda x: int(x), v)
	return(v)


#Make pairwise measurement of parameters of line
# line= union of all parameters on position of *.tab files
# return stat parameters for two filelines


# Function for conversion second array element to int
def intconv(x): return([x[0],int(x[1])]+x[2:])

# Function for conversion second array element to str
def strconv(x): return([x[0],str(x[1])]+x[2:])

# Creation of hashtable and list of file pointers
# rewrite!!
# def hashcreate(filenamelist,index,logf):
# 	hashtab = []
# 	files = []
# 	for i,filename in enumerate(filenamelist):
# 		hashtab.append([])
# 		logf.write('Load '+ filename+'\n')
# 		f_temp = open(filename,'r')
# 		tab = [line.strip().split() for line in f_temp]
# 		j= 0
# 		length = len(tab)
# 		for pos in index:
# 			hashtab[i].append(pos==tab[j][0:2])
# 			if (pos==tab[j][0:2]) and (j<(length-1)):
# 				j += 1
# 		f_temp.close()

# 		files.append(open(filename,'r'))#[:-6]+'tab','r'))
# 	return([files,hashtab])

def distMatrixCalculate(indexfilename,treshold,filenamelist,logf,chroms,filenum,path):

	crossh_dict = dict_ini(chroms,filenum)
	mutualcov_dict = dict_ini(chroms,filenum)
	commoncov_dict = dict_ini(chroms,filenum)
	diffcov_dict = dict_ini(chroms,filenum)


	t = time()
	indexFile = open(indexfilename,'r')
	logf.write('Creating stat tables: \n')
	filenum = len(filenamelist)
	files = [open(filename) for filename in filenamelist]
	indicator = [1]*filenum
	myLine = ['0']*filenum
	heterogenInd = [0.0]*filenum
	totalcov = [0]*filenum
	goodcov = [0]*filenum
	bigLine = ['0']*(4*filenum)
	els = [[]]*filenum
	N = 0
	for k,line in enumerate(indexFile):
		if (k % 100000) == 0: logf.write(str(k*100000)+'th step\n')
		pos = line.strip().split()
		if k == 0:
			prev = pos[0]
			f = open(path+'/bigtab.tab','w')
			f2 = open(path+'/heterogen.list','w')
		if prev != pos[0]:
			for l,num in enumerate(heterogenInd):
				f2.write('\t'.join(map(lambda x: str(x),[heterogenInd[l],totalcov[l],goodcov[l]]))+'\n')
			f.close()
			f2.close()
			prev = pos[0]
			N += 1
			f = open(path+'/bigtab'+str(N)+'.tab','w')
			f2 = open(path+'/heterogen'+str(N)+'.list','w')
			heterogenInd = [0.0]*filenum
			totalcov = [0]*filenum
			goodcov = [0]*filenum
			#totalcov = [0]*filenum
			#goodcov = [0]*filenum
		t2 = time()
		for i in xrange(filenum):
			if indicator[i] == 1:
				el = files[i].readline().strip().split()
				els[i] = el
				if el[:2]!=pos:
					indicator[i] = 0
					myLine[i] = '0'
					#bigLine[i*4:(4*i+4)] = ['0']*4
				else:
					totalcov[i] += 1
					vector = map(lambda x: int(x), el[2:])
					if sum(vector)>=treshold:
						heterogenInd[i] += heterogen(vector)
						myLine[i] = str(getmax(vector)+1)
						goodcov[i] +=1
					else:
						myLine[i]='0'

					#bigLine[i*4:(4*i+4)] = el[2:]
			else:
				el = els[i]
				if el[:2]!=pos:
					myLine[i] = '0'
					#bigLine[i*4:(4*i+4)] = ['0']*4
				else:
					indicator[i]=1
					vector = map(lambda x: int(x), el[2:])
					totalcov[i] += 1
					if sum(vector)>=treshold:
						heterogenInd[i] += heterogen(vector)
						myLine[i] = str(getmax(vector)+1)
						goodcov[i] += 1 
					else:
						myLine[i]='0'
					#bigLine[i*4:(4*i+4)] = el[2:]
		f.write('\t'.join([pos[1]]+myLine)+'\n')
		logf.write('BigLine is done! T='+str(time()-t2)+'\n')
		t2 = time()
	for l,num in enumerate(heterogenInd):
		f2.write('\t'.join(map(lambda x: str(x),[heterogenInd[l],totalcov[l],goodcov[l]]))+'\n')
		prev = pos[0]

		###################################################
		# for i in xrange(filenum-1):
		# 	for j in xrange(i+1,filenum):
		# 		print i,j
		# 		#vector1 = vectorcut(bigLine,i)
		# 		#vector2 = vectorcut(bigLine,j)
		# 		vector1 = [0,0,3,0]
		# 		vector2=[0,0,2,0]
		# 		crossh = crossheterogen(vector1,vector2)
		# 		dist = distance(vector1,vector2,treshold)
		# 		crossh_dict[pos[0]][i][j] += crossh
		# 		mutualcov_dict[pos[0]][i][j] += dist[0]
		# 		commoncov_dict[pos[0]][i][j] += dist[1]
		# 		diffcov_dict[pos[0]][i][j] += dist[2]
		##################################################
		# logf.write('Stat is done! T='+str(time()-t2)+'\n')
		
	#return([crossh_dict,mutualcov_dict,commoncov_dict,diffcov_dict])

		
				


	# temp = hashcreate(filenamelist,index,logf)
	# files = temp[0]
	# hashtab = temp[1]
	# logf.write('Done! T='+str(time()-t)+'\n')

	
	# for k, pos in enumerate(index):
	# 	if (k % 10000) == 0: logf.write(str(k*10000)+'th step\n')
	# 	line = strconv(pos)
	# 	for i,f in enumerate(files):
	# 		if hashtab[i][k]:
	# 			line = line + f.readline().strip().split()[2:6]
	# 		else:
	# 			line = line + ['0']*4
	# 	M = measure(line,filenum,treshold)
	# 	crossh_dict[line[0]][i][j] += M[0]
	# 	mutualcov_dict[line[0]][i][j] += M[1]
	# 	commoncov_dict[line[0]][i][j] += M[2]
	# 	diffcov_dict[line[0]][i][j] += M[3]

	##################################################
	# logf.write('Done! T='+str(time()-t)+'\n')
	
	# logf.write('Writing dumps\n')
	# t = time()

	# outf = open('json.dumps.tresh'+ str(treshold) + '.txt','w')
	# outf.write(json.dumps(crossh_dict)+'\n')
	# outf.write(json.dumps(mutualcov_dict)+'\n')
	# outf.write(json.dumps(commoncov_dict)+'\n')
	# outf.write(json.dumps(diffcov_dict)+'\n')
	# logf.write('Done! T='+str(time()-t)+'\n')
	##################################################


def main(indexfilename,filenamelistfile,chromlistfilename,strtreshold,path):
	logf = open('logfile.log','w')
	#try:
	treshold = int(strtreshold)
	
	#get list of *.tab filenames
	listFile = open(filenamelistfile,'r')
	filenamelist = [line.strip() for line in listFile]
	listFile.close()
	filenum = len(filenamelist)
	logf.write ('File number: '+str(filenum)+'\n')
	
	#get list of contignames 
	chromlistFile = open(chromlistfilename,'r')
	chroms = Set([line.strip() for line in chromlistFile])
	chromlistFile.close()

	#t = time()
	
	# logf.write ('Reading indexfile\n')
	# #rewrite!!
	# #for line in open(indexfilename,'r'):
	# indexFile = open(indexfilename,'r')
	# index = [line.strip().split() for line in indexFile]
	# indexFile.close()
	# logf.write('Done! T='+str(time()-t)+'\n')
	# t=time()
		
	# logf.write('\nSstatistics\n')
	distMatrixCalculate(indexfilename,treshold,filenamelist,logf,chroms,filenum,path)
	# except:
	# 	error = sys.exc_info()
	# 	print(error)
	# 	#logf.write(error)
	# finally:
	logf.close()
	
if __name__ == '__main__':
	flag = 'none'
	if (len(sys.argv) < 2):
		print 'usage:\tBigTab2.py <in: sorted index file> <in: sorted *.tab filename list> <in: contig list> <in: treshold>'
		sys.exit(0)
	else:
		main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])





