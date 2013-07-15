import sys
import json
import numpy

print sys.argv[1]
path = sys.argv[1]
distance_f = open(path + 'distance.json','r')
data = json.loads(distance_f.readline().strip())
distance_f.close()

def summarize(matrix_dict):
	values = (numpy.matrix(value) for key, value in matrix_dict.iteritems())
	res = sum(values)
	return(res.tolist())

names = data[0]
crosshet = {'chrom': summarize(data[1])}
comm_pos = {'chrom': summarize(data[2][0])}
diff_pos = {'chrom': summarize(data[3][0])}

outdata = [names, crosshet, [comm_pos], [diff_pos]]
outline = json.dumps(outdata)
distance_f = open(path + 'distance.json','w')
distance_f.write(outline + '\n')
distance_f.close()


stat_f = open(path + 'stat.json', 'r')
data = json.loads(stat_f.readline().strip())
stat_f.close()

cover = {'chrom': summarize(data[0][0])[0]}
hetero = {'chrom': summarize(data[1])[0]}
outdata = [[cover], hetero]
outline = json.dumps(outdata)
stat_f = open(path + 'stat.json', 'w')
stat_f.write(outline + '\n')
stat_f.close()
