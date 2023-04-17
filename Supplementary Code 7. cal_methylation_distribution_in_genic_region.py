import re
import gzip
import argparse
import numpy as np
from collections import defaultdict

"""
	要求输入的methl_counts, 一条染色体单独一个文件
"""

parser = argparse.ArgumentParser()
parser.add_argument('ingff')
parser.add_argument('methl_counts')
parser.add_argument('--chromosome-length', dest='chromosomeLength', type=int) #单独一个染色体的长度
parser.add_argument('--chromosome')	#染色体的名称
parser.add_argument('--gene-list', dest='listFile') #基因list文件，一个行一个基因，若有该文件，则只统计文件中的基因
parser.add_argument('--type', dest='methyl_type', default='CG')
parser.add_argument('-o', '--outfile', default='genic_region_methylation_level.csv')
args = parser.parse_args()

limit = 2000
interval=100

idRegrex = re.compile(r'ID=([^;]+)')

if args.listFile:
	genes = {}
	with open(args.listFile, 'r') as LF:
		for line in LF:
			gene = line.rstrip().split()[0]
			genes[gene] = True
else:
	genes = None


genePosition = []
with open(args.ingff, 'r') as IG:
	for line in IG:
		if line.startswith('#'):
			continue
		else:
			linelist = line.rstrip().split('\t')
			if linelist[0] != args.chromosome:
				continue

			if linelist[2] == 'mRNA':
				mrna_id = idRegrex.search(linelist[-1]).group(1)
				start, end = list(map(int, linelist[3:5]))
				strand = linelist[6]
				if strand == '+':
					start, end = sorted([start, end])
				else:
					start, end = reversed([start, end])

				if genes == None:
					genePosition.append((start, end))
				else:
					if genes.get(mrna_id, False):
						genePosition.append((start, end))


if args.methl_counts.endswith(('gz', 'GZ')):
	MC = gzip.open(args.methl_counts, 'rt')
else:
	MC = open(args.methl_counts, 'r')

GM = open(args.outfile, 'w')

position_array = np.zeros(args.chromosomeLength, dtype=bool)
methyl_level_array = np.zeros(args.chromosomeLength)

for line in MC:
	linelist = line.rstrip().split('\t')
	methyl_type = linelist[3]
	index = int(linelist[1]) -1
	methyl_num = float(linelist[6])
	unmethyl_num = float(linelist[7])
	
	if methyl_num + unmethyl_num < 5:
		continue

	if methyl_type != args.methyl_type:
		continue

	methyl_level = methyl_num/(methyl_num + unmethyl_num)
	methyl_level_array[index] = methyl_level
	position_array[index] = True

MC.close()

def get_relative_level(methyl_level, c_num):
	if methyl_level == 0 and c_num == 0:
		return 0
	else:
		return methyl_level/c_num

#print(genePosition)
for start, end in genePosition:
	strand = '+' if start < end else '-'
	start, end = sorted([start, end])

	upstream = []
	if start > limit:
		for i in range(start-limit-1,start-1,interval):
			relative_level = get_relative_level(np.sum(methyl_level_array[i:i+interval]), np.sum(position_array[i:i+interval]))
			upstream.append(relative_level)
	else:
		if (start-1)%interval == 0:
			for i in range(int(limit-(start-1)/interval)):
				upstream.append(0)
			for i in range(0, start-1, interval):
				relative_level = get_relative_level(np.sum(methyl_level_array[i:i+interval]), np.sum(position_array[i:i+interval]))
				upstream.append(relative_level)
		else:
			for i in range((limit-(start-1))//interval):
				upstream.append(0)

			remainder = (start-1)%interval
			relative_level = get_relative_level(np.sum(methyl_level_array[0:remainder]), np.sum(position_array[0:remainder]))
			upstream.append(relative_level)

			for i in range(remainder, start-1, interval):
				relative_level = get_relative_level(np.sum(methyl_level_array[i:i+interval]), np.sum(position_array[i:i+interval]))
				upstream.append(relative_level)
	if len(upstream) != 20:
		print(start, end)
		print('length of upstream is ',len(upstream))
	
	genebody = []
	geneInterval = (end-start)//20 + 1
	
	for i in range(19):
		v = start-1 + i*geneInterval
		t = start-1 + (i+1)*geneInterval
		relative_level = get_relative_level(np.sum(methyl_level_array[v:t]), np.sum(position_array[v:t]))
		genebody.append(relative_level)
	
	v = int(start-1+ 19*geneInterval)
	t = end
	relative_level = get_relative_level(np.sum(methyl_level_array[v:t]), np.sum(position_array[v:t]))
	genebody.append(relative_level)
	
	if len(genebody) != 20:
		print(start, end)
		print('length of genebody is ',len(genebody))
	
	downstream = []
	#print('#downstream')
	if args.chromosomeLength-end >= limit:
		for i in range(end, end+limit, interval):
			#print('i',i)
			#print('a',np.sum(methyl_level_array[i:i+interval]))
			#print('b', np.sum(position_array[i:i+interval]))
			relative_level = get_relative_level(np.sum(methyl_level_array[i:i+interval]), np.sum(position_array[i:i+interval]))
			downstream.append(relative_level)
			#print('relative_level', relative_level)
	else:
		for i in range(end, args.chromosomeLength-1, interval):
			j = i + interval if i + interval<= args.chromosomeLength-1 else args.chromosomeLength-1
			relative_level = get_relative_level(np.sum(methyl_level_array[i:j]), np.sum(position_array[i:j]))
			downstream.append(relative_level)
			#print('a',np.sum(methyl_level_array[i:i+interval]))
			#print('b', np.sum(position_array[i:i+interval]))
			#print('relative_level', relative_level)
		
		for i in range((limit-(args.chromosomeLength-end))//interval):
			downstream.append(0)
	
	if len(downstream) != 20:
		print(end, args.chromosomeLength)
		print('length of downstream is ',len(downstream))
	
	genicRegion = upstream + genebody + downstream
	if strand == '-':
		genicRegion = reversed(genicRegion)
	
	GM.write(','.join(map(str,genicRegion)) + '\n')

GM.close()
