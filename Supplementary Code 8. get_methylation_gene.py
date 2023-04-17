import re
import gzip
import argparse
import numpy as np
from scipy.stats import binom
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('methl_sites_list')
parser.add_argument('ingff')
parser.add_argument('--chromosome-length', dest='lengthFile', required=True)
parser.add_argument('--type', dest='methyl_type', choices=('CG', 'CHG', 'CHH'), default='CG')
parser.add_argument('--methylation-proportion', dest='methyl_prop', required=True,\
					type=float, help='the proportion of methylated cytosine residues at CG sites across the whole genome')
args = parser.parse_args()

lengthMap = {}
with open(args.lengthFile, 'r') as LF:
	for line in LF:
		chromosome, length = line.rstrip().split()
		lengthMap[chromosome] = int(length)


id_regrex = re.compile(r'ID=([^;]+)')

geneMap = defaultdict(list)
with open(args.ingff, 'r') as IG:
	for line in IG:
		if line.startswith('#'):
			continue
		linelist = line.rstrip().split('\t')
		if linelist[2] == 'mRNA':
			chromosome = linelist[0]
			mrna_id = id_regrex.search(linelist[-1]).group(1)
			start ,end = sorted(map(int, linelist[3:5]))
			geneMap[chromosome].append([mrna_id, start ,end])

with open(args.methl_sites_list, 'r') as MSL:
	for l in MSL:
		methl_sites, chromosome = l.rstrip().split()
		
		if methl_sites.endswith(('gz', 'GZ')):
			MS = gzip.open(methl_sites, 'rt')
		else:
			MS = open(methl_sites, 'r')
		
		cytosine_array = np.zeros(lengthMap[chromosome], dtype=bool)
		methyl_array = np.zeros(lengthMap[chromosome], dtype=bool)
		for line in MS:
			linelist = line.rstrip().split('\t')
			position = int(linelist[1])
			methyl_type = linelist[3]
			coverage = int(linelist[6]) + int(linelist[7])
			
			if methyl_type !=args.methyl_type or coverage < 5:
				continue
			
			cytosine_array[position-1] = True
			if linelist[-1] == 'Y':
				methyl_array[position-1] = True

		for item in geneMap[chromosome]:
			mrna_id, start, end = item
			cytosine_site_num = np.sum(cytosine_array[start-1:end])
			methyl_site_num = np.sum(methyl_array[start-1:end])
			if cytosine_site_num == 0:
				methyl_rate = 0
				p = 1
			else:
				methyl_rate = float(methyl_site_num)/cytosine_site_num
				p = 1 - binom.cdf(methyl_site_num-1, cytosine_site_num, args.methyl_prop)
			print('{}\t{}\t{}\t{}\t{}'.format(mrna_id, methyl_site_num, cytosine_site_num, methyl_rate, p))
