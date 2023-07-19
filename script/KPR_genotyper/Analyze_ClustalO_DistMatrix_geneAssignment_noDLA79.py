#!/usr/bin/python

import sys
import re

file0 = sys.argv[1]

with open(file0,'r') as f:
	file=f.read()

lst = file.split('\n')[1:-1]

out = open(file0 + '_Summary','w')

allele_dic = {}
for i in range(len(lst)):
	rec = re.sub(' +','\t',lst[i])
	allele = rec.split('\t')[0]
	allele_dic[i] = allele

for i in range(len(lst)):
	rec = re.sub(' +','\t',lst[i])
	allele = rec.split('\t')[0]
	if ('Align' in allele) or ('contig' in allele):
		dist = rec.split('\t')[1:]
		dist = map(float,dist)
		flag = True
		while flag:
			idx = dist.index(max(dist))
			best_allele = allele_dic[idx]
			if 'DLA-88' in best_allele:
				flag = False
				string = allele + '\t' + 'DLA-88' + '\t' + str(max(dist)) + '\n'
			elif 'DLA-12' in best_allele:
                                flag = False
                                string = allele + '\t' + 'DLA-12' + '\t' + str(max(dist)) + '\n'
			elif 'DLA-64' in best_allele:
                                flag = False
                                string = allele + '\t' + 'DLA-64' + '\t' + str(max(dist)) + '\n'
			elif 'HLA-A' in best_allele:
				flag = False
                                string = allele + '\t' + 'HLA-A' + '\t' + str(max(dist)) + '\n'
			elif 'HLA-B' in best_allele:
                                flag = False
                                string = allele + '\t' + 'HLA-B' + '\t' + str(max(dist)) + '\n'
			elif 'HLA-C' in best_allele:
                                flag = False
                                string = allele + '\t' + 'HLA-C' + '\t' + str(max(dist)) + '\n'
			else:
				dist[idx] = 0.0
		out.write(string)

out.close()


