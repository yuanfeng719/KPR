#!/usr/bin/python

import sys
import re

#############################################
def sequence_dic(seq):
        dic = {}
        for i in range(len(seq)):
                name = seq[i].split('\n')[0]
                sequence = ''.join(seq[i].split('\n')[1:])
                dic[name] = sequence
        return dic

def write_file(dic,lst,gp):
	out = open(gp + '.txt','w')
	for i in range(len(lst)):
		name = lst[i]
		seq = dic[name]
		string = '>' + name + '\n' + seq + '\n'
		out.write(string)
	out.close()

############################################3
if __name__ == '__main__':
	file0 = sys.argv[1]
	seq = sys.argv[2]
	with open(file0,'r') as f:
        	file=f.read()
	with open(seq,'r') as f:
        	seq = f.read()
	lst = file.split('\n')[1:-1]
	seq = seq.split('>')[1:]
	# data processing
	allele_dic = {}
	for i in range(len(lst)):
		rec = re.sub(' +','\t',lst[i])
		allele = rec.split('\t')[0]
		allele_dic[i] = allele
	X50 = []
	D88norm = []
	D12 = []
	D64 = []
	for i in range(len(lst)):
		rec = re.sub(' +','\t',lst[i])
		allele = rec.split('\t')[0]
		if 'Align' in allele:
			dist = rec.split('\t')[1:]
			dist = map(float,dist)
			flag = True
			while flag:
				idx = dist.index(max(dist))
				best_allele = allele_dic[idx]
				if 'DLA-88*5' in best_allele or 'DLA-88*J' in best_allele:
					flag = False
					X50.append(allele)
					print allele + '\tDLA88-50X_contigs\t' + best_allele
				elif 'DLA-88' in best_allele:
                                	flag = False
                                	D88norm.append(allele)
					print allele + '\tDLA88-norm_contigs\t' + best_allele
				elif 'DLA-12' in best_allele:
					flag = False
					D12.append(allele)
					print allele + '\tDLA12_contigs\t' + best_allele
				elif 'DLA-64' in best_allele:
					flag = False
					D64.append(allele)
					print allele + '\tDLA64_contigs\t' + best_allele
				else:
					dist[idx] = 0.0
	# output files
	dic = sequence_dic(seq)
	write_file(dic,X50,'DLA88-50X_contigs')
	write_file(dic,D88norm,'DLA88-norm_contigs')
	write_file(dic,D12,'DLA12_contigs')
	write_file(dic,D64,'DLA64_contigs')

