#!/usr/bin/python

import sys
from os import path

if __name__ == '__main__':
	file0 = sys.argv[1]
	geneAssignment = sys.argv[2]
	with open(file0,'r') as f:
		file = f.read()
	lst = file.split('\n')[:-1]
	allele_dic = {}
	if path.isfile(geneAssignment):
		with open(geneAssignment,'r') as f:
			geneAssignment = f.read()
		geneAssignment = geneAssignment.split('\n')[:-1]
		for i in range(len(geneAssignment)):
			allele = geneAssignment[i].split('\t')[0]
			gene = geneAssignment[i].split('\t')[1]
			allele_dic[allele] = gene
	out = open(file0.split('.')[0] + '_withGeneAssignment.txt','w')
	# adjust total
        totalPerc = 0.0
        for i in range(len(lst)):
                perc = float(lst[i].split('\t')[2])
                totalPerc += perc
        if (totalPerc == 0.0) and (len(lst) > 0):
                indPerc = 1.0 / len(lst)
        else:
                indPerc = 0.0
        for i in range(len(lst)):
		info = lst[i].split('\t')
		alleleName = info[1]
                perc = float(info[2])
                if totalPerc == 0.0:
                        perc = indPerc
		allele = info[3]
		gp = info[4]
		seq = info[5]
		if alleleName in allele_dic.keys():
			gene = allele_dic[alleleName]
		elif allele != '':
			gene = allele.split('*')[0]
		elif gp != '':
			gene = gp.split('*')[0]
		else:
			gene = ''
		if ('50X' not in alleleName) and ('Norm' not in alleleName):
			if len(seq) > 182:
				tag = '50X'
			else:
				tag = 'Norm'
		alleleName = tag + alleleName
		string = info[0] + '\t' + alleleName + '\t' + str(perc) + '\t' + '\t'.join(info[3:5]) + '\t' + gene + '\t' + '\t'.join(info[5:]) + '\n'
		out.write(string)
	out.close()
