#!/usr/bin/python

import sys
from os import path

if __name__ == '__main__':
	file0 = sys.argv[1]
	if path.isfile(file0):
		with open(file0,'r') as f:
			file = f.read()
		lst = file.split('\n')[:-1]
		# output files
		out = open(file0 + '_newAllele.fa','w')
		# data processing
		for i in range(len(lst)):
			info = lst[i].split('\t')
			alleleName = info[1]
			allele = info[3]
			gp = info[4]
			seq = info[5]
			if allele == '' and gp == '':
				string = '>' + alleleName + '\n' + seq + '\n'
				out.write(string)
		out.close()


