#!/usr/bin/python

import sys
import re

if __name__ == '__main__':
	fa = sys.argv[1]
	with open(fa,'r') as f:
		file = f.read()
	lst = file.split('>')[1:]
	dic = {}
	for i in range(len(lst)):
		name = lst[i].split('\n')[0]
		seq = ''.join(lst[i].split('\n')[1:])
		dic[name] = seq
	seqLst = dic.values()
	poly = {}
	for seq in seqLst:
		for i in range(len(seq)):
			if i not in poly.keys():
				poly[i] = []
			nuc = seq[i]
			if (nuc != '-') and (nuc not in poly[i]):
				poly[i].append(nuc)
	out = open(fa.split('.')[0] + '_PolymorphicSites','w')
	out2 = open(fa.split('.')[0] + '_ContigLinkages','w')
	siteLst = poly.keys()
	siteLst.sort()
	real = []
	for i in range(len(siteLst)):
		nucLst = poly[siteLst[i]]
		nucLst.sort()
		if len(nucLst) >= 2:
			real.append(siteLst[i])
			string = str(siteLst[i]) + '\t' + '\t'.join(nucLst) + '\n'
			out.write(string)
	contigLst = dic.keys()
	for contig in contigLst:
		seq = dic[contig]
		result = []
		for site in real:
			nuc = seq[site] 
			result.append(str(site) + nuc)
		string = '>' + contig + '\n' + '\t'.join(result) + '\n'
		out2.write(string)
	out.close()
	out2.close()

