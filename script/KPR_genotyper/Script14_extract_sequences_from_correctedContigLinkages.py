#!/usr/bin/python

import sys

if __name__ == '__main__':
	correctedContigLinkages0 = sys.argv[1]
	seq0 = sys.argv[2]
	# read in
	with open(correctedContigLinkages0,'r') as f:
		correctedContigLinkages = f.read()
	with open(seq0,'r') as f:
		seq = f.read()
	# get correct contig header
	info = correctedContigLinkages.split('>')
	contigs = []
	if len(info) > 1:
		lst = info[1:]
		for i in range(len(lst)):
			header = lst[i].split('\n')[0]
			contigs.append(header)
	# get sequences
	out = open(seq0 + '_validated','w')
	if len(seq.split('>')) > 1:
		lst = seq.split('>')[1:]
		for i in range(len(lst)):
			header = lst[i].split('\n')[0]
			if header in contigs:
				string = '>' + lst[i]
				out.write(string)
	out.close()	

