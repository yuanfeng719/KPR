#!/usr/bin/python

import sys

if __name__ == '__main__':
	contigs = sys.argv[1]
	tag = sys.argv[2]
	with open(contigs,'r') as f:
		file=f.read()
	lst = file.split('>')[1:]
	# output files
	out = open(contigs + '_renamed','w')
	for i in range(len(lst)):
		name = tag + '_' + str(i+1)
		seq = ''.join(lst[i].split('\n')[1:])
		string = '>' + name + '\n' + seq + '\n'
		out.write(string)
	out.close()

