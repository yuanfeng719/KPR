#!/usr/bin/python

import sys

if __name__ == '__main__':
	file0 = sys.argv[1]
	with open(file0,'r') as f:
		file = f.read()
	header = file.split('>')[1].split('\n')[0]
	seq = ''.join(file.split('>')[1].split('\n')[1:])
	if 'N' not in seq and (len(seq) >= 546):
		print header


