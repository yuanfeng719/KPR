#!/usr/bin/python

import sys

if __name__ == '__main__':
	file0 = sys.argv[1]
	with open(file0,'r') as f:
		file=f.read()
	header = file.split('>')[1].split('\n')[0]
	seq = ''.join(file.split('>')[1].split('\n')[1:])
	length = len(seq)
	string = '>' + header + '\n' + seq + '\n'
	if length == 546:
		out = open(header + '_normDLA' , 'w')
		out.write(string)
	elif length == 549:
		out = open(header + '_50XDLA' , 'w')
		out.write(string)
	out.close()

