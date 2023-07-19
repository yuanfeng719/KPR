#!/usr/bin/python

import sys

if __name__ == '__main__':
	sam0 = sys.argv[1]
	clustalO_path = sys.argv[2]
	X50 = sys.argv[3]
	with open(sam0,'r') as f:
		sam = f.read()
	with open(X50,'r') as f:
		X50 = f.read()
	lst = sam.split('\n')[:-1]
	X50_lst = X50.split('\n')[:-1]
	print ';'.join(X50_lst)
	out = open(sam0 + '_adjusted','w')
	for i in range(len(lst)):
		rec = lst[i].split('\t')
		contig = rec[2]
		pos = int(rec[3])
		with open(clustalO_path + contig + '.clustalO' ,'r') as f:
			alignment = f.read()
		contig_seq = ''.join(alignment.split('>')[1].split('\n')[1:])
		norm_seq = ''.join(alignment.split('>')[2].split('\n')[1:])
		gap_length = 0
		flag = True # contig longer than ref norm
		while flag and gap_length < len(norm_seq):
			if norm_seq[gap_length] == '-':
				gap_length += 1
			else:
				flag = False
		new_pos = pos - gap_length + 73
		if gap_length == 0: #contig shorter than ref norm
			flag = True
			gap_length = 0
			while flag and gap_length < len(contig_seq):
				if contig_seq[gap_length] == '-':
					gap_length += 1
				else:
					flag = False
			new_pos = pos + gap_length + 73
		rec[3] = str(new_pos)
		print contig + '\t' + str(gap_length)
		if contig in X50_lst:
			new_contig = 'DLA-88*50X'
		else:
			new_contig = 'DLA-88*0XX'
		rec[2] = new_contig
		string = '\t'.join(rec) + '\n'
		out.write(string)
	out.close()

