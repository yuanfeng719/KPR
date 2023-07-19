#!/usr/bin/python

import sys

#########################
def trim_seq_ends(lst,length):
	result = []
	for i in range(len(lst)):
		seq = lst[i]
		if len(seq) > 2*length:
			new_seq = seq[length:-length]
			result.append(new_seq)
		else:
			result.append(seq)
	return result

def remove_values_from_list(the_list, val):
	return [value for value in the_list if value != val]

def merge_seq(seq_lst):
	lst = []
	for i in range(len(seq_lst)):
		lst.append(seq_lst[i])
	# eliminate contained and identical sequences
	for i in range(len(lst)-1):
		seq1 = lst[i]
		for j in range(i+1,len(lst)):
			seq2 = lst[j]
			if seq2 in seq1:
				lst[j] = ''
			elif seq1 in seq2:
				lst[i] = ''
	lst = remove_values_from_list(lst,'')
	return lst

#########################
if __name__ == '__main__':
	fasta = sys.argv[1]
	length = int(sys.argv[2])
	trim = sys.argv[3]
	with open(fasta,'r') as f:
        	file=f.read()
	lst = file.split('>')[1:]
	# data processing
	# create sequence list
	seq_lst = []
	for i in range(len(lst)):
		seq = ''.join(lst[i].split('\n')[1:])
		seq_lst.append(seq)
	#trim sequences in list
	if trim == 'trim':
		seq_lst = trim_seq_ends(seq_lst,length)
	# merge sequence in list
	merged = merge_seq(seq_lst)
	# output files
	out = open(fasta+'_merged_and_renamed','w')
	# write results
	for i in range(len(merged)):
		header = 'Align_' + str(i+1)
		seq = merged[i]
		if len(seq) > 10:
			string = '>' + header + '\n' + seq + '\n'
			out.write(string)
	out.close()

