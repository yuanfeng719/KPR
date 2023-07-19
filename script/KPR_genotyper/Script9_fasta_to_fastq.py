#!/usr/bin/python

import sys
import re

###################################
# search characters other than atcg
def special_match(strg, search=re.compile(r'[^ATGC]').search):
	return not bool(search(strg))

def complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	bases = list(seq) 
	bases = [complement[base] for base in bases] 
	return ''.join(bases)

def reverse_complement(s):
	return complement(s[::-1])

def intersection(lst1, lst2): 
	lst3 = [value for value in lst1 if value in lst2] 
	return lst3 

##################################
if __name__ == '__main__':
	fasta_fwd = sys.argv[1]
	fasta_rev = sys.argv[2]
	qual = sys.argv[3]
	with open(fasta_fwd,'r') as f:
        	file_fwd=f.read()
	with open(fasta_rev,'r') as f:
        	file_rev = f.read()
	output_fwd = fasta_fwd.split('.')[0] + '.fq'
	output_rev = fasta_rev.split('.')[0] + '.fq'
	out_fwd = open(output_fwd,'w')
	out_rev = open(output_rev,'w')
	lst_fwd = file_fwd.split('>')[1:]
	lst_rev = file_rev.split('>')[1:]
	# data processing
	header_fwd = {}
	for i in range(len(lst_fwd)):
        	header = lst_fwd[i].split('\n')[0]
        	if len(lst_fwd[i].split('\n')) > 1:
                	seq = lst_fwd[i].split('\n')[1]
                	if special_match(seq):
                        	header_fwd[header] = seq
	header_rev = {}
	for i in range(len(lst_rev)):
        	header = lst_rev[i].split('\n')[0]
        	if len(lst_rev[i].split('\n')) > 1:
                	seq = lst_rev[i].split('\n')[1]
                	if special_match(seq):
                        	new_seq = reverse_complement(seq)
                        	header_rev[header] = new_seq
	shared = intersection(header_fwd.keys(),header_rev.keys())
	for i in range(len(shared)):
		header = shared[i]
		seq_fwd = header_fwd[header]
		seq_rev = header_rev[header]
		qual_fwd = qual * len(seq_fwd)
		qual_rev = qual * len(seq_rev)
		string = '@' + header + '\n' + seq_fwd + '\n' + '+' + '\n' + qual_fwd + '\n'
		out_fwd.write(string)
		string = '@' + header + '\n' + seq_rev + '\n' + '+' + '\n' + qual_rev + '\n'
		out_rev.write(string)
	out_fwd.close()
	out_rev.close()

