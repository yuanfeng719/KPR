#!/usr/bin/python

import sys

#################################################
def Index_pairing(query,ref):
	dic = {}
	idxQ = 0
	idxR = 0
	for i in range(len(query)):
		if query[i] != '-':
			dic[idxQ] = idxR
			idxQ += 1
		idxR += 1
	return (dic)

def reverse(string): 
	string = string[::-1] 
	return string

# replace '-' from 2 tails with 'N'
def replace_tail_char(seq,orig,new):
	counter = 0
	while counter < len(seq) and seq[counter] == orig:
		counter += 1
	new_seq = new * counter + seq[counter:]
	seq = reverse(new_seq)
	counter = 0
        while counter < len(seq) and seq[counter] == orig:
                counter += 1
        new_seq = new * counter + seq[counter:]
        seq = reverse(new_seq)
	return (seq)
		
#################################################
if __name__ == '__main__':
	file0 = sys.argv[1]
	with open(file0,'r') as f:
        	file = f.read()
	lst = file.split('>')[1:]
	# data processing
	assembly = {}
	for i in range(len(lst)):
        	allele = lst[i].split('\n')[0]
        	seq = ''.join(lst[i].split('\n')[1:])
        	if 'DLA' not in allele:
                	assembly[allele] = seq
        	elif 'DLA-88*5' in allele:
                	D88_50X = seq
       		else:
                	D88_norm = seq
	# output files
	out = open(file0 + '_Exon2N3','w')
	patterns = assembly.keys()
	for i in range(len(patterns)):
        	dic_50X = Index_pairing(D88_50X,assembly[patterns[i]])
        	dic_norm = Index_pairing(D88_norm,assembly[patterns[i]])
        	consensus = assembly[patterns[i]]
		# use DLA88*50101 as reference
		d88_50101 = D88_50X.replace('-','')
		length = len(d88_50101)
		exon2N3 = consensus[dic_50X[0]:dic_50X[length-1]+1]
		# replace '-' from 2 tails with 'N'
		exon2N3 = replace_tail_char(exon2N3,'-','N')	
		exon2N3 = exon2N3.replace('-','')
		string = '>' + patterns[i] + '\n' + exon2N3 + '\n'
		out.write(string)
	out.close()

