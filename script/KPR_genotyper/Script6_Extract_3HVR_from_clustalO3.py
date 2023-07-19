#!/usr/bin/python

import sys

######################################################
# Index pairing between consensus and references based on clustalO result
def Index_pairing(query,ref):
        dic = {}
        idxQ = 0
        idxR = 0
	flag = True
        for i in range(len(query)):
                if query[i] != '-':
			flag = True
                        dic[idxQ] = idxR
                        idxQ += 1
                        idxR += 1
                elif (query[i] == '-' and ref[i] != 'X') and flag:
                        idxR += 1
		elif query[i] == '-' and ref[i] == 'X':
			flag = False
        return (dic)

def remove_tail_char(string,char):
        flag = True
        i = 0
        while flag and i < len(string):
                if string[i] == char:
                        i += 1
                else:
                        flag = False
                        string = string[i:]
        string = string[::-1]
        flag = True
        j = 0
        while flag and j < len(string):
                if string[j] == char:
                        j += 1
                else:
                        flag = False
                        string = string[j:]
        string = string[::-1]
        return (string,i,j)

def extract_seq(consensus,start,end):
	seq = consensus[start:end]
	seq2,i,j = remove_tail_char(seq,'-')
	while start > 0 and i > 0:
		start -= 1
		if consensus[start] != '-':
			i -= 1
	while end < len(consensus)-1 and j > 0:
		end += 1
		if consensus[end] != '-':
			j -= 1
	seq = consensus[start:end]
	return seq

#######################################################
if __name__ == '__main__':
	file0 = sys.argv[1]
	with open(file0,'r') as f:
        	file=f.read()
	lst = file.split('>')[1:]
	# output files
	out = open(file0 + '_3HVR','w')
	out2 = open(file0 + '_Exon2N3','w')
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
	patterns = assembly.keys()
	for i in range(len(patterns)):
		dic_50X = Index_pairing(D88_50X,assembly[patterns[i]])
		dic_norm = Index_pairing(D88_norm,assembly[patterns[i]])
		consensus = assembly[patterns[i]]
		# HVR1
		HVR1 = ''
		for j in range(61,77):
			if consensus[dic_norm[j]] != '-':
				HVR1 += consensus[dic_norm[j]]
		# HVR2
		HVR2 = ''
		for j in range(90,116):
			if consensus[dic_norm[j]] != '-':
				HVR2 += consensus[dic_norm[j]]
		# HVR3
		#HVR3 = consensus[dic_50X[151]:dic_50X[158]]
		HVR3 = extract_seq(consensus,dic_50X[151],dic_50X[158])
		HVR3 = HVR3.replace('-','')
		string = '>' + patterns[i] + '\n' + HVR1 + HVR2 + HVR3 + '\n'
		out.write(string)
		# Exon 2 and 3
		whole = ''
		first_idx = dic_50X[0]
		last_idx = dic_50X[182]
		#whole = consensus[first_idx:last_idx+1]
		whole = extract_seq(consensus,first_idx,last_idx+1)
		whole = whole.replace('-','')
		string = '>' + patterns[i] + '\n' + whole + '\n'
		out2.write(string)
	out.close()
	out2.close()

