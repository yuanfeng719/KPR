#!/usr/bin/python

import sys
import re

##########################################################
def Index_pairing(query,ref):
        dic = {}
        idxQ = 0
        idxR = 0
        for i in range(len(query)):
                if query[i] != '-':
                        dic[idxQ] = idxR
                        idxQ += 1
                if ref[i] != '-':
                        idxR += 1
        return (dic)

def Reformat_pairing(lst):
        for i in range(len(lst)):
                header = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                if header == 'ConsensusSequence':
                        query_seq = seq
                elif 'DLA-88*5' in header:
                        seq_50X = seq
                elif 'DLA-88*0' in header:
                        seq_norm = seq
                else:
                        seq_D12 = seq
        # pairing between query and reference sequences
        dic_norm = Index_pairing(query_seq,seq_norm)
        dic_50X = Index_pairing(query_seq,seq_50X)
        dic_D12 = Index_pairing(query_seq,seq_D12)
        return (dic_norm,dic_50X,dic_D12,query_seq,seq_50X)

# analyze if seq_50X contains gap
# remove '-' from two tails of a string(sequence)
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

def check_seq_50X(seq_50X,query_seq):
	seq,head,tail = remove_tail_char(seq_50X,'-')
	seq2,head2,tail2 = remove_tail_char(query_seq,'-')
	flag = False
	all_idx_lst_inRef_orig = []
	if '-' in seq:
		flag = True
		all_idx_lst = [m.start() for m in re.finditer('-', seq)]
		all_idx_lst_inRef = [x+head for x in all_idx_lst]
		all_idx_lst_inRef_orig = [x-head2 for x in all_idx_lst_inRef]
		for idx in all_idx_lst_inRef:
			query_seq = query_seq[:idx] + '#' + query_seq[idx+1:]
		query_seq = query_seq.replace('#','')
		query_seq,h,t = remove_tail_char(query_seq,'-')
	return (query_seq,all_idx_lst_inRef_orig,flag)

# search the gaps near the identified places and delete them in the original alignment
def search_closest_gap(seq,lst):
	flag2 = True
	for i in lst:
		flag = True
		j = 0
		length = len(seq)
		while flag and (i+j < length or i-j >= 0):
			if seq[i+j] == '-' and i+j < length:
				seq = seq[:i+j] + '#' + seq[i+j+1:]
				flag = False
			elif seq[i-j] == '-' and i-j >= 0:
                                seq = seq[:i-j] + '#' + seq[i-j+1:]
                                flag = False
			else:
				j += 1
		if flag:
			flag2 = False
	seq = seq.replace('#','')
	return seq,flag2
	
def modify_clustalO_alignment(clustalO_lst,all_idx_lst_inRef_orig,flag):
        if flag:
		new_clustalO_lst = []
		for i in range(len(clustalO_lst)):
			name = clustalO_lst[i].split('\n')[0]
			seq = ''.join(clustalO_lst[i].split('\n')[1:])
			new_seq,flag2 = search_closest_gap(seq,all_idx_lst_inRef_orig)
			if flag2:
				string = name + '\n' + new_seq + '\n'
				new_clustalO_lst.append(string)
		return (new_clustalO_lst)
	else:
		return (clustalO_lst)

def adjust_contig_gaps(ref_lst,new_clustalO_lst):
        refDic = {}
        for i in range(len(ref_lst)):
                allele = ref_lst[i].split('\n')[0]
                seq = ''.join(ref_lst[i].split('\n')[1:])
                refDic[allele] = seq
        norm_allele = refDic['DLA-88*001:01']
        consensus = refDic['ConsensusSequence']
        # counter number of '-' at the begining of norm allele
        counter = 0
        flag = True
        while counter < len(norm_allele) and flag:
                if norm_allele[counter] == '-':
                        counter += 1
                else:
                        flag = False
        # counter number of '-' at the begining of consensus allele
        counter2 = 0
        flag = True
        while counter2 < len(consensus) and flag:
                if consensus[counter2] == '-':
                        counter2 += 1
                else:
                        flag = False
        # remove '-' in the middle and locate the 3 nuc insertion
        norm_allele_simp = norm_allele.replace('-','')
        location = norm_allele_simp.index('GAGCAC') + 6 + counter - counter2
        # adjust contig alignment gaps
        result = []
        for i in range(len(new_clustalO_lst)):
                contig = new_clustalO_lst[i].split('\n')[0]
                seq = ''.join(new_clustalO_lst[i].split('\n')[1:])
                seq_trim,begin,end = remove_tail_char(seq,'-')
                if '---' in seq_trim:
                        idx = seq_trim.index('---')
                        if location != idx + begin:
                                seq_trim = seq_trim.split('---')[0] + '---'.join(seq_trim.split('---')[1:])
                                seq_trim = seq_trim[:location-begin] + '---' + seq_trim[location-begin:]
                string = contig + '\n' + '-'*begin + seq_trim + '-'*end + '\n'
                result.append(string)
                print string
        return result

#######################################################
if __name__ == '__main__':
	ref_alignment = sys.argv[1]
	consensus = sys.argv[2]
	clustalO = sys.argv[3]
	with open(ref_alignment,'r') as f:
        	ref_alignment = f.read()
	with open(consensus,'r') as f:
        	consensus = f.read()
	with open(clustalO,'r') as f:
        	clustalOF = f.read()
	ref_lst = ref_alignment.split('>')[1:]
	clustalO_lst = clustalOF.split('>')[1:]
	# data processing
	dic_norm,dic_50X,dic_D12,query_seq,seq_50X = Reformat_pairing(ref_lst)
	new_query_seq,all_idx_lst_inRef_orig,flag = check_seq_50X(seq_50X,query_seq)
	new_clustalO_lst = modify_clustalO_alignment(clustalO_lst,all_idx_lst_inRef_orig,flag)
        new_clustalO_lst_gapAdjust = adjust_contig_gaps(ref_lst,new_clustalO_lst)
	# output files
	out = open(clustalO + '_corrected','w')
	for i in range(len(new_clustalO_lst_gapAdjust)):
		string = '>' + new_clustalO_lst_gapAdjust[i]
		out.write(string)
	out.close()


