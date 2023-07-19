#!/usr/bin/python

import sys
import re

#################################################
#sequence name from clustalO is truncated, and only 60 bp in each sequence row
def rename_and_reformat(lst):
        dic = {}
        contig_lst = []
        for i in range(len(lst)):
                header = lst[i].split("\n")[0]
                seq = ''.join(lst[i].split('\n')[1:])
                dic[header] = seq
                contig_lst.append(header)
        return contig_lst,dic

#Old definition changing the Align's number
# def rename_and_reformat(lst):
#         dic = {}
#         contig_lst = []
#         for i in range(len(lst)):
#                 header = lst[i].split("\n")[0]
#                 seq = ''.join(lst[i].split('\n')[1:])
#                 dic[header] = seq
#                 contig_lst.append(header)
#         return contig_lst,dic

# remove all occurance of a character from a list
def remove_values_from_list(the_list, val):
        while val in the_list:
                the_list.remove(val)

# identify the element with highest occurance in a list
def find_majority(k):
        myMap = {}
        maximum = ( '', 0 ) # (occurring element, occurrences)
        for n in k:
                if n in myMap: myMap[n] += 1
                else: myMap[n] = 1
                if myMap[n] > maximum[1]: maximum = (n,myMap[n])
        return maximum

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
        i = 0
        while flag and i < len(string):
                if string[i] == char:
                        i += 1
                else:
                        flag = False
                        string = string[i:]
        string = string[::-1]
	return (string)
		

# Generate consensus sequence based on clustalO alignment given contigs
def consensus_seq_given_contigs(dic,lst):
	seq_lst = []
	#print(dic.keys())
	#print(lst)
	for i in range(len(lst)):
		#print(dic.keys())
		seq = dic[lst[i]]
		seq_lst.append(seq)
	# sequence length
        length = len(seq_lst[0])
        consensus = ''
        for i in range(length):
                lst = []
                for j in range(len(seq_lst)):
                        lst.append(seq_lst[j][i])
                # remove all '-'
                remove_values_from_list(lst,'-')
		if lst == []:
			lst = ['N']
                # element with highest occurance
                consensus += find_majority(lst)[0]
	#consensus = remove_tail_char(consensus,'N')
	return consensus

# generate consensus sequences for groups
def group_consensus(dic,group_lst):
	out = {}
	used_contigs = []
	ungrouped_contigs = []
	for i in range(len(group_lst)):
		pattern = group_lst[i].split('\n')[0]
		contigs = group_lst[i].split('\n')[1].split(';')
		if pattern == 'Ungrouped':
			ungrouped_contig = contigs
		else:
			if contigs != ['']:
				for cont in contigs:
					used_contigs.append(cont)
				consensus = consensus_seq_given_contigs(dic,contigs)
				out[pattern] = consensus
	used_contigs = list(set(used_contigs))
	return (out,used_contigs,ungrouped_contigs)

# sequence dictionary
def build_seq_dic(lst):
        dic = {}
        for i in range(len(lst)):
                header = lst[i].split('\n')[0]
                seq = lst[i].split('\n')[1]
                dic[header] = seq
        return dic

# check dictionary
def check_dic(dic,lst):
	out = []
	for i in range(len(lst)):
		result = dic[int(lst[i])]
		out.append(result)
	return out

#
def consensus_including_unused(unused_contigs,dic,consensus_dic,forwd_dic,rev_dic,pattern_contig_dic,shared,dic_norm,dic_50X,dic_D12):
	patterns = pattern_contig_dic.keys()
	if 'Ungrouped' in patterns:
		patterns.remove('Ungrouped')
	result = {}
	for i in range(len(patterns)):
		sites = []
		ref_sites_norm = []
		ref_sites_50X = []
		ref_sites_D12 = []
		consensus = consensus_dic[patterns[i]]
		sites = re.findall(r'\d+',patterns[i])
		ref_sites_norm = check_dic(dic_norm,sites)
		ref_sites_50X = check_dic(dic_50X,sites)
		ref_sites_D12 = check_dic(dic_D12,sites)
		unused = []
		for j in range(len(unused_contigs)):
			unused_seq = dic[unused_contigs[j]]
			flag = False
			for k in range(len(shared)):
				forwd = forwd_dic[shared[k]]
				forwd_allele = forwd[0]
				forwd_pos = forwd[1]
				forwd_seq = forwd[2]
				forwd_ref_sites = []
				if 'DLA-88*5' in forwd_allele:
					forwd_ref_sites = ref_sites_50X
				elif 'DLA-12' in forwd_allele:
					forwd_ref_sites = ref_sites_D12
				else:
					forwd_ref_sites = ref_sites_norm
				rev = rev_dic[shared[k]]
				rev_allele = rev[0]
				rev_pos = rev[1]
				rev_seq = rev[2]
				rev_ref_sites = []
				if 'DLA-88*5' in rev_allele:
					rev_ref_sites = ref_sites_50X
				elif 'DLA-12' in rev_allele:
					rev_ref_sites = ref_sites_D12
				else:
					rev_ref_sites = ref_sites_norm
				flag2 = False
				for l in range(len(forwd_ref_sites)):
					if int(forwd_pos)-1 <= int(forwd_ref_sites[l]) and int(forwd_pos)-1+len(forwd_seq) >= int(forwd_ref_sites[l]):
						flag2 = True
				for l in range(len(rev_ref_sites)):
					if int(rev_pos)-1 <= int(rev_ref_sites[l]) and int(rev_pos)-1+len(rev_seq) >= int(rev_ref_sites[l]):
						flag2 = True
				if flag2 and (forwd_seq in consensus and rev_seq in unused_seq) or (rev_seq in consensus and forwd_seq in unused_seq):
					flag = True
			if flag:
				unused.append(unused_contigs[j])
		# combine pattern contig and unused contigs
		New_lst = []
		pats = pattern_contig_dic[patterns[i]]
		for p in pats:
			New_lst.append(p)
		for p in unused:
			New_lst.append(p)
		result[patterns[i]] = New_lst
	return (result)

# Build sequence dic for forwd and reverse reads
def sequence_allele_pos_dic(lst):
	dic = {}
	for i in range(len(lst)):
		name = lst[i].split('\t')[0]
		allele = lst[i].split('\t')[2]
		pos = lst[i].split('\t')[3]
		seq = lst[i].split('\t')[9]
		if name not in dic.keys():
			dic[name] = [allele,pos,seq]
	return dic

# Index pairing between consensus and references based on clustalO result
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
        return (dic_norm,dic_50X,dic_D12)

def filter_contigs_using_polySite(pattern,contigs,dic):
        result = []
        polySiteLst = pattern.split('\t')
        for contig in contigs:
                seq = dic[contig]
                flag = True
                for polySite in polySiteLst:
                        if flag:
                                pos = int(re.findall(r'\d+', polySite)[0])
                                nuc = re.findall(r'\D+', polySite)[0]
                                if (seq[pos] != nuc) and (seq[pos] != '-'):
                                        flag = False
                if flag:
                        result.append(contig)
        return result

##################################################
if __name__ == '__main__':
	alignment = sys.argv[1]
	contig_groups = sys.argv[2]
	forwd = sys.argv[3]
	rev = sys.argv[4]
	db_alignment = sys.argv[5]
	with open(alignment,'r') as f:
        	alignment = f.read()
	with open(contig_groups,'r') as f:
        	groups = f.read()
	with open(forwd,'r') as f:
        	forwd = f.read()
	with open(rev,'r') as f:
        	rev = f.read()
	with open(db_alignment,'r') as f:
        	db_alignment = f.read()
	lst = alignment.split('>')[1:]
	group_lst = groups.split('>')[1:]
	forwd = forwd.split('\n')[:-1]
	rev = rev.split('\n')[:-1]
	db_alignment = db_alignment.split('>')[1:]
	# output files
	out = open(contig_groups + '_PatternConsensusSequences','w')
	out2 = open(contig_groups + '_UnusedContigs','w')
	out3 = open(contig_groups + '_ConsensusSequencesWithUnusedContigs','w')
	# data processing
	contig_lst,dic = rename_and_reformat(lst)
	consensus_dic,used_contigs,ungrouped_contigs = group_consensus(dic,group_lst)
	# consensus_dic contains consensus sequences from all linkage pattern groups
	keys = consensus_dic.keys()
	for key in keys:
		string = '>' + key + '\n' + consensus_dic[key] + '\n'
		out.write(string)
	out.close()
	# remove all used contigs from total and extract unused contigs
	for i in range(len(used_contigs)):
		remove_values_from_list(contig_lst, used_contigs[i])
	for i in range(len(contig_lst)):
		seq = dic[contig_lst[i]]
		# remove '-' from string
		seq = re.sub('-', '', seq)
		string = '>' + contig_lst[i] + '\n' + seq + '\n'
		out2.write(string)
	out2.close()
	# contig_lst contains the names of all unused contigs
	# build sequence dictinoary for both forwd and rev reads
	forwd_dic = sequence_allele_pos_dic(forwd)
	rev_dic = sequence_allele_pos_dic(rev)
	shared = list(set(forwd_dic.keys()).intersection(rev_dic.keys())) # paired reads in both forwd and rev files
	dic_norm,dic_50X,dic_D12 = Reformat_pairing(db_alignment)
	# pattern and contigs dictionary
	pattern_contig_dic = {}
	for i in range(len(group_lst)):
		pattern = group_lst[i].split('\n')[0]
		contigs = group_lst[i].split('\n')[1].split(';')
		if contigs != ['']:
			pattern_contig_dic[pattern] = contigs
	pattern_contigLst_with_unused = consensus_including_unused(contig_lst,dic,consensus_dic,forwd_dic,rev_dic,pattern_contig_dic,shared,dic_norm,dic_50X,dic_D12)
	patterns = pattern_contigLst_with_unused.keys()
        for i in range(len(patterns)):
		pattern = patterns[i]
		contigs = pattern_contigLst_with_unused[pattern]
                contigs_checked = filter_contigs_using_polySite(pattern,contigs,dic)
		consensus = consensus_seq_given_contigs(dic,contigs_checked)
                # consensus seq may contain artificial error. Check linkage again
                flag = False
                for contig in contigs_checked:
                        seq = dic[contig][5:-6] # terminals are more likely to contain sequencing errors, remove 5 nucs
                        if seq in consensus:
                                flag = True
                if flag:
		        string = '>' + pattern + '\n' + consensus + '\n'
		        out3.write(string)
	out3.close()

