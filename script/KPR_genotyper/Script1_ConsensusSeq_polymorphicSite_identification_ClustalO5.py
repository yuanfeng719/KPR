#!/usr/bin/python

import sys
import re

# sequence name from clustalO is truncated, and only 60 bp in each sequence row
def reformat(lst):
	dic = {}
	seq_lst = []
	for i in range(len(lst)):
		header = lst[i].split('\n')[0]
		seq = ''.join(lst[i].split('\n')[1:])
		dic[header] = seq
		seq_lst.append(seq)
	return seq_lst,dic

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

# check more-than-once occurance
def percentMore_occur(lst,percentage):
        remove_values_from_list(lst, '-')
        remove_values_from_list(lst, 'N')
        unique = list(set(lst))
        length = len(lst)
        out = []
        for nuc in unique:
                count = 0
                for nuc2 in lst:
                        if nuc == nuc2:
                                count += 1
                if float(count)/length > 0:
                        out.append(nuc)
        if len(out) >= 2:
                flag = True
        else:
                flag = False
        return flag,out

#trim sequence to remove the first and last few nuc from polymorphic sites calling
def trim_aligned_seq(lst,cutoff):
        result = []
        for j in range(len(lst)):
                string = lst[j]
                dic = {}
                real = 0
                for i in range(len(string)):
                        if string[i] != '-':
                                dic[real] = i
                                real += 1
                range_lst = range(min(dic.keys())+cutoff,max(dic.keys())-cutoff+1)
                new_lst = []
                for i in range_lst:
                        new_lst.append(dic[i])
                new = ''
                for i in range(len(string)):
                        if i not in new_lst:
                                new += '-'
                        else:
                                new += string[i]
                result.append(new)
        return (result)

# generate consensus sequence and identify polymorphic sites
def consensus_polymorphic(seq_lst):
        # sequence length
        length = len(seq_lst[0])
        consensus = ''
        polymorphic_site = {}
        new_lst = trim_aligned_seq(seq_lst,5)
        for i in range(length):
                lst = []
                lst2 = [] # without removal to make consensus
                for j in range(len(new_lst)):
                        lst.append(new_lst[j][i])
                        lst2.append(seq_lst[j][i])
                # remove all '-'
                remove_values_from_list(lst,'-')
                remove_values_from_list(lst2,'-')
                # element with highest occurance
                consensus += find_majority(lst2)[0]
                # check if polymorphic #cutoff >= 5%
                flag,poli = percentMore_occur(lst,0.05)
                if flag:
                        polymorphic_site[i] = poli
        return consensus,polymorphic_site

# linkage -> site dictionary
def create_site_dic(lst):
	dic = {}
	for i in range(len(lst)):
		location = re.findall(r'\d+',lst[i])[0]
		nuc = re.findall(r'\D+',lst[i])[0]
		dic[location] = nuc
	return dic

def intersection(lst1, lst2): 
	lst3 = [value for value in lst1 if value in lst2] 
	return lst3 

# merge linkages to groups
def merge_linkages(contig_sites,percentage):
	lst = contig_sites.values()
	length = -1
	while len(lst) != length and len(lst) > 1:
		length = len(lst)
		for i in range(len(lst)-1):
			for j in range(i+1,len(lst)):
				linkage1 = lst[i]
				linkage2 = lst[j]
				if linkage1 != '' and linkage2 != '':
					dic1 = create_site_dic(linkage1)
					dic2 = create_site_dic(linkage2)
					site1 = dic1.keys()
					site2 = dic2.keys()
					intersect = intersection(site1,site2)
					if len(intersect) > 0:
						flag = True
						for k in range(len(intersect)):
							if dic1[intersect[k]] != dic2[intersect[k]]:
								flag = False
						if flag and (float(len(intersect))/len(site1) >= percentage or float(len(intersect))/len(site2) >= percentage):
							merged_linkage = list(set().union(linkage1,linkage2))
							linkage1 = merged_linkage
							lst[i] = merged_linkage
							linkage2 = ''
							lst[j] = ''
		remove_values_from_list(lst,'')
		remove_values_from_list(lst,[])
	lst.sort()
	return lst
	

# polymorphic site linkage
def linkage(dic,polymorphic_site,contig_linkages,percentage):
	contigs = dic.keys()
	sites = polymorphic_site.keys()
	sites.sort()
	contig_sites = {}
	for key in contigs:
		seq = dic[key]
		contig_sites[key] = []
		string = ''
		for i in range(len(sites)):
			nuc = seq[int(sites[i])]
			if nuc != '-':
				substring = str(sites[i]) + nuc + '\t'
				string += substring
				contig_sites[key].append(str(sites[i]) + nuc)
			else:
				string += '\t'
		if 'A' in string or 'T' in string or 'C' in string or 'G' in string:
			print string
			contig_linkages.write('>'+key+'\n')
			contig_linkages.write(string+'\n')
	linkageGroups = merge_linkages(contig_sites,percentage)
	linkage_header_dic = {}
	linkage_seq_dic = {}
	ungrouped_contig = []
	for rec in contigs:
		ungrouped_contig.append(rec)
	for i in range(len(linkageGroups)):
		linkage = linkageGroups[i]
		linkage.sort()
		linkage_string = '\t'.join(linkage)
		linkage_header_dic[linkage_string] = []
		linkage_seq_dic[linkage_string] = []
		for j in range(len(contigs)):
			contig = contigs[j]
			sites = contig_sites[contig]
			if set(sites).issubset(set(linkage)) and sites != []:
				linkage_header_dic[linkage_string].append(contig)
				linkage_seq_dic[linkage_string].append(dic[contig])
				if contig in ungrouped_contig:
					ungrouped_contig.remove(contig)
	return (linkage_header_dic,linkage_seq_dic,ungrouped_contig)
	

#####################################################################3 
if __name__ == '__main__':
	file0 = sys.argv[1]
	percentage = sys.argv[2]
	with open(file0,'r') as f:
        	file=f.read()
	lst = file.split('>')[1:]
	percentage = float(percentage)
	# output files
	cons_seq = open(file0 + '_ConsensusSequence','w')
	poli_sites = open(file0 + '_PolymorphicSites','w')
	contig_linkages = open(file0 + '_ContigLinkages','w')
	linkage_groups = open(file0 + '_ContigGroupsAndOutliers','w')
	# data processing
	seq_lst,dic = reformat(lst)
	consensus,polymorphic_site = consensus_polymorphic(seq_lst)
	print seq_lst
	# print and outputs
	print ('>ConsensusSequence')
	print consensus
	cons_seq.write('>ConsensusSequence\n')
	cons_seq.write(consensus+'\n')
	cons_seq.close()
	print '\n\n'
	keys = polymorphic_site.keys()
	keys.sort()
	for key in keys:
		string = str(key) + '\t' + '\t'.join(polymorphic_site[key])
		print string
		poli_sites.write(string + '\n')
	poli_sites.close()
	print '\n\n'
	linkage_header_dic,linkage_seq_dic,ungrouped_contig = linkage(dic,polymorphic_site,contig_linkages,percentage)
	contig_linkages.close()
	linkages = linkage_header_dic.keys()
	for i in range(len(linkages)):
		linkage = linkages[i]
		string = '>' + linkage + '\n' + ';'.join(linkage_header_dic[linkage]) + '\n'
		linkage_groups.write(string)
	string = '>Ungrouped' + '\n' + ';'.join(ungrouped_contig) + '\n'
	linkage_groups.write(string)
	linkage_groups.close()

