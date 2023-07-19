#!/usr/bin/python

import sys
import re
import os

###########################################
def build_read_pairwise_linkage_dic(read_pairwise_linkage):
	dic = {}
	for i in range(len(read_pairwise_linkage)):
		key = read_pairwise_linkage[i].split('\t')[0]
		value = int(read_pairwise_linkage[i].split('\t')[1])
		dic[key] = value
	return dic

def remove_values_from_list(the_list, val):
	return [value for value in the_list if value != val]

def linkage_dic(contig_linkage):
	dic = {}
	for i in range(len(contig_linkage)):
		contig = contig_linkage[i].split('\n')[0]
		linkage_lst = contig_linkage[i].split('\n')[1].split('\t') # sorted
		linkage_lst = remove_values_from_list(linkage_lst,'')
		dic[contig] = linkage_lst
	return dic

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def find_most_diverse_pair(linkage_lst,pairwise_linkage_dic,pos_linkage_dic,exclude_pair_lst):
        # this function returns index pairs with highest coverage # this pair should have reads spaning it
        result = []
        diversity = 0
        count = 0
	idx = []
        if len(linkage_lst) >= 2:
                length = len(linkage_lst[0])
                for i in range(length-1):
                        for j in range(i+1,length):
                                site_lst = []
				pair = str(re.findall(r'\d+',linkage_lst[0][i])[0]) + '_' + str(re.findall(r'\d+',linkage_lst[0][j])[0])
				if pair not in exclude_pair_lst:
					flag = True
                             		for k in range(len(linkage_lst)):
                                        	rec = linkage_lst[k]
						pair2 = linkage_lst[k][i] + '_' + linkage_lst[k][j]
						if pair2 not in pairwise_linkage_dic.keys():
							flag = False
                                        	tmp = [rec[i],rec[j]]
                                        	site_lst.append(tmp)
                                	unique_site_lst = [list(x) for x in set(tuple(x) for x in site_lst)]
                                	if len(unique_site_lst) > diversity and flag:
                                        	pair = str(re.findall(r'\d+',linkage_lst[0][i])[0]) + '_' + str(re.findall(r'\d+',linkage_lst[0][j])[0])
                                        	count = pos_linkage_dic[pair]
                                       		result = [re.findall(r'\d+',linkage_lst[0][i])[0],re.findall(r'\d+',linkage_lst[0][j])[0]]
						idx = [i,j]
                                        	diversity = len(unique_site_lst)
                                	elif len(unique_site_lst) == diversity and flag:
                                        	pair = str(re.findall(r'\d+',linkage_lst[0][i])[0]) + '_' + str(re.findall(r'\d+',linkage_lst[0][j])[0])
                                        	if pos_linkage_dic[pair] > count:
                                                	count = pos_linkage_dic[pair]
                                                	result = [re.findall(r'\d+',linkage_lst[0][i])[0],re.findall(r'\d+',linkage_lst[0][j])[0]]
							idx = [i,j]
        return result,idx

def assign_frequency(linkage_lst,pairwise_linkage_dic,pos_linkage_dic,total_percentage,exclude_pair_lst,result_dic):
        most_diverse_pair,idx = find_most_diverse_pair(linkage_lst,pairwise_linkage_dic,pos_linkage_dic,exclude_pair_lst)
        pos1 = int(most_diverse_pair[0])
        pos2 = int(most_diverse_pair[1])
	idx1 = int(idx[0])
	idx2 = int(idx[1])
        freq_count = {}
	linkage_dic = {}
	percentage_count = {}
        if len(linkage_lst) >= 2:
                for i in range(len(linkage_lst)):
                        linkage = linkage_lst[i]
                        pair = linkage[idx1] + '_' + linkage[idx2]
			pos_pair = str(re.findall(r'\d+',linkage[idx1])[0]) + '_' + str(re.findall(r'\d+',linkage[idx2])[0])
			exclude_pair_lst = [pos_pair]
			if pair in pairwise_linkage_dic.keys():
				freq = pairwise_linkage_dic[pair]
			else:
				freq = 0.0
			if pair not in freq_count.keys():
				freq_count[pair] = freq
				linkage_dic[pair] = [linkage]
			else:
				linkage_dic[pair].append(linkage)
		# convert freq to percentage
		pairs = linkage_dic.keys()
		total = sum(freq_count.values())
		for i in range(len(pairs)):
			pair = pairs[i]
			percentage = float(freq_count[pair]) * total_percentage / total
			percentage_count[pair] = percentage
		# analyze unique and dup pairs
		pairs = linkage_dic.keys()
		for i in range(len(pairs)):
			pair = pairs[i]
			if len(linkage_dic[pair]) == 1:
				total_percentage -= percentage_count[pair]
		for i in range(len(pairs)):
                        pair = pairs[i]
                        if len(linkage_dic[pair]) == 1:
                                link = '_'.join(linkage_dic[pair][0])
                                result_dic[link] = percentage_count[pair]
			else:
				new_linkage_lst = linkage_dic[pair]
				total_percentage = percentage_count[pair]
				result_dic = assign_frequency(new_linkage_lst,pairwise_linkage_dic,pos_linkage_dic,total_percentage,exclude_pair_lst,result_dic)
	return (result_dic)

def nucFreq_dic(nuc_freq):
        dic = {}
        header = nuc_freq.split('\n')[0].split('\t')[1:]
        lst = nuc_freq.split('\n')[1:-1]
        for i in range(len(lst)):
                pos,c1,c2,c3,c4 = lst[i].split('\t')
                dic[int(pos)] = {}
                dic[int(pos)][header[0]] = int(c1)
                dic[int(pos)][header[1]] = int(c2)
                dic[int(pos)][header[2]] = int(c3)
                dic[int(pos)][header[3]] = int(c4)
        return dic

def identify_div_site(nuc_freq_dic,excludeDic,remainContig):
        siteLst = nuc_freq_dic.keys()
        siteMax = -1
        depthMax = -1
	excludeLst = excludeDic[remainContig[0]]
        for site in siteLst:
                if site not in excludeLst:
                        depth = sum(nuc_freq_dic[site].values())
                        if depth > depthMax:
                                depthMax = depth
                                siteMax = site
        return siteMax,depthMax

def perc_assignment(totalPerc,pos0,contig_linkage_dic,contigLst,nuc_freq_dic):
        dic = {}
        for contig in contigLst:
                linkage = contig_linkage_dic[contig]
                for site in linkage:
                        pos = int(re.findall(r'\d+',site)[0])
                        nuc = re.findall(r'\D+',site)[0]
                        if pos == pos0:
                                if nuc not in dic.keys():
                                        dic[nuc] = []
                                dic[nuc].append(contig)
        percDic = {}
        totalDepth = 0
        moveOn = False
        nucLst = dic.keys()
        for nuc in nucLst:
                contigLst = dic[nuc]
                totalDepth += nuc_freq_dic[pos0][nuc]
                if len(contigLst) > 1:
                        moveOn = True
        for nuc in nucLst:
                perc = totalPerc * float(nuc_freq_dic[pos0][nuc]) / totalDepth
                percDic[nuc] = perc
        return dic,percDic,moveOn

def dic_all_key_append(dic,char):
	lst = dic.keys()
	for item in lst:
		if char not in dic[item]:
			dic[item].append(char)
	return dic

def contig_freq_assignment_singleSite(contig_linkage_dic,nuc_freq_dic,totalPerc,excludeDic,moveOn,contigLst,resultDic,remainContig):
	siteLst = [] # make sure sites in nuc_freq_dic match contig linkage sites
	for site in contig_linkage_dic[contigLst[0]]:
		pos = int(re.findall(r'\d+',site)[0])
		siteLst.append(pos)
        while moveOn:
		siteMax,depthMax = identify_div_site(nuc_freq_dic,excludeDic,remainContig)
		while (siteMax not in siteLst) and (siteMax != -1):
			excludeDic = dic_all_key_append(excludeDic,siteMax) # remove disagreed site
                	siteMax,depthMax = identify_div_site(nuc_freq_dic,excludeDic,remainContig) # repeat most div site identification
                if siteMax == -1:
			moveOn = False
		else:
	                sepDic,percDic,moveOn = perc_assignment(totalPerc,siteMax,contig_linkage_dic,contigLst,nuc_freq_dic)
			for contig in remainContig: # site used for these contigs
				excludeDic[contig].append(siteMax)
                	if len(sepDic.keys()) > 1:
                        	nucLst = sepDic.keys()
                        	for nuc in nucLst:
                                	contigLst = sepDic[nuc]
                                	totalPerc = percDic[nuc]
					remainContig = []
                                	for contig in contigLst:
						remainContig.append(contig)
                                        	if contig not in resultDic.keys():
                                                	resultDic[contig] = {}
                                        	resultDic[contig][siteMax] = totalPerc
                                	resultDic = contig_freq_assignment_singleSite(contig_linkage_dic,nuc_freq_dic,totalPerc,excludeDic,moveOn,contigLst,resultDic,remainContig)
					if len(resultDic[contig]) >= len(siteLst):
						moveOn = False
	return resultDic

def single_contig_freq_assignment(resultDic,totalPerc,contigLst):
        siteMax = 'NA'
        for contig in contigLst:
                resultDic[contig] = {}
                resultDic[contig][siteMax] = totalPerc
        return resultDic

def sort_sites(resultDic,nuc_freq_dic,contig_linkage_dic):
	# reformat contig_linkage_dic
        contig_linkage_dic_tmp = {}
        contigLst = contig_linkage_dic.keys()
        for contig in contigLst:
                if contig not in contig_linkage_dic_tmp.keys():
                        contig_linkage_dic_tmp[contig] = {}
                siteLst = contig_linkage_dic[contig]
                for site in siteLst:
                        pos = int(re.findall(r'\d+',site)[0])
                        nuc = re.findall(r'\D+',site)[0]
                        contig_linkage_dic_tmp[contig][pos] = nuc
        dic = {}
	percDic = {}
	contigLst = resultDic.keys()
	tmp = {}
	for contig in contigLst:
		if contig not in tmp.keys():
			tmp[contig] = {}
		siteLst = resultDic[contig].keys()
		percLst = []
		for site in siteLst:
                        if site in resultDic[contig].keys():
			        perc = resultDic[contig][site] 
			        tmp[contig][perc] = site
			        percLst.append(perc)
		percLst.sort(reverse = True)
		dic[contig] = []
		for perc in percLst:
			site = tmp[contig][perc]
                        if site in contig_linkage_dic_tmp[contig].keys():
                                nuc = contig_linkage_dic_tmp[contig][site]
                                count = float(sum(nuc_freq_dic[site].values())) * perc
			        string = str(site) + '(' + str(perc) + ';' + str(count) + ')'
			        dic[contig].append(string)
		perc = percLst[-1]
		percDic[contig] = perc
	return dic,percDic

def count_readFile(file0):
        if os.path.isfile(file0):
                with open(file0,'r') as f:
                        file = f.read()
                if '>' in file:
                        lst = file.split('>')[1:]
                        result = len(lst)
                else:
                        result = 0
        else:
                result = 0
        return result

def prot_classification(prot_seq):
        normContig = []
        X50Contig = []
        if os.path.isfile(prot_seq):
                with open(prot_seq) as f:
                        file = f.read()
        else:
                file = []
        if '>' in file:
                lst = file.split('>')[1:]
                for i in range(len(lst)):
                        name = lst[i].split('\n')[0]
                        seq = lst[i].split('\n')[1]
                        if len(seq) > 182:
                                X50Contig.append(name)
                        else:
                                normContig.append(name)
        return normContig,X50Contig

#########################################
if __name__ == '__main__':
	contig_linkage = sys.argv[1]
        nuc_freq = sys.argv[2]
        prot_seq = sys.argv[3]
        norm_read = sys.argv[4]
        X50_read = sys.argv[5]
	with open(contig_linkage,'r') as f:
        	contig_linkage = f.read()
        with open(nuc_freq,'r') as f:
                nuc_freq = f.read()
	contig_linkage = contig_linkage.split('>')[1:]
        norm = count_readFile(norm_read)
        X50 = count_readFile(X50_read)
        normContig,X50Contig = prot_classification(prot_seq)
	# output files
	out = open('Correct_contigs_percentageScores.txt','w')
	out2 = open('Correct_contigs_percentage_siteChain.txt','w')
	# data processing
	contig_linkage_dic = linkage_dic(contig_linkage)
        nuc_freq_dic = nucFreq_dic(nuc_freq)

        totalPerc = 1.0
        if norm + X50 != 0:
                normPerc = totalPerc * norm / (norm + X50)
                X50Perc = totalPerc * X50 / (norm + X50)
        else:
                normPerc = 0.0
                X50Perc = 0.0
	resultDic = {}
	
        # norm
        if normContig != []:
                excludeDic = {} # indicate the excluded sites for each contig
                contigLst = normContig
                moveOn = True
	        remainContig = normContig
	        for contig in contigLst:
		        excludeDic[contig] = []
                
                if len(normContig) > 1:
	                resultDic = contig_freq_assignment_singleSite(contig_linkage_dic,nuc_freq_dic,normPerc,excludeDic,moveOn,contigLst,resultDic,remainContig)
                else:
                        resultDic = single_contig_freq_assignment(resultDic,normPerc,contigLst)
        # X50
        if X50Contig != []:
                excludeDic = {} # indicate the excluded sites for each contig
                contigLst = X50Contig
                moveOn = True
                remainContig = X50Contig
                for contig in contigLst:
                        excludeDic[contig] = []
                
                if len(X50Contig) > 1:
                        resultDic = contig_freq_assignment_singleSite(contig_linkage_dic,nuc_freq_dic,X50Perc,excludeDic,moveOn,contigLst,resultDic,remainContig)
                else:
                        resultDic = single_contig_freq_assignment(resultDic,X50Perc,contigLst)
       
	siteChainDic,percDic = sort_sites(resultDic,nuc_freq_dic,contig_linkage_dic)

	contigLst = siteChainDic.keys()
	contigLst.sort(key=natural_keys)
	for contig in contigLst:
		perc = percDic[contig]
		string = '>' + contig + '\n' + str(perc) + '\n' + '\t'.join(contig_linkage_dic[contig]) + '\n'
		out.write(string)
		siteChain = siteChainDic[contig]
		string = contig + '\t' + ','.join(siteChain) + '\n'
		out2.write(string)
	out.close()
	out2.close()

