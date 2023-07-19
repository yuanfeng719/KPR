#!/usr/bin/python

import sys
import re

def linkage_dic(contigLinkage):
        dic = {}
        with open(contigLinkage,'r') as f:
                contigLinkage = f.read()
        lst = contigLinkage.split('>')[1:]
        for i in range(len(lst)):
                contig = lst[i].split('\n')[0]
                linkage = lst[i].split('\n')[1].split('\t')
                linkage2 = [i for i in linkage if i != '']
                dic[contig] = linkage2
        return dic

def natural_key(string_):
        return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_) if s]

def sam_analyzer(Sam,polymorphicSite,startPos): # BWA reference includes E1,4-8 but polySite inclueds E2-3 only
	with open(polymorphicSite,'r') as f:
		polySite = f.read()
	with open(Sam,'r') as f:
		sam = f.read()
	# check if norm and 50X contigs exist simutaniously
	flag = False
	if ('50XAlign' in sam) and ('NormAlign' in sam):
		flag = True
	# poly site dic
	polySite = polySite.split('\n')[:-1]
	polySiteDic = {}
	for i in range(len(polySite)):
		pos = int(polySite[i].split('\t')[0])
		nucLst = polySite[i].split('\t')[1:]
		polySiteDic[pos] = nucLst
	polySiteLst = polySiteDic.keys()
	polySiteLst.sort()
	# read pair poly site matching
	sam = sam.split('\n')[:-1]
	samDic = {}
	for i in range(len(sam)):
		info = sam[i].split('\t')
		contig = info[2]
		read = info[0]
		pos = int(info[3])
		seq = info[9]
		if 'NM:i:0' in sam[i]:
			if read not in samDic.keys():
				samDic[read] = []
			if flag:
				if (pos + len(seq)) >= (startPos + 464) and pos <= (startPos + 464):
					if 'NormAlign' in contig:
						seq = seq[:(startPos + 464) - pos] + '---' + seq[(startPos + 464) - pos:]
				if pos > startPos + 464:
					if 'NormAlign' in contig:
						pos += 3
			for site in polySiteLst:
				if (site >= pos - startPos) and (site < pos + len(seq) - startPos):
					nuc = seq[site - (pos - startPos)]
					string = str(site) + nuc
					samDic[read].append(string)
	# dedup (read pair)
	readLst = samDic.keys()
	for read in readLst:
		siteLst = list(set(samDic[read]))
		siteLst = sorted(siteLst,key=natural_key)
		samDic[read] = siteLst
	return samDic

def summerize_pairwise_linkage(samDic):
	posDic = {}
	linkDic = {}
	freqDic = {}
	readLst = samDic.keys()
        for read in readLst:
                linkLst = samDic[read]
		if len(linkLst) >= 2:
			for i in range(len(linkLst)-1):
				site1 = linkLst[i]
				p1 = int(re.findall(r'\d+',site1)[0])
				for j in range(i+1,len(linkLst)):
					site2 = linkLst[j]
					p2 = int(re.findall(r'\d+',site2)[0])
					if p1 < p2:
						linkPair = site1 + '_' + site2
						posPair = str(p1) + '_' + str(p2)
					else:
						linkPair = site2 + '_' + site1
						posPair = str(p2) + '_' + str(p1)
					if posPair not in posDic.keys():
						posDic[posPair] = 0
					posDic[posPair] += 1
					if linkPair not in linkDic.keys():
						linkDic[linkPair] = 0
					linkDic[linkPair] += 1
		for link in linkLst:
			pos = int(re.findall(r'\d+',link)[0])
			if pos not in freqDic.keys():
				freqDic[pos] = {}
				freqDic[pos]['A'] = 0
				freqDic[pos]['C'] = 0
				freqDic[pos]['G'] = 0
				freqDic[pos]['T'] = 0
			nuc = link.split(str(pos))[1]
			freqDic[pos][nuc] += 1
	return posDic,linkDic,freqDic
					

if __name__ == '__main__':
	polymorphicSite = sys.argv[1]
	contigLinkage = sys.argv[2]
	Sam = sys.argv[3]
	# read in
	linkageDic = linkage_dic(contigLinkage)
	samDic = sam_analyzer(Sam,polymorphicSite,74)
	# result
	out = open('ReadPair_polySite_linkages.txt','w')
	out2 = open('Pairwise_position_counts.txt','w')
	out3 = open('Pairwise_linkage.txt','w')
	out4 = open('PolymorphicSites_nucleotide_freq.txt','w')
	readLst = samDic.keys()
	for read in readLst:
		linkLst = samDic[read]
		string = '>' + read + '\n' + '\t'.join(linkLst) + '\n'
		out.write(string)
	out.close()
	posDic,linkDic,freqDic = summerize_pairwise_linkage(samDic)
	posLst = posDic.keys()
	linkLst = linkDic.keys()
	for pos in posLst:
		count = str(posDic[pos])
		string = pos + '\t' + count + '\n'
		out2.write(string)
	out2.close()
	for link in linkLst:
		count = str(linkDic[link])
		string = link + '\t' + count + '\n'
		out3.write(string)
	out3.close()
	posLst = freqDic.keys()
	header = '\tA\tC\tG\tT\n'
	out4.write(header)
	for pos in posLst:
		string = str(pos) + '\t' + str(freqDic[pos]['A']) + '\t' + str(freqDic[pos]['C'])+ '\t' + str(freqDic[pos]['G'])+ '\t' + str(freqDic[pos]['T']) + '\n'
		out4.write(string)
	out4.close()


