#!/usr/bin/python

import sys
import glob
import re


def pos_dic(posReadCount):
	dic = {}
	with open(posReadCount,'r') as f:
		posReadCount = f.read()
	lst = posReadCount.split('\n')[1:-1]
	header = posReadCount.split('\n')[0].split('\t')[1:]
	for i in range(len(lst)):
		pos,c1,c2,c3,c4 = lst[i].split('\t')
		dic[pos + header[0]] = c1
		dic[pos + header[1]] = c2
		dic[pos + header[2]] = c3
		dic[pos + header[3]] = c4
	return dic

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

def seq_dic(nucSeq):
	dic = {}
	with open(nucSeq,'r') as f:
                nucSeq = f.read()
        lst = nucSeq.split('>')[1:]
        for i in range(len(lst)):
                contig = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:]).replace('-','')
                dic[contig] = seq
        return dic

def read_assignment(readLinkage,linkageDic):
	dic = {}
	with open(readLinkage,'r') as f:
		readLinkage = f.read()
	lst = readLinkage.split('>')[1:]
	contigLst = linkageDic.keys()
	for contig in contigLst:
		if contig not in dic.keys():
			dic[contig] = {}
		linkage = linkageDic[contig]
		for i in range(len(lst)):
			read_lst = lst[i].split('\n')[1].split('\t')
			readKey = ';'.join(read_lst)
			if readKey in dic[contig].keys():
				dic[contig][readKey] += 1
			else:
				flag = True
				for pos in read_lst:
					if pos not in linkage:
						flag = False
				if flag:
					dic[contig][readKey] = 1
	return dic

def natural_key(string_):
	return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_) if s]

def pairwise_pos_dic(file0):
	dic = {}
	with open(file0,'r') as f:
		file = f.read()
	lst = file.split('\n')[:-1]
	for i in range(len(lst)):
		pair,count = lst[i].split('\t')
		if pair != '':
			dic[pair] = int(count)
	return dic

def compare_2seq(seq1,seq2,cutoff):
        flag = True
        counter = 0
        if len(seq1) <= len(seq2):
                length = len(seq1)
        else:
                length = len(seq2)
        while (counter < length) and flag:
                if seq1[counter] != seq2[counter]:
                        flag = False
                counter += 1
        flag2 = True
        counter2 = 0
        seq1_rev = seq1[::-1]
        seq2_rev = seq2[::-1]
        while (counter2 < length) and flag:
                if seq1_rev[counter2] != seq2_rev[counter2]:
                        flag2 = False
                counter2 += 1
        chimeric = False
        if (counter >= cutoff) or (counter2 >= cutoff):
                chimeric = True
        return chimeric

def identical_head_tail(seqDic,cutoff):
        problem_contig = []
        contigLst = seqDic.keys()
        if len(contigLst) >= 2:
                for i in range(len(contigLst)-1):
                        contig1 = contigLst[i]
                        seq1 = seqDic[contig1]
                        for j in range(i+1,len(contigLst)):
                                contig2 = contigLst[j]
                                seq2 = seqDic[contig2]
                                chimeric = compare_2seq(seq1,seq2,cutoff)
                                if chimeric == True:
                                        problem_contig.append(contig1)
                                        problem_contig.append(contig2)
        problem_contig_result = list(set(problem_contig))
        return problem_contig_result

def contig_validation(linkageDic,pairPosDic,casePosDic,readAssignDic,DLA88_depthCutoff,DLA12_depthCutoff,seqDic,lengthCutoff,genoDic):
	dic1 = {} # depth >= cutoff and count == 0
	dic2 = {} # 0 < depth < cutoff and count == 0
	tagDic = {}
	contigLst = linkageDic.keys()
        print contigLst
        print genoDic
	for contig in contigLst:
                genoInfo = genoDic[contig]
                if 'DLA-88' in genoInfo:
                        cutoff = DLA88_depthCutoff
                else:
                        cutoff = DLA12_depthCutoff
		linkage = linkageDic[contig]
		flag1 = True
		flag2 = True
		if len(linkage) >= 2:
			for idx1 in range(len(linkage)-1):
				pos1 = linkage[idx1]
				p1 = int(re.findall(r'\d+',pos1)[0])
				for idx2 in range(idx1+1,len(linkage)):
					pos2 = linkage[idx2]
					p2 = int(re.findall(r'\d+',pos2)[0])
					if p1 < p2:
						pair = str(p1) + '_' + str(p2)
						casePair = pos1 + '_' + pos2
					else:
						pair = str(p2) + '_' + str(p1)
						casePair = pos2 + '_' + pos1
					if pair in pairPosDic.keys():
						depth = pairPosDic[pair]
					else:
						depth = 0
					if casePair in casePosDic.keys():
						count = casePosDic[casePair]
					else:
						count = 0
					if depth >= cutoff:
						if count == 0:
							flag1 = False
							if contig not in dic1.keys():
								dic1[contig] = []
							dic1[contig].append(casePair)
					else:
						if count == 0:
							flag2 = False
							if contig not in dic2.keys():
                                                                dic2[contig] = []
                                                        dic2[contig].append(casePair)
		if contig not in dic1.keys():
			dic1[contig] = []
		if contig not in dic2.keys():
			dic2[contig] = []
		if flag1 and flag2:
			tagDic[contig] = 'Confident'
		elif flag1:
			tagDic[contig] = 'Possibly chimeric (Pairwise linkage validation failed)'
		elif flag2:
			tagDic[contig] = 'Chimeric'
		else:
			tagDic[contig] = 'Chimeric'
		# check if a single read pair covers the entire linkage (High confidence)
		readLinkage = readAssignDic[contig].keys()
		if ';'.join(linkage) in readLinkage:
			tagDic[contig] = 'Highly confident'
        ## special cases
        # [1] onyly one contig (without polymorphic sites)
        if len(contigLst) == 1:
                for contig in contigLst:
                        tagDic[contig] = 'Unable to determine'
        # [2] two contigs shares completely identical begining and tail, likely caused by chimeric assembly (length >= cutoff)
        problem_contig = identical_head_tail(seqDic,lengthCutoff)
        for contig in problem_contig:
                if tagDic[contig] == 'Possibly chimeric (Pairwise linkage validation failed)':
                        tagDic[contig] = 'Possibly chimeric (Pairwise linkage validation failed & Long identical head/tail)'
                elif tagDic[contig] != 'Chimeric':
                        tagDic[contig] = 'Possibly chimeric (Long identical head/tail)'
        # known alleles
        contigLst = genoDic.keys()
        for contig in contigLst:
                if contig in tagDic.keys():
                        info = genoDic[contig].split('\t')
                        allele = info[3]
                        if allele != '':
                                tagDic[contig] = 'Known allele'
	return tagDic,dic1,dic2

def evidence_dic(linkageDic,pairPosDic,casePosDic):
	dic = {}
	contigLst = linkageDic.keys()
	for contig in contigLst:
		linkage = linkageDic[contig]
		if contig not in dic.keys():
			dic[contig] = []
		if len(linkage) >= 2:
			for i in range(len(linkage)-1):
				site1 = linkage[i]
				p1 = int(re.findall(r'\d+',site1)[0])
				for j in range(i+1,len(linkage)):
					site2 = linkage[j]
					p2= int(re.findall(r'\d+',site2)[0])
					if p1 < p2:
                                                posPair = str(p1) + '_' + str(p2)
                                                linkPair = site1 + '_' + site2
                                        else:
                                                posPair = str(p2) + '_' + str(p1)
                                                linkPair = site2 + '_' + site1
					if posPair in pairPosDic.keys():
						posCount = pairPosDic[posPair]
					else:
						posCount = 0
					if linkPair in casePosDic.keys():
						linkCount = casePosDic[linkPair]
					else:
						linkCount = 0
					string = linkPair + '(' + str(linkCount) + '/' + str(posCount) + ')'
					dic[contig].append(string)
	return dic

def adjust_relativeAF_basedOn_tag(genoDic,tagDic):
        X50_total = 0.0
        X50_failed = 0.0
        norm_total = 0.0
        norm_failed = 0.0
        result = {}
        origPerc = {} # still record the RAF before adjustment
        contigLst = genoDic.keys()
        for contig in contigLst:
                info = genoDic[contig].split('\t')
                contig = info[1]
                perc = info[2]
                origPerc[contig] = perc
                if 'DLA64' in genoDic[contig]:
                        result[contig] = genoDic[contig]
                elif '50X' in genoDic[contig]:
                        contig = info[1]
                        perc = float(info[2])
                        X50_total += perc
                        if contig in tagDic.keys():
                                if tagDic[contig] == 'Chimeric':
                                        X50_failed += perc
                                        tmp = '\t'.join(info[:2]) + '\t' + str(0.0) + '\t' + '\t'.join(info[3:])
                                        result[contig] = tmp
                else:
                        contig = info[1]
                        perc = float(info[2])
                        norm_total += perc                        
                        if contig in tagDic.keys():
                                if tagDic[contig] == 'Chimeric':
                                        norm_failed += perc
                                        tmp = '\t'.join(info[:2]) + '\t' + str(0.0) + '\t' + '\t'.join(info[3:])
                                        result[contig] = tmp
        for contig in contigLst:
                info = genoDic[contig].split('\t')
                if 'DLA64' not in genoDic[contig]:
                        if '50X' in genoDic[contig]:
                                contig = info[1]
                                perc = float(info[2])
                                if contig in tagDic.keys():
                                        if tagDic[contig] != 'Chimeric':
                                                if X50_total != X50_failed:
                                                        perc_new = perc * X50_total / (X50_total - X50_failed)
                                                else:
                                                        perc_new = '0.0'
                                                tmp = '\t'.join(info[:2]) + '\t' + str(perc_new) + '\t' + '\t'.join(info[3:])
                                                result[contig] = tmp
                        else:
                                contig = info[1]
                                perc = float(info[2])
                                if contig in tagDic.keys():
                                        if tagDic[contig] != 'Chimeric':
                                                if norm_total != norm_failed:
                                                        perc_new = perc * norm_total / (norm_total - norm_failed)
                                                else:
                                                        perc_new = '0.0'
                                                tmp = '\t'.join(info[:2]) + '\t' + str(perc_new) + '\t' + '\t'.join(info[3:])
                                                result[contig] = tmp
        # get sum RAF of all adjusted non-DLA64 contigs
        nonDLA64 = 0.0
        contigLst = result.keys()
        for contig in contigLst:
                if 'DLA64' not in contig:
                        RAF = float(result[contig].split('\t')[2])
                        nonDLA64 += RAF
        return result,origPerc,nonDLA64

if __name__ == '__main__':
	genotypingResult0 = sys.argv[1]
	resultPath = sys.argv[2]
        DLA64nucSeq = sys.argv[3]
	DLA88_depthCutoff = int(sys.argv[4]) 
        DLA12_depthCutoff = int(sys.argv[5])
        lengthCutoff = int(sys.argv[6]) # 270
	DLAsam = sys.argv[7]
        with open(genotypingResult0,'r') as f:
		genotypingResult = f.read()
        with open(DLAsam,'r') as f:
                DLAsam = f.read()
        lst = genotypingResult.split('\n')[:-1]
        DLAcount = len(DLAsam.split('\n')[:-1]) # read count instead of read pair
	normLst = []
	DLA64 = []
	genoDic = {}
	for i in range(len(lst)):
		contigName = lst[i].split('\t')[1]
		genoDic[contigName] = lst[i]
		contig = contigName
		if 'DLA64' not in lst[i]:
			normLst.append(contig)
		else:
			DLA64.append(contig)
	# check if norm 50X contigs exist
	out = open(genotypingResult0.split('.')[0] + '.stat','w')
	header = 'Sample\tContigName\tRelativeExpLevel\tAllele\tAlleleGroup\tGene\tProtSeq\tProtHVR\tNucSequence\tConfidence\tActuralDLAexp\tRelativeExpLevelBeforeAjustment\tActuralDLAexpBeforeAdjustment\tEvidence\tPolySiteLinkageAndSupportingReadCount\tReadPairLinkages\n'
	out.write(header)
	if normLst != []:
		contigLinkage = resultPath + '/NonDLA64_E23NucSeqAlignment_ContigLinkages'
		posReadCount = resultPath + '/PolymorphicSites_nucleotide_freq.txt'
		nucSeq = resultPath + '/NonDLA64_E23NucSeqAlignment.fa'
		readLinkage = resultPath + '/ReadPair_polySite_linkages.txt'
		pairPos = resultPath + '/Pairwise_position_counts.txt'
		casePos = resultPath + '/Pairwise_linkage.txt'
		posCountDic = pos_dic(posReadCount)
		linkageDic = linkage_dic(contigLinkage)
		seqDic = seq_dic(nucSeq)
		readAssignDic = read_assignment(readLinkage,linkageDic)
		pairPosDic = pairwise_pos_dic(pairPos)
		casePosDic = pairwise_pos_dic(casePos)
		tagDic,dic1,dic2 = contig_validation(linkageDic,pairPosDic,casePosDic,readAssignDic,DLA88_depthCutoff,DLA12_depthCutoff,seqDic,lengthCutoff,genoDic) # validation
		evidenceDic = evidence_dic(linkageDic,pairPosDic,casePosDic)
                genoDic,origPerc,nonDLA64RAF = adjust_relativeAF_basedOn_tag(genoDic,tagDic)
		for contig in normLst:
			if contig in seqDic.keys():
				seq = seqDic[contig]
			else:
				seq = ''
			if contig in linkageDic.keys():
				linkage = linkageDic[contig]
			else:
				linkage = []
			countLst = []
			for pos in linkage:
				if pos in posCountDic.keys():
					count = posCountDic[pos]
				else:
					count = 'NA'
				tmp = pos + ':' + count
				countLst.append(tmp)
			contigName = contig
			geno = genoDic[contigName]
			readLst = readAssignDic[contig].keys()
			readLinkageCount = []
			for read in readLst:
				count = readAssignDic[contig][read]
				tmp2 = read + '(' + str(count) + ')'
				readLinkageCount.append(tmp2)
			readLinkageCount = sorted(readLinkageCount,key=natural_key)
			tag = tagDic[contig]
			evid =','.join(sorted(evidenceDic[contig],key=natural_key))
                        origRAF = origPerc[contig]
                        RAF = float(geno.split('\t')[2])
                        if nonDLA64RAF != 0.0:
                                EXP = round(RAF * DLAcount / nonDLA64RAF)
                                origEXP = round(float(origRAF) * DLAcount / nonDLA64RAF)
                        else:
                                EXP = 0.0
                                origEXP = 0.0
                        string = geno + '\t' + seq + '\t' + tag + '\t' + str(EXP) + '\t' + str(origRAF) + '\t' + str(origEXP) + '\t' + evid + '\t' + ','.join(countLst) + '\t' + ','.join(readLinkageCount) + '\n'
			out.write(string)
	if DLA64 != []:
		seqDic = seq_dic(DLA64nucSeq)
		for contig in DLA64:
			geno = genoDic[contig]
			seq = seqDic[contig.split('NormDLA64')[1]]
                        origRAF = origPerc[contig]
                        allele = geno.split('\t')[3]
                        if allele != '':
                                tag = 'Known allele'
                        else:
                                tag = 'New DLA-64'
                        RAF = float(geno.split('\t')[2])
                        try:
                                nonDLA64RAF
                        except NameError:
                                nonDLA64RAF = 0.0
                        if nonDLA64RAF != 0.0:
                                EXP = round(RAF * DLAcount / nonDLA64RAF)
                                origEXP = round(float(origRAF) * DLAcount / nonDLA64RAF)
                        else:
                                EXP = 0.0
                                origEXP = 0.0
                        string = geno + '\t' + seq + '\t' + tag + '\t' + str(EXP) + '\t' + str(origRAF) + '\t' + str(origEXP) + '\t\t\t\n'
                        out.write(string)
	out.close()

