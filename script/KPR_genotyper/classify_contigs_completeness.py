#!/usr/bin/python

import sys

#############
def correct_lst(correct_contigs):
	result = []
	dic = {}
	for i in range(len(correct_contigs)):
		contig = correct_contigs[i].split('\n')[0]
		info = correct_contigs[i].split('\n')[1]
		result.append(contig)
		dic[contig] = info
	return (result,dic)

def seq_dic(contigs):
	dic = {}
	for i in range(len(contigs)):
		contig = contigs[i].split('\n')[0]
		seq = ''.join(contigs[i].split('\n')[1:])
		dic[contig] = seq
	return (dic)

def classify_completeness(corrected,dic,complete):
	complete_dic = {}
	incomplete_dic = {}
	for i in range(len(corrected)):
		contig = corrected[i]
		seq = dic[contig]
		if contig in complete:
			complete_dic[contig] = seq
		else:
			incomplete_dic[contig] = seq
	return (complete_dic,incomplete_dic)

def contigGroups_dic(contigGroupsAndOutliers):
	dic = {}
	for i in range(len(contigGroupsAndOutliers)):
		linkage = contigGroupsAndOutliers[i].split('\n')[0]
		contig = contigGroupsAndOutliers[i].split('\n')[1]
		dic[contig] = linkage
	return dic

###########
if __name__ == '__main__':
	correct_contigs0 = sys.argv[1]
	contigs = sys.argv[2]
	complete = sys.argv[3]
	contigGroupsAndOutliers0 = sys.argv[4]
	with open(correct_contigs0,'r') as f:
        	correct_contigs = f.read()
	with open(contigs,'r') as f:
        	contigs = f.read()
	with open(complete,'r') as f:
        	complete = f.read()
	with open(contigGroupsAndOutliers0,'r') as f:
        	contigGroupsAndOutliers = f.read()
	correct_contigs = correct_contigs.split('>')[1:]
	contigs = contigs.split('>')[1:]
	complete = complete.split('\n')[:-1]
	contigGroupsAndOutliers = contigGroupsAndOutliers.split('>')[1:]
	# output files
	out = open('Complete_contigs.fa','w')
	out2 = open('Incomplete_contigs.fa','w')
	out3 = open(correct_contigs0 + '_incomplete','w')
	out4 = open(contigGroupsAndOutliers0 + '_incomplete','w')
	# data processing
	corrected,inc_dic = correct_lst(correct_contigs)
	dic = seq_dic(contigs)
	complete_dic,incomplete_dic = classify_completeness(corrected,dic,complete)
	contigGroupsAndOutliers_dic = contigGroups_dic(contigGroupsAndOutliers)
	complete_contigs = complete_dic.keys()
	for i in range(len(complete_contigs)):
		contig = complete_contigs[i]
		seq = complete_dic[contig]
		string = '>' + contig + '\n' + seq + '\n'
		out.write(string)
	out.close()
	incomplete_contigs = incomplete_dic.keys()
	for i in range(len(incomplete_contigs)):
        	contig = incomplete_contigs[i]
        	seq = incomplete_dic[contig]
        	string = '>' + contig + '\n' + seq + '\n'
        	out2.write(string)
	out2.close()
	for i in range(len(incomplete_contigs)):
        	contig = incomplete_contigs[i]
		info = inc_dic[contig]
		string = '>' + contig + '\n' + info + '\n'
		out3.write(string)
	out3.close()
	linkage_contigs = contigGroupsAndOutliers_dic.keys()
	for i in range(len(linkage_contigs)):
		contig_lst = linkage_contigs[i].split(';')
		flag = False
		for j in range(len(incomplete_contigs)):
			contig = incomplete_contigs[j]
			if contig in contig_lst: # list contains incomplete contigs
				flag = True
		if flag:
			linkage = contigGroupsAndOutliers_dic[linkage_contigs[i]]
			string = '>' + linkage + '\n' + linkage_contigs[i] + '\n'
			out4.write(string)
	out4.close()

