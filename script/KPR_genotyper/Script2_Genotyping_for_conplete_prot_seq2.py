#!/usr/bin/python

import sys

######################################################
def database_dic(lst):
	dic = {}
	for i in range(len(lst)):
		allele = lst[i].split('\n')[0]
		seq = ''.join(lst[i].split('\n')[1:])
		dic[seq] = allele
	return dic

def contig_dic(lst):
	dic = {}
	for i in range(len(lst)):
		contig = lst[i].split('\n')[0]
		seq = ''.join(lst[i].split('\n')[1:])
		dic[contig] = seq
	return dic

def score_dic(lst):
	dic = {}
	for i in range(len(lst)):
		contig = lst[i].split('\n')[0]
		score = lst[i].split('\n')[1]
		dic[contig] = score
	return dic

######################################################
if __name__ == '__main__':
	HVR = sys.argv[1]
	exon = sys.argv[2]
	HVRdb = sys.argv[3]
	exondb = sys.argv[4]
	score = sys.argv[5]
	sampleID = sys.argv[6]
	# output files
	out = open('Genotyping_result.txt','w')
	with open(HVR,'r') as f:
        	HVR = f.read()
	with open(exon,'r') as f:
        	exon = f.read()
	with open(HVRdb,'r') as f:
        	HVRdb = f.read()
	with open(exondb,'r') as f:
        	exondb = f.read()
	with open(score,'r') as f:
        	score = f.read()
	HVR = HVR.split('>')[1:]
	exon = exon.split('>')[1:]
	HVRdb = HVRdb.split('>')[1:]
	exondb = exondb.split('>')[1:]
	score = score.split('>')[1:]
	# data processing
	HVR_dic = database_dic(HVRdb)
	exon_dic = database_dic(exondb)
	contig_HVR_dic = contig_dic(HVR)
	contig_exon_dic = contig_dic(exon)
	score_dic = score_dic(score)
	contigs = contig_exon_dic.keys()
	for i in range(len(contigs)):
		contig = contigs[i]
		contig_HVR = contig_HVR_dic[contig]
		contig_exon = contig_exon_dic[contig]
		if contig_exon in exon_dic.keys():
			allele = exon_dic[contig_exon]
		else:
			allele = ''
		if contig_HVR in HVR_dic.keys():
			alleleGroup = HVR_dic[contig_HVR].split(':')[0]
		else:
			alleleGroup = ''
		if allele != '':
			alleleGroup = allele.split(':')[0]
		if contig in score_dic.keys():
			contig_score = score_dic[contig]
		else:
			contig_score = '0.0'
		string = sampleID + '\t' + contig + '\t' + str(contig_score) + '\t' + allele + '\t' + alleleGroup + '\t' + contig_exon + '\t' + contig_HVR + '\n'
		out.write(string)
	out.close()

