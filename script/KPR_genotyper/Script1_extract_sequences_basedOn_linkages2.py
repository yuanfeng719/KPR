#!/usr/bin/python

import sys
import re

############################
def contig_linkage_lst(linkage_to_contig):
	dic = {}
	for i in range(len(linkage_to_contig)):
		header = linkage_to_contig[i].split('\n')[0]
		lst = linkage_to_contig[i].split('\n')[1].split('\t')
		lst.sort()
		dic[header] = lst
	return dic

def contig_sequence(contig_seq):
	dic = {}
	for i in range(len(contig_seq)):
		header = contig_seq[i].split('\n')[0]
		seq = ''.join(contig_seq[i].split('\n')[1:])
		dic[header] = seq
	return dic

def remove_values_from_list(the_list, val):
	return [value for value in the_list if value != val]

def frequency_dic(nucleotide_frequencies):
	dic = {}
	for i in range(len(nucleotide_frequencies)):
		pos = int(nucleotide_frequencies[i].split('\t')[0])
		A = int(nucleotide_frequencies[i].split('\t')[1])
		C = int(nucleotide_frequencies[i].split('\t')[2])
		G = int(nucleotide_frequencies[i].split('\t')[3])
		T = int(nucleotide_frequencies[i].split('\t')[4])
		total = float(A + C + G + T)
		lst = [A/total,C/total,G/total,T/total]
		dic[pos] = lst
	return dic

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	return [ atoi(c) for c in re.split(r'(\d+)', text) ]

################################################################
if __name__ == '__main__':
	contig = sys.argv[1]
	contig_seq = sys.argv[2]
	nucleotide_frequencies = sys.argv[3]
	with open(contig,'r') as f:
        	contig = f.read()
	with open(contig_seq,'r') as f:
        	contig_seq = f.read()
	with open(nucleotide_frequencies,'r') as f:
        	nucleotide_frequencies = f.read()
	out = open('Correct_contigs.fa','w')
	out2 = open('Correct_contigs_scores.txt','w')
	contig = contig.split('>')[1:]
	contig_seq = contig_seq.split('>')[1:]
	nucleotide_frequencies = nucleotide_frequencies.split('\n')[1:-1]
	contig_linkage_dic = contig_linkage_lst(contig)
	contig_seq_dic = contig_sequence(contig_seq)
	nuc_freq_dic = frequency_dic(nucleotide_frequencies)
	contigs = []
	for i in range(len(contig)):
		align = contig[i].split('\n')[0]
		contigs.append(align)
	for i in range(len(contigs)):
		contig = contigs[i]
		print contig
		lst = contig_linkage_dic[contig]
		print '\t'.join(lst)
		lst = remove_values_from_list(lst,'')
		seq = contig_seq_dic[contig]
		# score for each linkage
		score = 0.0
		for k in range(len(lst)):
			rec = lst[k]
			pos = int(re.findall(r'\d+',rec)[0])
			nuc = re.findall(r'\D+',rec)[0]
			if pos in nuc_freq_dic.keys():
				freq_lst = nuc_freq_dic[pos]
				if nuc == 'A':
					score += freq_lst[0]
				elif nuc == 'C':
					score += freq_lst[1]
				elif nuc == 'G':
					score += freq_lst[2]
				elif nuc == 'T':
					score += freq_lst[3]
		string = '>' + contig + '\n' + seq + '\n'
		out.write(string)
		link = lst.sort(key=natural_keys)
		string2 = '>' + contig + '\n' + str(score) + '\n' + '\t'.join(lst) + '\n'
		out2.write(string2)
	out.close()
	out2.close()

