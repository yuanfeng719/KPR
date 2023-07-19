#!/usr/bin/python

import sys
import random

#############################################
def assign_scores_for_row(last_row,residue,sequence,match,mismatch,gap):
        output = [0]*(len(sequence)+1)
        output[0] = last_row[0] + gap
        for i in range(1,len(sequence)+1):
                score = []
                if residue == sequence[i-1]:
                        diag = last_row[i-1] + match
                else:
                        diag = last_row[i-1] + mismatch
                score.append(diag)
                vert = last_row[i] + gap
                score.append(vert)
                horiz = output[i-1] + gap
                score.append(horiz)
                output[i] = max(score)
        return (output)

def generate_matrix(seq1,seq2,match,mismatch,gap):
        matrix=[]
        first_row = [0]*(len(seq1)+1)
        for i in range(len(first_row)):
                first_row[i] = first_row[i] + gap*i
        matrix.append(first_row)
        for j in range(len(seq2)):
                lst = assign_scores_for_row(matrix[-1],seq2[j],seq1,match,mismatch,gap)
                matrix.append(lst)
        return (matrix)

def trace_back(matrix,seq1,seq2,match,mismatch,gap):
        ind1 = len(seq1)
        ind2 = len(seq2)
        seq1_out = seq1
        seq2_out = seq2
        while ind1 != 0 or ind2 != 0:
                pointer = matrix[ind2][ind1]
                vert = matrix[ind2-1][ind1]
                horiz = matrix[ind2][ind1-1]
                diag = matrix[ind2-1][ind1-1]
                path = []
                if ind1 > 0 and ind2 > 0:
                        seq1_residue = seq1[ind1-1]
                        seq2_residue = seq2[ind2-1]
                        ### check ###
                        if vert == pointer - gap:
                                path.append('vert')
                        if horiz == pointer - gap:
                                path.append('horiz')
                        if seq1_residue == seq2_residue and diag == pointer - match:
                                path.append('diag')
                        if seq1_residue != seq2_residue and diag == pointer - mismatch:
                                path.append('diag')
                        final_path = path[random.randint(0,len(path)-1)]
                        if final_path == 'vert':
                                seq1_out = seq1_out[:ind1]+'-'+seq1_out[ind1:]
                                ind2 = ind2 - 1
                        elif final_path == 'horiz':
                                seq2_out = seq2_out[:ind2]+'-'+seq2_out[ind2:]
                                ind1 = ind1 - 1
                        elif final_path == 'diag':
                                ind1 = ind1 - 1
                                ind2 = ind2 - 1
                elif ind1 == 0:
                        seq1_out = '-'*ind2 + seq1_out
                        ind2 = 0
                elif ind2 == 0:
                        seq2_out = '-'*ind1 + seq2_out
                        ind1 = 0
        return (seq1_out,seq2_out)

def block_size(string):
	maximum = 0
	size = 0
	for i in range(len(string)):
		if string[i] != '-':
			size += 1
		else:
			if size > maximum:
				maximum = size
			size = 0
	return (maximum)

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
	return complement(s[::-1])

def longest_frame(prot):
	lst = prot.split('*')
	seq = ''
	length = 0
	for i in range(len(lst)):
		if len(lst[i]) > length:
			length = len(lst[i])
			seq = lst[i]
	return (seq)

def translation(seq,dic):
	lst = []
	start = [0,1,2]
	for ind in start:
		prot = ''
		for j in range((len(seq)-ind)/3):
			codon = seq[3*j+ind:3*j+3+ind]
			AA = dic[codon]
			prot += AA
		frac = longest_frame(prot)
		lst.append(frac)
	seq = reverse_complement(seq)
	for ind in start:
                prot = ''
                for j in range((len(seq)-ind)/3):
                        codon = seq[3*j+ind:3*j+3+ind]
                        AA = dic[codon]
                        prot += AA
		frac = longest_frame(prot)
		lst.append(frac)
	for seq in lst:
		if ('GSHSL' in seq) or ('GSHT' in seq):
			lst = [seq]
	return (lst)

#########################################
if __name__ == '__main__':
	file1=sys.argv[1]
	with open(file1,'r') as f:
        	nuc = f.read()
	contigs = nuc.split('>')[1:]
	out = open(file1+'_ProteinSeq','w')
	# codon table
	bases = ['T', 'C', 'A', 'G']
	codons = [a+b+c for a in bases for b in bases for c in bases]
	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	codon_table = dict(zip(codons, amino_acids))
	# 'X' to represent codons with 'N'
	nuc = ['A','T','C','G','N']
	codons2 = [a+b+c for a in nuc for b in nuc for c in nuc]
	codons3 = []
	for i in range(len(codons2)):
		if 'N' in codons2[i]:
			codons3.append(codons2[i])
	for i in range(len(codons3)):
		codon_table[codons3[i]] = 'X'
	# reference DLA-88*001:01 DLA-64*02:3
	RefSequence = 'MEVVMPRALLVLLSAALAVTLTRAGSHSLRYFYTSVSRPGRGDPRFIAVGYVDDTQFVRFDSDAATGRMEPRAPWMEQEGPEYWDRETRTAKETAQRYRVDLDTLRGYYNQSEAGSHTRQTMYGCDLGPGGRLLRGYSQDAYDGADYIALNEDLRSWTAADTAAQITRRKWEAAGTAEHDRNYLETTCVEWLRRYLEMGKETLLRAEPPSTRVTRHPVSDREVTLRCWALGFYPEEITLTWQRDGEDQTQDTEVVDTRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLAEPVTRRWEPSPLSTIVIVSIAALVLLVVAGVIGAVIWRKQRSGGKGPGYSHAARDDSAQGSDVSLTAPRV'
	RefSequence2 = 'MEVVMPRALLVLLSAALAVTLTRAGSHSLRFFHTAVSRPGRGDPLYISVGYVDDTQFLRFNSDAASPKVEPRARWMEQEGPEFWEEQTEIAKVHAQTSRSNLQTALGYYNQSEAGSHTFQWTSGCDVGPDGRLLRGYEQFAYDGADYLALDEDLRSWTAADAAAQITRRKWEAAGAAQYYRVYLQGECVQSLLKYLERGKETLQRTDPPKIYLTRHPISDHEVTLRCWALGFYPAEITLTWQRDGEDQTQDTEVVDTRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLAEPVTRRWEPSPLSTIVIVSIAALVLLVVAGVIGAVIWRKQRSGGKGPGYSHAARDDSAQGSDVSLTAPRV'
	for i in range(len(contigs)):
		header = contigs[i].split('\n')[0]
		seq = ''.join(contigs[i].split('\n')[1:])
		seq_lst = translation(seq,codon_table)
		size_lst = []
		print i
		print header
		print seq
		print '\t'.join(seq_lst)
		for prot in seq_lst:
			matr = generate_matrix(prot,RefSequence,5,0,-20)
			tup = trace_back(matr,prot,RefSequence,5,0,-20)
			size = block_size(tup[0])
			size_lst.append(size)
		print '### size 1###'
		print '\t'.join(str(b) for b in size_lst)
		size_lst2 = []
		for prot in seq_lst:
       	        	matr = generate_matrix(prot,RefSequence2,5,0,-20)
                	tup = trace_back(matr,prot,RefSequence2,5,0,-20)
                	size = block_size(tup[0])
                	size_lst2.append(size)
		print '### size 2###'
		print '\t'.join(str(b) for b in size_lst2)
		if len(size_lst) > 0 and len(size_lst2) > 0:
			if max(size_lst) >= max(size_lst2):
				maximum_ind = [i for i, j in enumerate(size_lst) if j == max(size_lst)]
				max_seq = seq_lst[maximum_ind[0]]
				if len(max_seq) >= 2:
					if max_seq[0] == 'S':
						max_seq = 'G' + max_seq
				string = '>' + header + '\n' + max_seq + '\n'
				out.write(string)
				print 'both'
			else:
				maximum_ind = [i for i, j in enumerate(size_lst2) if j == max(size_lst2)]
                		max_seq = seq_lst[maximum_ind[0]]
				if len(max_seq) >= 2:
                                	if max_seq[0] == 'S':
                                        	max_seq = 'G' + max_seq
                		string = '>' + header + '\n' + max_seq + '\n'
                		out.write(string)
				print 'both'
		elif len(size_lst) > 0:
			maximum_ind = [i for i, j in enumerate(size_lst) if j == max(size_lst)]
			max_seq = seq_lst[maximum_ind[0]]
			if len(max_seq) >= 2:
                		if max_seq[0] == 'S':
                        		max_seq = 'G' + max_seq
			string = '>' + header + '\n' + max_seq + '\n'
			out.write(string)
			print 'lst1'
		elif len(size_lst2) > 0:
			maximum_ind = [i for i, j in enumerate(size_lst2) if j == max(size_lst2)]
			max_seq = seq_lst[maximum_ind[0]]
			if len(max_seq) >= 2:
                        	if max_seq[0] == 'S':
                                	max_seq = 'G' + max_seq
			string = '>' + header + '\n' + max_seq + '\n'
			out.write(string)
			print 'lst2'
	out.close()

