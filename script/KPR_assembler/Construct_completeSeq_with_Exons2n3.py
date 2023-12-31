#!/usr/bin/python

### this script takes the genotyping result <*.stat> as input and generate complete nucleotide sequences for each alleles ("Chimeric" excluded)

import sys

if __name__ == '__main__':
	geno0 = sys.argv[1]
	outputPath = sys.argv[2]

	out = open(outputPath + '/' + geno0.split('/')[-1].split('.')[0] + '_completeNucSeq.fa','w')

	with open(geno0,'r') as f:
		geno = f.read()
	if 'Align' in geno:
		lst = geno.split('\n')[1:-1]
		for i in range(len(lst)):
			if "Chimeric" not in lst[i]:
				info = lst[i].split('\t')
				contigName = info[1]
				allele,gp,gene = info[3:6]
				nuc = info[8]
				if allele != '':
					name = allele
				else:
					name = gene + '*' + contigName
				if 'DLA-12' in lst[i]:
					begin = 'ATGGTGCCCGGAACCCTAGCCCTGCTGCTGTCGGGGGCCCTGGCCGTGACCCTGACCCGGGCGG'
					end = 'AACCCCCCAGCACACGTGTGACCCGCCACCCCGTCTCTGACCATGAGGTCACCCTGAGGTGCTGGGCGCTGGGCTTCTACCCTGCGGAGATCACCCTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACAGAGGTTGTGGACACAAGGCCTGCAGGAGATGGGACCTTCCAGAAGTGGGCGGCCGTGGTGGTGCCCTCTGGACAGGAGCAGAGATACACGTGCCACGTGCAGCATGAGGGGCTGGCGGAGCCTGTCACGCGGAGATGGGAGCCTTCCCCTCTGTCCACCATTGTCATCGTCAGCATTGCTGCTCTGGTTCTCCTCGTGGTCGCTGGGGTGATTGGAGCTGTGATCTGGAGGAAGCAGCGCTCAGGAGGAAAAGGACCAGGCTACTCTCATGCTGCACGCGATGACAGTGCCCAGGGCTCTGATGTGTCTCTGACAGCTCCTAGAGTGTGA'
				elif 'DLA-64' in lst[i]:
					begin = 'ATGGAGGTGGTGATGCCGCGAGCCCTCCTCGTGCTGCTGTCGGCGGCCCTGGCCGTGACCCTGACCCGGGCGG'
					end = 'ATCCTCCAAAAATCTACCTGACCCGCCACCCCATCTCTGACCATGAGGTCACCCTGAGGTGCTGGGCGCTGGGCTTCTACCCTGCGGAGATCACCCTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACAGAGGTTGTGGACACCAGGCCTGCAGGAGATGGGACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCCTCTGGACAGGAGCAGAGATACACGTGCCACGTGCAGCATGAGGGGCTGGCGGAGCCTGTCACGCGGAGATGGGAGCCTTCCCCTCTGTCCACCATTGTCATCGTCAGCATTGCTGCTCTGGTTCTCCTCGTGGTCGCTGGGGTGATTGGAGCTGTGATCTGGAGGAAGCAGCGCTCGGGAGGAAAAGGACCAGGCTACTCTCATGCTGCACGTGATGACAGTGCCCAGGGCTCTGATGTGTCTCTGACAGCTCCTAGAGTGTGA'
				else:
					begin = 'ATGGAGGTGGTGATGCCGCGAGCCCTCCTCGTGCTGCTGTCGGCGGCCCTGGCCGTGACCCTGACCCGGGCGG'
					end = 'AACCCCCCAGCACACGTGTGACCCGCCACCCCGTCTCTGACCGTGAGGTCACCCTGAGGTGCTGGGCGCTGGGGTTCTACCCTGAGGAGATCACCCTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACAGAGGTTGTGGACACAAGGCCTGCAGGAGATGGGACCTTCCAGAAGTGGGCGGCCGTGGTGGTGCCCTCTGGACAGGAGCAGAGATACACGTGCCACGTGCAGCATGAGGGGCTGGCGGAGCCTGTCACGCGGAGATGGGAGCCTTCCCCTCTGTCCACCATAGTCATCGTCAGCATTGCTGCTCTGGTTCTCCTCGTGGTCGCTGGGGTGATTGGAGCTGTGATCTGGAGGAAGCAGCGCTCGGGAGGAAAAGGACCAGGCTACTCTCATGCTGCACGCGATGACAGTGCCCAGGGCTCTGATGTGTCTCTGACAGCTCCTAGAGTGTGA'
				completeSeq = begin + nuc + end
				string = '>' + name + '\n' + completeSeq + '\n'
				out.write(string)
	out.close()



