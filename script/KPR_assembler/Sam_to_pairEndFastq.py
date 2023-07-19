#!/usr/bin/python

import sys

def sam_dic(sam):
	dic = {}
	for i in range(len(sam)):
		info = sam[i].split('\t')
		read = info[0]
		tag = int(info[8])
		if read not in dic.keys():
			dic[read] = ['','']
		if tag > 0:
			dic[read][0] = sam[i]
		elif tag < 0:
			dic[read][1] = sam[i]
	return dic

if __name__ == '__main__':
	sam0 = sys.argv[1]
	outputPath = sys.argv[2]
	with open(sam0,'r') as f:
		sam = f.read()
	sam = sam.split('\n')[:-1]
	samDic = sam_dic(sam)
	# output
	out = open(outputPath + '/' + sam0.split('.')[0] + '_1.fastq','w')
	out2 = open(outputPath + '/' + sam0.split('.')[0] + '_2.fastq','w')
	readLst = samDic.keys()
	for read in readLst:
		forwd,rev = samDic[read]
		if (forwd != '') and (rev != ''):
			forwd_info = forwd.split('\t')
			rev_info = rev.split('\t')
			forwd_seq = forwd_info[9]
			rev_seq = rev_info[9]
			forwd_qual = forwd_info[10]
			rev_qual = rev_info[10]
			if (len(forwd_seq) == len(forwd_qual)) and (len(rev_seq) == len(rev_qual)):
				forwd_string = '@' + read + '\n' + forwd_seq + '\n' + '+' + '\n' + forwd_qual + '\n'
				rev_string = '@' + read + '\n' + rev_seq + '\n' + '+' + '\n' + rev_qual + '\n'
				out.write(forwd_string)
				out2.write(rev_string)
	out.close()
	out2.close()

