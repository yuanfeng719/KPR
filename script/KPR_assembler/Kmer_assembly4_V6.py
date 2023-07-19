#!/usr/bin/python

import sys
import random
import difflib

################################################
def build_kmer_dic(seqDic,kmer):
	dic = {}
	lst = seqDic.keys()
	for i in range(len(lst)):
		header = lst[i]
		seq = seqDic[header]
		for j in range(len(seq)-kmer+1):
			kmer_seq = seq[j:j+kmer]
			info = []
			info.append(header)
			info.append(j)
			if kmer_seq in dic:
				dic[kmer_seq].append(info)
			else:
				dic[kmer_seq] = []
				dic[kmer_seq].append(info)
	return (dic)

def check_identical(seq1,seq2,pos1,pos2,kmer):
	flag = True
	pos1_0 = pos1
	pos2_0 = pos2
	pos1_tail = pos1+kmer
	pos2_tail = pos2+kmer
	while pos1 >=0 and pos2 >=0 and flag == True:
		if seq1[pos1] != seq2[pos2]:
			flag = False
		else:
			pos1 = pos1 - 1
			pos2 = pos2 - 1
	while pos1_tail < len(seq1) and pos2_tail < len(seq2) and flag == True:
		if seq1[pos1_tail] != seq2[pos2_tail]:
			flag = False
		else:
			pos1_tail = pos1_tail + 1
			pos2_tail = pos2_tail + 1
	# assemble contig
	if flag:
		if pos1_0 >= pos2_0:
			head = seq1[:pos1_0]
		else:
			head = seq2[:pos2_0]
		body = seq1[pos1_0:pos1_0+kmer]
		if len(seq1)-(pos1_0+kmer) >= len(seq2)-(pos2_0+kmer):
			tail = seq1[pos1_0+kmer:]
		else:
			tail = seq2[pos2_0+kmer:]
		sequence = head + body + tail
	else:
		sequence = ''
	return (flag,sequence)

def build_seq_dic(lst):
	dic = {}
	for i in range(len(lst)):
		header = lst[i].split('\n')[0]
		seq = lst[i].split('\n')[1]
		dic[header] = seq
	return dic

def remove_element_from_dic_lst(value,name):
	out = []
	for i in range(len(value)):
		if value[i][0] not in name:
			out.append(value[i])
	return (out)

def Two_seq_expension(seq1_seq,seq2_seq,cutoff):
        # find match blocks
        s = difflib.SequenceMatcher(None,seq1_seq,seq2_seq,autojunk=False)
        if len(s.get_matching_blocks()) == 2:
                seq1_start,seq2_start,length = s.get_matching_blocks()[0]
                if length >= cutoff and (seq1_start+length==len(seq1_seq) or seq2_start+length==len(seq2_seq)):
                        if seq1_start==0 or seq2_start==0:
                                if seq1_start!=seq2_start and len(seq1_seq) > seq1_start+length and len(seq2_seq) > seq2_start+length and seq1_seq[seq1_start+length] == seq2_seq[seq2_start+length]:
                                        overlap = seq1_seq[seq1_start:seq1_start+length]
                                        if seq1_start >= seq2_start:
                                                front = seq1_seq[:seq1_start]
                                        else:
                                                front = seq2_seq[:seq2_start]
                                        if len(seq1_seq)-(seq1_start+length) >= len(seq2_seq)-(seq2_start+length):
                                                if len(seq1_seq)-(seq1_start+length) != 0:
                                                        back = seq1_seq[-(len(seq1_seq)-(seq1_start+length)):]
                                                else:
                                                        back = ''
                                        else:
                                                if len(seq2_seq)-(seq2_start+length) != 0:
                                                        back = seq2_seq[-(len(seq2_seq)-(seq2_start+length)):]
                                                else:
                                                        back = ''
                                        extend = front+overlap+back
                                        return (extend)
                                elif seq1_start!=seq2_start and (len(seq1_seq) == seq1_start+length or len(seq2_seq) == seq2_start+length):
                                        overlap = seq1_seq[seq1_start:seq1_start+length]
                                        if seq1_start >= seq2_start:
                                                front = seq1_seq[:seq1_start]
                                        else:
                                                front = seq2_seq[:seq2_start]
                                        if len(seq1_seq)-(seq1_start+length) >= len(seq2_seq)-(seq2_start+length):
                                                if len(seq1_seq)-(seq1_start+length) != 0:
                                                        back = seq1_seq[-(len(seq1_seq)-(seq1_start+length)):]
                                                else:
                                                        back = ''
                                        else:
                                                if len(seq2_seq)-(seq2_start+length) != 0:
                                                        back = seq2_seq[-(len(seq2_seq)-(seq2_start+length)):]
                                                else:
                                                        back = ''
                                        extend = front+overlap+back
                                        return (extend)
                                elif seq1_start == seq2_start:
					if len(seq1_seq) >= len(seq2_seq):
                                        	extend = seq1_seq
					else:
						extend = seq2_seq
                                        return (extend)
                                else:
                                        return('NA')
                        else:
                                return('NA')
                else:
                        return('NA')
        else:
                return('NA')


def expand_contig(fwSeq,rvSeq,forwdD,kmer):
	# random starting
	header = fwSeq.keys()[random.randint(0,len(fwSeq.keys())-1)]
        contigF = fwSeq[header]
	contigR = rvSeq[header]
	name = []
	name.append(header)
	# expand the forward contig from both end (head, tail)
        # head
	flag = True
	contig_head_f='NA'
	contig_tail_f='NA'
	while flag:
		if contigF[:kmer] != contig_head_f:
                	contig_head_f = contigF[:kmer]
                	head_f_value = forwdD[contig_head_f]
                	head_f_value = remove_element_from_dic_lst(head_f_value,header)
                if len(head_f_value) > 0:
                        seedF = random.randint(0,len(head_f_value)-1)
                        seedF = head_f_value[seedF]
                        header = seedF[0]
                        pos = seedF[1]
                        seqF = fwSeq[header]
                        flagT,resultF = check_identical(contigF,seqF,0,pos,kmer)
                        # check reverse strand
                        if flagT:
                                seqR = rvSeq[header]
				resultR = Two_seq_expension(contigR,seqR,kmer)
				if resultR != 'NA':
					contigF = resultF
					contigR = resultR
					name.append(header)
			head_f_value = remove_element_from_dic_lst(head_f_value,header)
		else:
			flag = False
	#tail
	flag = True
	while flag:
		if contigF[-kmer:] != contig_tail_f:
			contig_tail_f = contigF[-kmer:]
			tail_f_value = forwdD[contig_tail_f]
			tail_f_value = remove_element_from_dic_lst(tail_f_value,header)
		if len(tail_f_value) > 0:
			seedF = random.randint(0,len(tail_f_value)-1)
			seedF = tail_f_value[seedF]
			header = seedF[0]
			pos = seedF[1]
			seqF = fwSeq[header]
			flagT,resultF = check_identical(contigF,seqF,len(contigF)-kmer,pos,kmer)
			#check reverse strand
			if flagT:
				seqR = rvSeq[header]
				resultR = Two_seq_expension(contigR,seqR,kmer)
				if resultR != 'NA':
					contigF = resultF
					contigR = resultR
					name.append(header)
			tail_f_value = remove_element_from_dic_lst(tail_f_value,header)
		else:
			flag = False
	name = list(set(name))
	name.sort()
	return (contigF,contigR,name)

# filter the dic with minimum count
def filter_seq_dic(fw,rv,fwKmer,revKmer,fwCutoff,rvCutoff):
	# record all read pairs that have problem and remove them from sequence dictionary
	problem_reads = []
	#forwd
	keys = fwKmer.keys()
	for i in range(len(keys)):
		key = keys[i]
		lst = fwKmer[key]
		if len(lst) <= fwCutoff:
			for j in range(len(lst)):
				read = lst[j][0]
				if read not in problem_reads:
					problem_reads.append(read)
	#rev
	keys = revKmer.keys()
        for i in range(len(keys)):
                key = keys[i]
                lst = revKmer[key]
		if len(lst) <= rvCutoff:
			for j in range(len(lst)):
                                read = lst[j][0]
                                if read not in problem_reads:
                                        problem_reads.append(read)
	# remove problem reads from original forwd and reverse sequence dictionaries
	for i in range(len(problem_reads)):
		read = problem_reads[i]
		del fw[read]
		del rv[read]
	return (fw,rv)

# identify the first peak
def distribution_and_cutoff_determination(dic,minimumKmer):
	# the input dictionary is the kmer count dictionary
	result = {}
	values = dic.values()
	for i in range(len(values)):
		value = values[i]
		length = len(value)
		if length not in result.keys():
			result[length] = 1
		else:
			result[length] += 1
	counts = result.keys()
	counts.sort()
	# determine the cutoff to remove the seq error peak in the distribution
	length = len(counts)
	i = 0
	flag = True
	outcome = minimumKmer
	previous_freq = result[counts[0]] + 1 # starting freq, meaningless
	while flag and i < length: 
		count = counts[i]
		freq = result[count]
		if (freq > previous_freq) and (i >= minimumKmer):
			flag = False
			outcome = counts[i-1]
		else:
			previous_freq = freq
			i += 1
	if minimumKmer > length:
		outcome = 3
	return (outcome)
		
	
################################################
if __name__ == '__main__':
	forward = sys.argv[1]
	reverse = sys.argv[2]
	cutoff = int(sys.argv[3])
	N_number = int(sys.argv[4]) # repeat times
	minimumKmer = int(sys.argv[5])
	with open(forward,'r') as f:
		forwd = f.read()
	with open(reverse,'r') as f:
		rev = f.read()
	forwd = forwd.split('>')[1:]
	rev = rev.split('>')[1:]
	### data processing
	# sequence dictionary
	fw = build_seq_dic(forwd)
	rv = build_seq_dic(rev)
	# kmer dictionary
	fwKmer = build_kmer_dic(fw,cutoff)
	revKmer = build_kmer_dic(rv,cutoff)
	# summerize Kmer count distribution
	out = open(forward + '_KmerCounts_step1.txt','w')
	keys = fwKmer.keys()
	for i in range(len(keys)):
		key = keys[i]
		count = len(fwKmer[key])
		string = key + '\t' + str(count) + '\n'
		out.write(string)
	out.close()
	# generate the distribution of Kmer counts and select the cutoff to remove the first peak caused by sequencing error
	fwCutoff = distribution_and_cutoff_determination(fwKmer,minimumKmer)
	rvCutoff = distribution_and_cutoff_determination(revKmer,minimumKmer)
	# filter the sequence dictionaty with the minimum Kmer limit 
	fw,rv = filter_seq_dic(fw,rv,fwKmer,revKmer,fwCutoff,rvCutoff)
	# new kmer dictionary
	fwKmer = build_kmer_dic(fw,cutoff)
	revKmer = build_kmer_dic(rv,cutoff)
	# summerize Kmer count distribution
	out2 = open(forward + '_KmerCounts_step2.txt','w')
	keys = fwKmer.keys()
	for i in range(len(keys)):
        	key = keys[i]
        	count = len(fwKmer[key])
       		string = key + '\t' + str(count) + '\n'
        	out2.write(string)
	out2.close()
	# expand contage
	NAME = []
	for i in range(N_number):
		contigF,contigR,name = expand_contig(fw,rv,fwKmer,cutoff)
		if name not in NAME:
			flag = True
			j = 0
			while j < len(NAME) and flag:
				if set(name).issubset(NAME[j]):
					flag = False
				j+=1
			if flag:
				NAME.append(name)
				print '>'+';'.join(name)
				print contigF
				print contigR

