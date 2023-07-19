#!/usr/bin/python

import sys
import difflib

#####################################################
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

def compare_n_merge_seq(seqf,seqr,cutoff):
        s = difflib.SequenceMatcher(None,seqf,seqr,autojunk=False)
	contig = 'Mismatch'
        if len(s.get_matching_blocks()) == 2:
                seqf_start,seqr_start,length = s.get_matching_blocks()[0]
                # 2 contained cases which overlap size doesnot metter
		if length == len(seqf):
			contig = seqr
		elif length == len(seqr):
			contig = seqf
		if length >= cutoff:
                        # consider 2 other different cases
                        if seqf_start + length == len(seqf):
                                contig = Two_seq_expension(seqf,seqr,cutoff)
                        elif seqr_start + length == len(seqr):
                                contig = Two_seq_expension(seqf,seqr,cutoff)
                else:
			contig = 'Insufficient_overlap'
        else:
		contig = 'Mismatch'
	return (contig)

###################################################################
if __name__ == '__main__':
	file_in = sys.argv[1]
	cutoff = int(sys.argv[2])
	with open(file_in,'r') as f:
		file=f.read()
	lst = file.split('>')[1:]
	# output files
	match = open(file_in + '-Match','w')
	mismatch = open(file_in + '-Mismatch','w')
	InsufficientOverlap = open(file_in + '-InsufficientOverlap','w')
	# data processing
	for i in range(len(lst)):
		name = lst[i].split('\n')[0]
		forwd = lst[i].split('\n')[1]
		rev = lst[i].split('\n')[2]
		contig = compare_n_merge_seq(forwd,rev,cutoff)
		if contig == 'Insufficient_overlap':
			string = '>' + name + '\n' + forwd + '\n' + rev + '\n'
			InsufficientOverlap.write(string)
		elif contig == 'Mismatch':
			string = '>' + name + '\n' + forwd + '\n' + rev + '\n'
			mismatch.write(string)
		else:
			string = '>' + name + '\n' + contig + '\n'
			match.write(string)
	match.close()
	mismatch.close()
	InsufficientOverlap.close()

