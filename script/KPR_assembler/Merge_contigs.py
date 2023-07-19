#!/usr/bin/python

import sys
from Step3_Compare_for_rev_V2 import *

def contig_merging(contig_dict, cutoff):
    """ """
    seq_list = contig_dict.values()
    seq_list0 = seq_list[:]
    flag_progress = True
    while flag_progress:
        flag_progress = False
        result = []
        index_list = range(len(seq_list))
        for i in range(len(seq_list)-1):
            seq1 = seq_list[i]
            for j in range(i+1, len(seq_list)):
                seq2 = seq_list[j]
                extended_seq = Two_seq_expension(seq1, seq2, cutoff)
                if extended_seq != 'NA':
                    result.append(extended_seq)
                    if i in index_list: index_list.remove(i)
                    if j in index_list: index_list.remove(j)
                    flag_progress = True
        for idx in index_list:
            result.append(seq_list[idx])
        result = list(set(result))
        seq_list = result
        unmerged = []
        for seq0 in seq_list0:
            flag = True
            for seq in seq_list:
                if seq0 in seq:
                    flag = False
            if flag: unmerged.append(seq0)
        seq_list = seq_list + unmerged
    return seq_list


if __name__ == '__main__':
    file0 = sys.argv[1]
    cutoff = int(sys.argv[2])
    with open(file0, 'r') as f:
        file = f.read()
    lst = file.split('>')[1:]
    contig_dict = {}
    for i in range(len(lst)):
        contig = lst[i].split('\n')[0]
        seq = ''.join(lst[i].split('\n')[1:])
        contig_dict[contig] = seq
    seq_list = contig_merging(contig_dict, cutoff)
    
    out = open(file0.split('.')[0] + '_merged.fa', 'w')
    idx = 1
    for seq in seq_list:
        string = '>Align_' + str(idx) + '\n' + seq + '\n'
        out.write(string)
        idx += 1
    out.close()

