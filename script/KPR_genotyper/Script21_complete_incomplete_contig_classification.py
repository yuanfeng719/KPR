#!/usr/bin/python

import sys
import glob
import os

def contig_seq_dic(file0):
        dic = {}
        if os.path.isfile(file0):
                with open(file0,'r') as f:
                        file = f.read()
                lst = file.split('>')[1:]
                for i in range(len(lst)):
                        name = lst[i].split('\n')[0]
                        seq = ''.join(lst[i].split('\n')[1:])
                        dic[name] = seq
        return dic

def E23_sequence_dic(E23_files):
        dic = {}
        for file0 in E23_files:
                with open(file0,'r') as f:
                        file = f.read()
                name = file.split('>')[1].split('\n')[0]
                seq = ''.join(file.split('>')[1].split('\n')[1:])
                dic[name] = seq
        return dic

if __name__ == '__main__':
        complete_contig = sys.argv[1]
        E23_path = sys.argv[2]
        nonDLA64 = sys.argv[3]
        DLA64 = sys.argv[4]
        with open(complete_contig,'r') as f:
                complete_contig = f.read()
        complete_contig = complete_contig.split('\n')[:-1]
        E23_files = glob.glob(E23_path + '/*.clustalO_Exon2N3')
        E23_seqDic = E23_sequence_dic(E23_files)
        nonDLA64Dic = contig_seq_dic(nonDLA64)
        DLA64Dic = contig_seq_dic(DLA64)
        # [1] complete nonDLA64 [2] incomplete nonDLA64 [3] complete DLA64
        out = open(nonDLA64.split('.')[0] + '_complete.txt','w')
        out2 = open(nonDLA64.split('.')[0] + '_incomplete.txt','w')
        out3 = open(DLA64.split('.')[0] + '_complete.txt','w')
        # output
        nonDLA64_lst = nonDLA64Dic.keys()
        for contig in nonDLA64_lst:
                seq = nonDLA64Dic[contig]
                if contig in complete_contig:
                        E23 = E23_seqDic[contig]
                        string = '>' + contig + '\n' + E23 + '\n'
                        out.write(string)
                else:
                        string = '>' + contig + '\n' + seq + '\n'
                        out2.write(string)
        DLA64_lst = DLA64Dic.keys()
        DLA64_seqLst = []
        for contig in DLA64_lst:
                if contig in complete_contig:
                        E23 = E23_seqDic[contig]
                        if E23 not in DLA64_seqLst:
                                DLA64_seqLst.append(E23)
                                string = '>' + contig + '\n' + E23 + '\n'
                                out3.write(string)
        out.close()
        out2.close()
        out3.close()


