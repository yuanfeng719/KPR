#!/usr/bin/python

import sys

def seq_dic(nucSeq):
        dic = {}
        with open(nucSeq,'r') as f:
                nucSeq = f.read()
        lst = nucSeq.split('>')[1:]
        for i in range(len(lst)):
                contig = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                dic[contig] = seq
        return dic

if __name__ == '__main__':
        genotypingResult0 = sys.argv[1]
        resultPath = sys.argv[2]
        output_path = sys.argv[3]
        with open(genotypingResult0,'r') as f:
                genotypingResult = f.read()
        lst = genotypingResult.split('\n')[:-1]
        normLst = []
        X50Lst = []
        genoDic = {}
        for i in range(len(lst)):
                contigName = lst[i].split('\t')[1]
                genoDic[contigName] = lst[i]
                if ('Norm' in contigName) and ('DLA-64' not in lst[i]): # exclude DLA-64 alleles
                        contig = contigName.split('Norm')[1]
                        normLst.append(contig)
                elif '50X' in contigName:
                        contig = contigName.split('50X')[1]
                        X50Lst.append(contig)
        # check if norm 50X contigs exist
        out = open(output_path + '/NonDLA64_WholeNucSeq.fa','w')
        out2 = open(output_path + '/NonDLA64_E23NucSeqAlignment.fa','w')
        if normLst != []:
                nucSeq = resultPath + '/nonDLA64/combined/Correct_contigs.fa'
                seqDic = seq_dic(nucSeq)
                for contig in normLst:
                        seq = seqDic[contig]
                        new_seq = 'ATGGAGGTGGTGATGCCGCGAGCCCTCCTCGTGCTGCTGTCGGCGGCCCTGGCCGTGACCCTGACCCGGGCGG' + seq + 'AACCCCCCAGCACACGTGTGACCCGCCACCCCGTCTCTGACCGTGAGGTCACCCTGAGGTGCTGGGCGCTGGGGTTCTACCCTGAGGAGATCACCCTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACAGAGGTTGTGGACACAAGGCCTGCAGGAGATGGGACCTTCCAGAAGTGGGCGGCCGTGGTGGTGCCCTCTGGACAGGAGCAGAGATACACGTGCCACGTGCAGCATGAGGGGCTGGCGGAGCCTGTCACGCGGAGATGGGAGCCTTCCCCTCTGTCCACCATAGTCATCGTCAGCATTGCTGCTCTGGTTCTCCTCGTGGTCGCTGGGGTGATTGGAGCTGTGATCTGGAGGAAGCAGCGCTCGGGAGGAAAAGGACCAGGCTACTCTCATGCTGCACGCGATGACAGTGCCCAGGGCTCTGATGTGTCTCTGACAGCTCCTAGAGTGTGA'
                        string = '>' + 'Norm' + contig + '\n' + new_seq + '\n'
                        out.write(string)
                        if X50Lst != []:
                                new_seq2 = seq[:464] + '---' + seq[464:]
                        else:
                                new_seq2 = seq
                        string2 = '>' + 'Norm' + contig + '\n' + new_seq2 + '\n'
                        out2.write(string2)
        if X50Lst != []:
                nucSeq = resultPath + '/nonDLA64/combined/Correct_contigs.fa'
                seqDic = seq_dic(nucSeq)
                for contig in X50Lst:
                        seq = seqDic[contig]
                        new_seq = 'ATGGAGGTGGTGATGCCGCGAGCCCTCCTCGTGCTGCTGTCGGCGGCCCTGGCCGTGACCCTGACCCGGGCGG' + seq + 'AACCCCCCAGCACACGTGTGACCCGCCACCCCGTCTCTGACCGTGAGGTCACCCTGAGGTGCTGGGCGCTGGGGTTCTACCCTGAGGAGATCACCCTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACAGAGGTTGTGGACACAAGGCCTGCAGGAGATGGGACCTTCCAGAAGTGGGCGGCCGTGGTGGTGCCCTCTGGACAGGAGCAGAGATACACGTGCCACGTGCAGCATGAGGGGCTGGCGGAGCCTGTCACGCGGAGATGGGAGCCTTCCCCTCTGTCCACCATAGTCATCGTCAGCATTGCTGCTCTGGTTCTCCTCGTGGTCGCTGGGGTGATTGGAGCTGTGATCTGGAGGAAGCAGCGCTCGGGAGGAAAAGGACCAGGCTACTCTCATGCTGCACGCGATGACAGTGCCCAGGGCTCTGATGTGTCTCTGACAGCTCCTAGAGTGTGA'
                        string = '>' + '50X' + contig + '\n' + new_seq + '\n'
                        out.write(string)
                        string2 = '>' + '50X' + contig + '\n' + seq + '\n'
                        out2.write(string2)
        out.close()

