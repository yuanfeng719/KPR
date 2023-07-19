#!/usr/bin/python

import sys
import os

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

def read_in_file(file0):
        result = []
        if os.path.isfile(file0):
                with open(file0,'r') as f:
                        file = f.read()
                if '>' in file:
                        result = file.split('>')[1:]
        return result

def rename_DLA64_contig(lst):
        result = []
        for i in range(len(lst)):
                name = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                name = 'DLA64' + name
                tag = name + '\n' + seq + '\n'
                result.append(tag)
        return result

def rename_DLA79_contig(lst):
        result = []
        for i in range(len(lst)):
                name = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                name = 'DLA79' + name
                tag = name + '\n' + seq + '\n'
                result.append(tag)
        return result


######################################################
if __name__ == '__main__':
	HVR = sys.argv[1]
	exon = sys.argv[2]
        DLA64HVR = sys.argv[3]
        DLA64exon = sys.argv[4]
        DLA79HVR = sys.argv[5]
        DLA79exon = sys.argv[6]
        norm = sys.argv[7]
        X50 = sys.argv[8]
        d64 = sys.argv[9]
        d79 = sys.argv[10]
	HVRdb = sys.argv[11]
	exondb = sys.argv[12]
	sampleID = sys.argv[13]
	# output files
	out = open('Genotyping_result.txt','w')
        HVR = read_in_file(HVR)
        exon = read_in_file(exon)
        DLA64HVR = read_in_file(DLA64HVR)
        DLA64exon = read_in_file(DLA64exon)
        DLA79HVR = read_in_file(DLA79HVR)
        DLA79exon = read_in_file(DLA79exon)
	
        norm = read_in_file(norm)
        X50 = read_in_file(X50)
        d64 = read_in_file(d64)
        d79 = read_in_file(d79)

        HVRdb = read_in_file(HVRdb)
        exondb = read_in_file(exondb)
        
        # rename DLA64 contigs
        DLA64HVR = rename_DLA64_contig(DLA64HVR)
        DLA64exon = rename_DLA64_contig(DLA64exon)
        DLA79HVR = rename_DLA79_contig(DLA79HVR)
        DLA79exon = rename_DLA79_contig(DLA79exon)

        # norm x50 dla64 perc
        totalRead = len(norm) + len(X50) + len(d64) + len(d79)
        if totalRead == 0:
                totalRead = 1

        norm_perc = float(len(norm)) / totalRead
        X50_perc = float(len(X50)) / totalRead
        d64_perc = float(len(d64)) / totalRead
        d79_perc = float(len(d79)) / totalRead
        
        if len(DLA64exon) >= 1:
                d64_perc_each = d64_perc / len(DLA64exon)
        else:
                d64_perc_each = d64_perc
        if len(DLA79exon) >= 1:
                d79_perc_each = d79_perc / len(DLA79exon)
        else:
                d79_perc_each = d79_perc
	# data processing
	HVR_dic = database_dic(HVRdb)
	exon_dic = database_dic(exondb)
        totalHVR = HVR + DLA64HVR + DLA79HVR
        totalExon = exon + DLA64exon + DLA79exon
	contig_HVR_dic = contig_dic(totalHVR)
	contig_exon_dic = contig_dic(totalExon)
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
                if len(contig_exon) > 182 and 'DLA79' not in contig:
                        perc = X50_perc
                elif 'DLA64' in contig:
                        perc = d64_perc_each
                elif 'DLA79' in contig:
                        perc = d79_perc_each
                else:
                        perc = norm_perc
		string = sampleID + '\t' + contig + '\t' + str(perc) + '\t' + allele + '\t' + alleleGroup + '\t' + contig_exon + '\t' + contig_HVR + '\n'
		out.write(string)
	out.close()



