#!/usr/bin/python

import sys

def read_other_dla_elements(file0):
    """ """
    with open(file0, 'r') as f:
        file = f.read()
    lst = file.split('>')[1:]
    seq_dict = {}
    for i in range(len(lst)):
        gene, element = lst[i].split('\n')[0].split('_')
        sequence = ''.join(lst[i].split('\n')[1:])
        if gene not in seq_dict.keys(): seq_dict[gene] = {}
        seq_dict[gene][element] = sequence.upper()
    return seq_dict

def db_allele_list(db):
    """ """
    with open(db, 'r') as f:
        file = f.read()
    lst = file.split('>')[1:]
    result = []
    for i in range(len(lst)):
        result.append(lst[i].split('\n')[0])
    return result

if __name__ == '__main__':
    fasta0 = sys.argv[1]
    other_elements = '/scratch/yf94402/Pan_cancer/script/HLA-HD/create_database/reference_sequences.fa'
    db = '/scratch/yf94402/Canine_MHC-I_landscape/source/KPR_DBupdated_07162023/reference/DLA_Nucleotide_DB_GenBank_AllExons.fa'
    with open(fasta0, 'r') as f:
        file = f.read()
    gene = fasta0.split('/')[-1].split('_')[0]
    ele_dict = read_other_dla_elements(other_elements)
    known_alleles = db_allele_list(db)

    out = open(fasta0.split('.')[0] + '_cleaned.fasta', 'w')
    lst = file.split('>')[1:]
    for i in range(len(lst)):
        genbank_id = lst[i].split('\n')[0].split(' ')[0]
        nuc = ''.join(lst[i].split('\n')[1:])

        if 'GCTCCCACT' in nuc:
            exon2_start = nuc.index('GCTCCCACT')
        elif 'GCTCGCACT' in nuc:
            exon2_start = nuc.index('GCTCGCACT')
        elif 'GTTCCCACT' in nuc:
            exon2_start = nuc.index('GTTCCCACT')
        else:
            exon2_start = 0
                
        if 'GGTCTCACA' in nuc:
            exon3_start = nuc.index('GGTCTCACA')
        elif 'GGGCTCACA' in nuc:
            exon3_start = nuc.index('GGGCTCACA')
        else:
            exon3_start = 'NA'
            print (lst[i])

        if 'GCGCGCAG' in nuc:
            exon3_end = nuc.index('GCGCGCAG')
        elif 'GCGCACAG' in nuc:
            exon3_end = nuc.index('GCGCACAG')
        elif 'GCTCACAG' in nuc:
            exon3_end = nuc.index('GCTCACAG')
        else:
            exon3_end = 'NA'
            print (lst[i])

        exon2_length = 270

        exon2 = nuc[exon2_start: exon2_start+exon2_length]
        exon3 = nuc[exon3_start: exon3_end+8]

        head = ['exon1']
        tail = ['exon4', 'exon5', 'exon6', 'exon7', 'exon8']
        if gene == 'DLA-88L': gene = 'DLA-88'
        head_seq = ''
        tail_seq = ''
        for ele in head:
            head_seq += ele_dict[gene][ele]
        for ele in tail:
            tail_seq += ele_dict[gene][ele]
        if genbank_id not in known_alleles:
            string = '>' + genbank_id + '\n' + head_seq + exon2 + exon3 + tail_seq + '\n'
            out.write(string)
    out.close()

