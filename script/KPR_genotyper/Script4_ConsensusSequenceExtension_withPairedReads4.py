#!/usr/bin/python

import sys
import re
import random
import difflib

######################################################
# generate pseudo alignment for reads
def pseudo_align(seq,start,total_len):
	count = 0
	result = ''
	while count < total_len:
		if count < start:
			result += '-'
		elif count >= start and count < start + len(seq):
			result += seq[count-start]
		elif count >= start + len(seq):
			result += '-'
		count += 1
	return (result)

# Index pairing between consensus and references based on clustalO result
def Index_pairing(query,ref):
        dic = {}
        idxQ = 0
        idxR = 0
        for i in range(len(query)):
                if query[i] != '-':
                        dic[idxQ] = idxR
                        idxQ += 1
                if ref[i] != '-':
                        idxR += 1
        return (dic)

def Reformat_pairing(lst):
        for i in range(len(lst)):
                header = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                if header == 'ConsensusSequence':
                        query_seq = seq
                elif 'DLA-88*5' in header:
                        seq_50X = seq
                elif 'DLA-88*0' in header:
                        seq_norm = seq
                else:
                        seq_D12 = seq
        # pairing between query and reference sequences
        dic_norm = Index_pairing(query_seq,seq_norm)
        dic_50X = Index_pairing(query_seq,seq_50X)
        dic_D12 = Index_pairing(query_seq,seq_D12)
        return (dic_norm,dic_50X,dic_D12)

# remove all occurance of a character from a list
def remove_values_from_list(the_list, val):
        while val in the_list:
                the_list.remove(val)

# identify the element with highest occurance in a list
def find_majority(k):
        myMap = {}
        maximum = ( '', 0 ) # (occurring element, occurrences)
        for n in k:
                if n in myMap: myMap[n] += 1
                else: myMap[n] = 1
                if myMap[n] > maximum[1]: maximum = (n,myMap[n])
        return maximum

# check more-than-once occurance
def percentMore_occur(lst,percentage):
        remove_values_from_list(lst, '-')
	remove_values_from_list(lst, 'N')
	unique = list(set(lst))
	length = len(lst)
        out = []
        for nuc in unique:
                count = 0
                for nuc2 in lst:
                        if nuc == nuc2:
                                count += 1
                if float(count)/length >= float(percentage):
                        out.append(nuc)
        if len(out) >= 2:
                flag = True
        else:
                flag = False
        return flag,out

#trim sequence to remove the first and last few nuc from polymorphic sites calling
def trim_aligned_seq(lst,cutoff):
	result = []
	for j in range(len(lst)):
		string = lst[j]
		dic = {}
		real = 0
		for i in range(len(string)):
			if string[i] != '-':
				dic[real] = i
				real += 1
		range_lst = range(min(dic.keys())+cutoff,max(dic.keys())-cutoff+1)
		new_lst = []
		for i in range_lst:
			new_lst.append(dic[i])
		new = ''
		for i in range(len(string)):
			if i not in new_lst:
				new += '-'
			else:
				new += string[i]
		result.append(new)
	return (result)

# generate consensus sequence and identify polymorphic sites
def consensus_polymorphic(seq_lst):
        # sequence length
        length = len(seq_lst[0])
        consensus = ''
        polymorphic_site = {}
	new_lst = trim_aligned_seq(seq_lst,5)
        for i in range(length):
                lst = []
		lst2 = [] # without removal to make consensus
                for j in range(len(new_lst)):
                        lst.append(new_lst[j][i])
			lst2.append(seq_lst[j][i])
                # remove all '-'
                remove_values_from_list(lst,'-')
		remove_values_from_list(lst2,'-')
		#if lst == []:
		#	lst = ['N']
                # element with highest occurance
                consensus += find_majority(lst2)[0]
                # check if polymorphic #cutoff >= 5%
                flag,poli = percentMore_occur(lst,0.05)
                if flag:
                        polymorphic_site[i] = poli
        return consensus,polymorphic_site

# remove '-' from two tails of a string(sequence)
def remove_tail_char(string,char):
        flag = True
        i = 0
        while flag and i < len(string):
                if string[i] == char:
                        i += 1
                else:
                        flag = False
                        string = string[i:]
        string = string[::-1]
        flag = True
        i = 0
        while flag and i < len(string):
                if string[i] == char:
                        i += 1
                else:
                        flag = False
                        string = string[i:]
        string = string[::-1]
        return (string)

# replace all occurance of char in string
def replace_all_in_string(string,char,replace):
	result = ''
	for i in range(len(string)):
		if string[i] == char:
			result += replace
		else:
			result += string[i]
	return result

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# analyzing reads mapped in uncovered regions and identify new polymorphic sites
def reads_polimorphicSites(forwd,rev,alignment,patterns,dic_norm,dic_50X,dic_D12,start_norm,start_50X,start_D12,end_norm,end_50X,end_D12,patternFreq):
	length = len(''.join(alignment[0].split('\n')[1:]))
	pattern_seq = ''.join(patterns[0].split('\n')[1:])
	range_norm = [dic_norm[0],dic_norm[len(pattern_seq)-1]]
	range_50X = [dic_50X[0],dic_50X[len(pattern_seq)-1]]
	range_D12 = [dic_D12[0],dic_D12[len(pattern_seq)-1]]
	# pattern Frequency
	patSorted = []
	patternFreq.sort(key=natural_keys,reverse = True)
	for i in range(len(patternFreq)):
		#freq = int(patternFreq[i].split('\t')[0])
		pat = '\t'.join(patternFreq[i].split('\t')[1:])
		patSorted.append(pat)
	for i in range(len(alignment)):
                header = alignment[i].split('\n')[0]
                seq = ''.join(alignment[i].split('\n')[1:])
                ref_seq = 'N'*length
                if 'DLA-88*5' in header:
                        seq_50X = seq
                elif 'DLA-88*0' in header:
                        seq_norm = seq
                else:
                        seq_D12 = seq
	ref_norm = Index_pairing(seq_norm,ref_seq)
        ref_50X = Index_pairing(seq_50X,ref_seq)
        ref_D12 = Index_pairing(seq_D12,ref_seq)	
	target_reads = [] # allpseudo aligned reads covering the two gap regions
	forwd_dic = {}
	rev_dic = {}
	read_lst = forwd
	for i in range(len(read_lst)):
		read = read_lst[i].split('\t')[0]
		allele = read_lst[i].split('\t')[2]
		start = int(read_lst[i].split('\t')[3])
		seq = read_lst[i].split('\t')[9]
		if 'DLA-88*5' in allele:
			if (start-1 <= range_50X[0] and start-2+len(seq) >= start_50X) or (start-1 <= end_50X and start-2+len(seq) >= range_50X[1]):
				aligned_read = pseudo_align(seq,ref_50X[start-1],length)
				target_reads.append(aligned_read)
				forwd_dic[read] = [ref_50X[start-1],aligned_read]
		elif 'DLA-12' in allele:
			if (start-1 <= range_D12[0] and start-2+len(seq) >= start_D12) or (start-1 <= end_D12 and start-2+len(seq) >= range_D12[1]):
				aligned_read = pseudo_align(seq,ref_D12[start-1],length)
                                target_reads.append(aligned_read)
				forwd_dic[read] = [ref_D12[start-1],aligned_read]
		else:
			if (start-1 <= range_norm[0] and start-2+len(seq) >= start_norm) or (start-1 <= end_norm and start-2+len(seq) >= range_norm[1]):
				aligned_read = pseudo_align(seq,ref_norm[start-1],length)
                                target_reads.append(aligned_read)
				forwd_dic[read] = [ref_norm[start-1],aligned_read]
	read_lst = rev
	for i in range(len(read_lst)):
                read = read_lst[i].split('\t')[0]
                allele = read_lst[i].split('\t')[2]
                start = int(read_lst[i].split('\t')[3])
                seq = read_lst[i].split('\t')[9]
                if 'DLA-88*5' in allele:
                        if (start-1 <= range_50X[0] and start-2+len(seq) >= start_50X) or (start-1 <= end_50X and start-2+len(seq) >= range_50X[1]):
                                aligned_read = pseudo_align(seq,ref_50X[start-1],length)
                                target_reads.append(aligned_read)
                                rev_dic[read] = [ref_50X[start-1],aligned_read]
				print '>read_50X-'+str(i)
				print replace_all_in_string(aligned_read,'-','N')
                elif 'DLA-12' in allele:
                        if (start-1 <= range_D12[0] and start-2+len(seq) >= start_D12) or (start-1 <= end_D12 and start-2+len(seq) >= range_D12[1]):
                                aligned_read = pseudo_align(seq,ref_D12[start-1],length)
                                target_reads.append(aligned_read)
                                rev_dic[read] = [ref_D12[start-1],aligned_read]
				print '>read_D12-'+str(i)
                                print replace_all_in_string(aligned_read,'-','N')
                else:
                        if (start-1 <= range_norm[0] and start-2+len(seq) >= start_norm) or (start-1 <= end_norm and start-2+len(seq) >= range_norm[1]):
                                aligned_read = pseudo_align(seq,ref_norm[start-1],length)
                                target_reads.append(aligned_read)
                                rev_dic[read] = [ref_norm[start-1],aligned_read]
				print '>read_norm-'+str(i)
                                print replace_all_in_string(aligned_read,'-','N')
	consensus0 = '' #useless 
	polymorphic_site = {} #polymorphic sites from gap regions
	if len(target_reads) > 1:
		consensus0,polymorphic_site = consensus_polymorphic(target_reads)
	forwd_reads = forwd_dic.keys()
	rev_reads = rev_dic.keys()
	shared = list(set(forwd_reads).intersection(rev_reads))
	# pattern linkages
	result = {}
	poly_sites = polymorphic_site.keys()
	new_linkage_patterns = {}
	for i in range(len(patterns)):
		pattern = patterns[i].split('\n')[0]
		consensus = ''.join(patterns[i].split('\n')[1:])
		old_site = re.findall(r'\d+',pattern)
		nuc = re.findall(r'[A-Za-z]',pattern)
		new_site = {}
		#sites = new_site.keys()
		for j in range(len(old_site)):
			new = ref_norm[dic_norm[int(old_site[j])]]
			new_site[new] = nuc[j]
		overall_sites = poly_sites + new_site.keys()
		aligned_consensus = pseudo_align(consensus,ref_norm[dic_norm[0]],length)
		pattern_reads = [aligned_consensus]
		linkages = []
		new_sites_linkage_reads = {}
		flag = True #flag if frequency based extension are necessary. False if read pairs spaning old and new polymorphic sites are found
		for j in range(len(shared)):
			read = shared[j]
			forwd_start = int(forwd_dic[read][0])
			forwd_seq = forwd_dic[read][1]
			raw_forwd_seq_len = len(remove_tail_char(forwd_seq,'-'))
			rev_start = int(rev_dic[read][0])
			rev_seq = rev_dic[read][1]
			raw_rev_seq_len = len(remove_tail_char(rev_seq,'-'))
			if (forwd_seq in aligned_consensus and forwd_start <= any(new_site.keys()) <= forwd_start + raw_forwd_seq_len and rev_start <= any(polymorphic_site.keys()) <= rev_start + raw_rev_seq_len) or (rev_seq in aligned_consensus and rev_start <= any(new_site.keys()) <= rev_start + raw_rev_seq_len and forwd_start <= any(polymorphic_site.keys()) <= forwd_start + raw_forwd_seq_len):
				pattern_reads.append(forwd_seq)
				pattern_reads.append(rev_seq)
				flag = False
			site_nucs = [] # Gap region polymorphic sites. Reads
			for k in range(len(overall_sites)):
				if forwd_seq[int(overall_sites[k])] != '-':
					info = str(overall_sites[k]) + forwd_seq[int(overall_sites[k])]
					site_nucs.append(info)
				if rev_seq[int(overall_sites[k])] != '-':
					info = str(overall_sites[k]) + rev_seq[int(overall_sites[k])]
					site_nucs.append(info)
			site_nucs = list(set(site_nucs))
			site_nucs.sort(key=natural_keys,reverse = True)
			if len(site_nucs) >= 2:
				linkages.append(site_nucs)
				l = '\t'.join(site_nucs)
				if l not in new_sites_linkage_reads.keys():
					new_sites_linkage_reads[l] = []
				new_sites_linkage_reads[l].append(forwd_seq)
				new_sites_linkage_reads[l].append(rev_seq)
		for j in range(len(forwd_reads)):
			if forwd_reads[j] not in shared:
				site_nucs = []
				forwd_seq = forwd_dic[forwd_reads[j]][1]
				for k in range(len(overall_sites)):
                                	if forwd_seq[int(overall_sites[k])] != '-':
                                        	info = str(overall_sites[k]) + forwd_seq[int(overall_sites[k])]
                                        	site_nucs.append(info)
                        	site_nucs.sort(key=natural_keys,reverse = True)
                        	if len(site_nucs) >= 2:
                                	linkages.append(site_nucs)
					l = '\t'.join(site_nucs)
					if site_nucs not in new_sites_linkage_reads.keys():
						new_sites_linkage_reads[l] = []
					new_sites_linkage_reads[l].append(forwd_seq)
		for j in range(len(rev_reads)):
			if rev_reads[j] not in shared:
				site_nucs = []
				rev_seq = rev_dic[rev_reads[j]][1]
				for k in range(len(overall_sites)):
					if rev_seq[int(overall_sites[k])] != '-':
						info = str(overall_sites[k]) + rev_seq[int(overall_sites[k])]
						site_nucs.append(info)
				site_nucs.sort(key=natural_keys,reverse = True)
				if len(site_nucs) >= 2:
                                	linkages.append(site_nucs)
					l = '\t'.join(site_nucs)
					if site_nucs not in new_sites_linkage_reads.keys():
						new_sites_linkage_reads[l] = []
					new_sites_linkage_reads[l].append(rev_seq)
		output,Freq = Merge_linkages_and_frequency(linkages)
		new_pattern = []
                keys = new_site.keys()
                for j in range(len(keys)):
                        new_pattern.append(str(keys[j])+new_site[keys[j]])
                new_pattern.sort()
                new_pattern = pattern + '\t(' + '\t'.join(new_pattern) + ')'
		new_linkage_patterns[new_pattern] = []
		new_sites_frequency = []
		for i in range(len(output)):
			output[i].sort(key=natural_keys,reverse = True)
			string = str(Freq[i]) + '\t' + '\t'.join(output[i]) + '\n'
			new_sites_frequency.append(string)
			new_linkage_patterns[new_pattern].append(string)
		new_sites_frequency.sort(key=natural_keys,reverse = True)
		if flag and pattern != 'Ungrouped':
			idx = patSorted.index(pattern)
			if idx >= len(new_sites_frequency):
				idx = 0
			if len(new_sites_frequency) != 0:
				p = '\t'.join(new_sites_frequency[idx].split('\n')[0].split('\t')[1:])
				reads = new_sites_linkage_reads[p]
				print pattern
				print p
				for read in reads:
					pattern_reads.append(read)
					print read
		pattern_consensus,polymorphic_site2 = consensus_polymorphic(pattern_reads)
		pattern_consensus = remove_tail_char(pattern_consensus,'N')
		result[pattern] = pattern_consensus
	return result,polymorphic_site,new_linkage_patterns

# merge linkages and summerize the frequency of each pattern
def Merge_linkages_and_frequency(lst):
        # generate the frequency list
        freq = []
        for i in range(len(lst)):
                freq.append(1)
        # pairwise merging
        for i in range(len(lst)-1):
                for j in range(i+1,len(lst)):
                        if lst[i] != 'NA' and lst[j] != 'NA':
                                if all(x in lst[i] for x in lst[j]): # j is a subset of i
                                        lst[j] = 'NA'
                                        freq[i] += freq[j]
                                        freq[j] = 0
                                elif all(x in lst[j] for x in lst[i]): # i is a subset of j
                                        lst[i] = 'NA'
                                        freq[j] += freq[i]
                                        freq[i] = 0
        # outputs
        Freq = []
        output = []
        for i in range(len(lst)):
                if lst[i] != 'NA':
                        output.append(lst[i])
                        Freq.append(freq[i])
        return output,Freq

# sequence dictionary
def build_seq_dic(lst):
        dic = {}
        for i in range(len(lst)):
                header = lst[i].split('\t')[0]
                seq = lst[i].split('\t')[9]
                dic[header] = seq
        return dic

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

# remove all occurance of a character from a list
def remove_values_from_list(the_list, val):
        while val in the_list:
                the_list.remove(val)

# consensus sequence extension using paired end reads
def consensus_extension(consensus_dic,forwd_dic,rev_dic,cutoff):
        patterns = consensus_dic.keys()
        result = {}
        for i in range(len(patterns)):
                consensus = consensus_dic[patterns[i]]
                reads = list(set(forwd_dic.keys()).intersection(rev_dic.keys()))
                length = -1
                while length != len(reads):
                        length = len(reads)
                        for j in range(len(reads)):
                                forwd = forwd_dic[reads[j]]
                                rev = rev_dic[reads[j]]
                                forwd_ext = Two_seq_expension(consensus,forwd,cutoff)
                                rev_ext = Two_seq_expension(consensus,rev,cutoff)
                                if forwd_ext != 'NA' and rev in consensus:
                                        consensus = forwd_ext
                                        reads[j] = 'NA'
                                elif rev_ext != 'NA' and forwd in consensus:
                                        consensus = rev_ext
                                        reads[j] = 'NA'
                        remove_values_from_list(reads,'NA')
                result[patterns[i]] = consensus
        return result

##########################################################
if __name__ == '__main__':
	forwd = sys.argv[1]
	rev = sys.argv[2]
	alignment = sys.argv[3]
	patterns0 = sys.argv[4]
	patternFreq = sys.argv[5]
	cutoff = int(sys.argv[6])
	with open(forwd,'r') as f:
        	forwd = f.read()
	with open(rev,'r') as f:
        	rev = f.read()
	with open(alignment,'r') as f:
        	alignment = f.read()
	with open(patterns0,'r') as f:
        	patterns = f.read()
	with open(patternFreq,'r') as f:
        	patternFreq = f.read()
	forwd = forwd.split('\n')[:-1]
	rev = rev.split('\n')[:-1]
	alignment = alignment.split('>')[1:]
	patterns = patterns.split('>')[1:]
	patternFreq = patternFreq.split('\n')[:-1]
	# output files
	out = open(patterns0 + '_extendedWithPairedReads','w')
	out2 = open(patterns0 + '_newPolymorphicSites_Linkage','w')
	out3 = open(patterns0 + '_ConsensusAndPattern','w')
	# data processing
	dic_norm,dic_50X,dic_D12 = Reformat_pairing(alignment)
	#Exon 2 and 3 in original sequences
	start_norm = 73
	start_50X = 73
	start_D12 = 64
	end_norm = 618 # [73,618+1]
	end_50X = 621
	end_D12 = 609
	# NEW polymorphic sites index are based on reference clustal O alignment !!!!!
	consensus_dic,polymorphic_sites,new_linkage_patterns = reads_polimorphicSites(forwd,rev,alignment,patterns,dic_norm,dic_50X,dic_D12,start_norm,start_50X,start_D12,end_norm,end_50X,end_D12,patternFreq)
	forwd_dic = build_seq_dic(forwd)
	rev_dic = build_seq_dic(rev)
	result = consensus_extension(consensus_dic,forwd_dic,rev_dic,cutoff)
	keys = result.keys()
	for i in range(len(keys)):
		string = '>Consensus' + str(i+1) + '\n' + result[keys[i]] + '\n'
		out.write(string)
		string2 = '>Consensus' + str(i+1) + '\n' + keys[i] + '\n'
		out3.write(string2)
	out.close()
	out3.close()
	sites = polymorphic_sites.keys()
	for site in sites:
		nuc = '\t'.join(polymorphic_sites[site])
		string = str(site) + '\t' + nuc
		print string
	pt = new_linkage_patterns.keys()
	for i in range(len(pt)):
		string = '>' + pt[i] + '\n' +  ''.join(new_linkage_patterns[pt[i]])
		out2.write(string)
	out2.close()

