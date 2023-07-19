#!/usr/bin/python

import sys
import re

########################################
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

# identify the nucleotide in given poltmorphic site
def extract_polymorphic_residue(rec,dic_norm,dic_50X,dic_D12,sites):
	read = rec.split('\t')[0]
	mapped_allele = rec.split('\t')[2]
	mapped_pos = int(rec.split('\t')[3])
	seq = rec.split('\t')[9]
	# get reference index range and check if the read covers polymorphic sites
	length = len(seq)
	start = mapped_pos - 1
	end = mapped_pos + length - 2
	site_nucs = []
	if 'DLA-88*5' in mapped_allele:
		for pos in sites:
			ref_pos = dic_50X[pos]
			if ref_pos >= start and ref_pos <= end:
				nuc = seq[ref_pos-mapped_pos+1] # real position = ref_pos - mapped_pos + 1
				info = str(pos) + nuc
				site_nucs.append(info)
	elif 'DLA-12' in mapped_allele:
                for pos in sites:
                        ref_pos = dic_D12[pos]
			if ref_pos >= start and ref_pos <= end:
                                nuc = seq[ref_pos-mapped_pos+1] # real position = ref_pos - mapped_pos + 1
                                info = str(pos) + nuc
                                site_nucs.append(info)
	else:
		for pos in sites:
			ref_pos = dic_norm[pos]
			if ref_pos >= start and ref_pos <= end:
                                nuc = seq[ref_pos-mapped_pos+1] # real position = ref_pos - mapped_pos + 1
                                info = str(pos) + nuc
                                site_nucs.append(info)
	if mapped_pos <= 540 and mapped_pos + length >= 538:
		if 'CTACAATGG' in seq:	
			tag = '50X'
                elif 'CAGTACTAC' in seq:
                        tag = 'DLA64'
                else:
			tag = 'norm'
	else:
		tag = 'NA'
	site_nucs.sort()
	return (site_nucs,tag)

# analyze forwd and rev sam files 
def analyze_samFiles(forwd,rev,dic_norm,dic_50X,dic_D12,sites):
	# get pair and single reads info
	forwd_reads = []
	rev_reads = []
	rev_recs = {}
	X50_seq = []
	norm_seq = []
        DLA64_seq = []
	for i in range(len(forwd)):
		read = forwd[i].split('\t')[0]
		forwd_reads.append(read)
	for i in range(len(rev)):
		read = rev[i].split('\t')[0]
		rev_reads.append(read)
		rev_recs[read] = rev[i]
	shared_reads = list(set(forwd_reads).intersection(rev_reads))
	result = []
	### forwd
	for i in range(len(forwd)):
		rec = []
		rec2 = []
		rec = forwd[i]
		read = rec.split('\t')[0]
		site_nucs = []
		site_nucs,tagF = extract_polymorphic_residue(rec,dic_norm,dic_50X,dic_D12,sites)
		if read in shared_reads:
			rec2 = rev_recs[read]
			site_nucs2 = []
			site_nucs2,tagR = extract_polymorphic_residue(rec2,dic_norm,dic_50X,dic_D12,sites)
			site_nucs = site_nucs + site_nucs2
			site_nucs = list(set(site_nucs))
			site_nucs.sort()
		if len(site_nucs) >= 1:
			result.append(site_nucs)
		if len(site_nucs) >= 2:
			#result.append(site_nucs)
			if tagF == '50X':
				X50_seq.append(rec.split('\t')[9])
			elif tagF == 'norm':
				norm_seq.append(rec.split('\t')[9])
                        elif tagF == 'DLA64':
                                DLA64_seq.append(rec.split('\t')[9])
			if tagR == '50X':
                                X50_seq.append(rec2.split('\t')[9])
			elif tagR == 'norm':
                                norm_seq.append(rec2.split('\t')[9])
                        elif tagR == 'DLA64':
                                DLA64_seq.append(rec2.split('\t')[9])
	### rev 
	for i in range(len(rev)):
		rec = rev[i]
		read = rec.split('\t')[0]
		if read not in shared_reads:
			site_nucs = []
			site_nucs,tag = extract_polymorphic_residue(rec,dic_norm,dic_50X,dic_D12,sites)
			site_nucs.sort()
			if len(site_nucs) >= 1:
				result.append(site_nucs)
			if len(site_nucs) >= 2:
				#result.append(site_nucs)
				if tag == '50X':
					X50_seq.append(rec.split('\t')[9])
				elif tag == 'norm':
					norm_seq.append(rec.split('\t')[9])
                                elif tag == 'DLA64':
                                        DLA64_seq.append(rec.split('\t')[9])
	# calculate the average coverage for each polymorphic site (with read pairs spaning at least 2 variants)
	total_variants = 0
	for i in range(len(result)):
		rec = result[i]
		num_of_variants = len(rec)
		total_variants += num_of_variants
	if len(sites) > 0:
		average_cov = float(total_variants) / len(sites)
	else:
		average_cov = float(total_variants)
	return (result,average_cov,X50_seq,norm_seq,DLA64_seq)

# compare two patterns and decide if one contains another or can be merged
def merge_patterns(lst,lst2):
	dic = {}
	dic2 = {}
	for i in range(len(lst)):
		location = re.findall(r'\d+',lst[i])[0]
		nuc = re.findall(r'[A-Za-z]',lst[i])[0]
		if location not in dic.keys():
			dic[location] = nuc
	for i in range(len(lst2)):
                location = re.findall(r'\d+',lst2[i])[0]
                nuc = re.findall(r'[A-Za-z]',lst2[i])[0]
                if location not in dic2.keys():
                        dic2[location] = nuc
	shared = list(set(dic.keys()).intersection(dic2.keys()))
	key1 = dic.keys()
	key2 = dic2.keys()
	flag = True
	if len(shared) < len(lst)*0.75 and len(shared) < len(lst2)*0.75:
		flag = False
	for i in range(len(shared)):
		location = shared[i]
		if dic[location] != dic2[location]:
			flag = False
	result = []
	if flag:
		for i in range(len(shared)):
			key = shared[i]
			string = str(key) + dic[key]
			result.append(string)
		for i in range(len(key1)):
			key = key1[i]
			if key not in shared:
				string = str(key) + dic[key]
				result.append(string)
		for i in range(len(key2)):
			key = key2[i]
			if key not in shared:
				string = str(key) + dic2[key]
				result.append(string)
	result.sort()
	return (flag,result)

# merge linkages and summerize the frequency of each pattern
def Merge_linkages_and_frequency(lst):
	# generate the frequency list
	freq = []
	for i in range(len(lst)):
		freq.append(1)
	# pairwise merging
	NA_counter = 0
	previous_NA_count = -1
	while NA_counter != previous_NA_count:
		previous_NA_count = NA_counter
		NA_counter = 0
		for i in range(len(lst)-1):
			for j in range(i+1,len(lst)):
				if lst[i] != 'NA' and lst[j] != 'NA':
					flag,new_lst = merge_patterns(lst[i],lst[j])
					if flag:
						lst[i] = new_lst
						lst[j] = 'NA'
						freq[i] += freq[j]
						freq[j] = 0
						NA_counter += 1
	# outputs
	Freq = []
	output = []
	for i in range(len(lst)):
		if lst[i] != 'NA':
			output.append(lst[i])
			Freq.append(freq[i])
	return output,Freq

# remove all occurance of a character from a list
def remove_values_from_list(the_list, val):
        while val in the_list:
                the_list.remove(val)

# validate original contigs and linkage groups with paired-end reads (polymorphic site based).
# analyze the linkage of every two sites in all contigs and llinkage groups. [1] if reads found supporting the linkage -- PASS; [2] if no reads found, but the number number of reads spaning the region is smaller than the average coverage -- PASS (no enough evidence to reject); [3] if no reads found, and spanning read pairs greater than average -- FAIL (the linkage in contig/linkageGroup is likely caused by error)
def read_pairwise_linkage_counts(lst):
	result = {} # 26A_27T:5
	pos_count = {} # 26_27:10
	# count all pairwise linkages in paired-end raads
	for i in range(len(lst)):
		read_linkage_lst = lst[i]
		if len(read_linkage_lst) >= 2:
			for j in range(len(read_linkage_lst)-1):
				site1 = read_linkage_lst[j]
				loc1 = int(re.findall(r'\d+',site1)[0])
				for k in range(j+1,len(read_linkage_lst)):
					site2 = read_linkage_lst[k]
					loc2 = int(re.findall(r'\d+',site2)[0])
					key = ''
					key2 = ''
					if loc1 < loc2:
						key = site1 + '_' + site2
						key2 = str(loc1) + '_' + str(loc2)
					elif loc1 > loc2:
						key = site2 + '_' + site1
						key2 = str(loc2) + '_' + str(loc1)
					if key not in result.keys():
						result[key] = 1
					else:
						result[key] += 1
					if key2 not in pos_count.keys():
						pos_count[key2] = 1
					else:
						pos_count[key2] += 1
	return (result,pos_count)

# validate each contig/linkageGroup linkage and adjust the contig grouping result
def validate_contig_linkageGroup(read_linkage_dic,contig_linkages,ContigGroupsAndOutliers,cutoff,pos_count):
	problem_contig = []
	remaining_contig = {}
        problem_evidence = {}
	# analyze each contig
	for i in range(len(contig_linkages)):
		header = contig_linkages[i].split('\n')[0]
		site_lst = contig_linkages[i].split('\n')[1].split('\t')
		remove_values_from_list(site_lst,'')
		remaining_contig[header] = contig_linkages[i].split('\n')[1]
		if len(site_lst) >= 2:
			for j in range(len(site_lst)-1):
				site1 = site_lst[j]
				loc1 = int(re.findall(r'\d+',site1)[0])
				for k in range(j+1,len(site_lst)):
					site2 = site_lst[k]
					loc2 = int(re.findall(r'\d+',site2)[0])
					key = ''
					key2 = ''
					if loc1 < loc2:
                                                key = site1 + '_' + site2
						key2 = str(loc1) + '_' + str(loc2)
                                        elif loc1 > loc2:
                                                key = site2 + '_' + site1
						key2 = str(loc2) + '_' + str(loc1)
					if (key not in read_linkage_dic.keys()) and (key2 in pos_count.keys()):
						# check coverage of the same positions
						count = pos_count[key2]
						if count >= cutoff:
							problem_contig.append(header)
                                                        if header not in problem_evidence.keys():
                                                                problem_evidence[header] = []
                                                        problem_evidence[header].append(key + '(0/' + str(count) + ')')
	problem_contig = list(set(problem_contig))
	for i in range(len(problem_contig)):
		contig = problem_contig[i]
		del remaining_contig[contig]	
	# analyze linkage pattern and contig file
	result = []
	for i in range(len(ContigGroupsAndOutliers)):
		linkage_lst = ContigGroupsAndOutliers[i].split('\n')[0].split('\t')
		contigs_lst = ContigGroupsAndOutliers[i].split('\n')[1].split(';')
		flag = True
		if len(linkage_lst) >= 2:
			for j in range(len(linkage_lst)-1):
                                site1 = linkage_lst[j]
                                loc1 = int(re.findall(r'\d+',site1)[0])
                                for k in range(j+1,len(linkage_lst)):
                                        site2 = linkage_lst[k]
                                        loc2 = int(re.findall(r'\d+',site2)[0])
                                        key = ''
                                        key2 = ''
                                        if loc1 < loc2:
                                                key = site1 + '_' + site2
                                                key2 = str(loc1) + '_' + str(loc2)
                                        elif loc1 > loc2:
                                                key = site2 + '_' + site1
                                                key2 = str(loc2) + '_' + str(loc1)
                                        if (key not in read_linkage_dic.keys()) and (key2 in pos_count.keys()):
                                                # check coverage of the same positions
                                                count = pos_count[key2]
                                                if count >= cutoff:
							flag = False
		if flag: #linkage group is right
			# remove all problem contigs
			for j in range(len(problem_contig)):
				contig = problem_contig[j]
				remove_values_from_list(contigs_lst,contig)
			string = '\t'.join(linkage_lst) + '\n' + ';'.join(contigs_lst) + '\n'
			result.append(string)
	return(result,problem_contig,remaining_contig,problem_evidence)

# summerize the amount of ATGC in each polymorphic site
def summerize_site_nuc_frequency(lst):
	dic = {} # polymorphic sites and counts dictionary
	for i in range(len(lst)):
		rec = lst[i]
		for j in range(len(rec)):
			site = rec[j]
			pos = int(re.findall(r'\d+',site)[0])
			nuc = re.findall(r'\D+',site)[0]
			if pos not in dic.keys():
				dic[pos] = [0,0,0,0] # A,C,G,T
			if nuc == 'A':
				dic[pos][0] += 1
			elif nuc == 'C':
				dic[pos][1] += 1
			elif nuc == 'G':
				dic[pos][2] += 1
			elif nuc == 'T':
				dic[pos][3] += 1
	return (dic)	

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	return [ atoi(c) for c in re.split(r'(\d+)', text) ]


##########################################################
if __name__ == '__main__':
	forwd = sys.argv[1]
	rev = sys.argv[2]
	clustalO = sys.argv[3]
	poli_sites = sys.argv[4]
	contig_linkages = sys.argv[5]
	ContigGroupsAndOutliers = sys.argv[6]
	cutoff = int(sys.argv[7]) # frequency cutoff, small than which a problem pair linkage will not be removed
	with open(forwd,'r') as f:
        	forwd = f.read()
	with open(rev,'r') as f:
        	rev = f.read()
	with open(clustalO,'r') as f:
        	alignment = f.read()
	with open(poli_sites,'r') as f:
        	poli_sites = f.read()
	with open(contig_linkages,'r') as f:
        	contig_linkages = f.read()
	with open(ContigGroupsAndOutliers,'r') as f:
        	ContigGroupsAndOutliers = f.read()
	alignment = alignment.split('>')[1:]
	poli_sites = poli_sites.split('\n')[:-1]
	forwd = forwd.split('\n')[:-1]
	rev = rev.split('\n')[:-1]
	contig_linkages = contig_linkages.split('>')[1:]
	ContigGroupsAndOutliers = ContigGroupsAndOutliers.split('>')[1:]
	# output files
	out = open(clustalO+'_PolymorphicSites_Linkage','w')
	out2 = open(clustalO+'_correctedContigGroupsAndOutliers','w')
	out3 = open(clustalO+'_correctedContigLinkages','w')
	# data processing
	dic_norm,dic_50X,dic_D12 = Reformat_pairing(alignment)
	sites = []
	for i in range(len(poli_sites)):
		sites.append(int(poli_sites[i].split('\t')[0]))
	result,average_cov,X50_seq,norm_seq,DLA64_seq = analyze_samFiles(forwd,rev,dic_norm,dic_50X,dic_D12,sites)
	outT = open('Script2_readPairs_polymorphicSites_linkages','w')
	for i in range(len(result)):
		rec = result[i]
		rec.sort(key=natural_keys)
		string = '\t'.join(rec) + '\n'
		outT.write(string)
	outT.close()
	outT = open('Script2_50X_sequences','w')
	for i in range(len(X50_seq)):
		string = '>50X_' + str(i+1) + '\n' + X50_seq[i] + '\n'
		outT.write(string)
	outT.close()
	outT = open('Script2_norm_sequences','w')
	for i in range(len(norm_seq)):
        	string = '>norm_' + str(i+1) + '\n' + norm_seq[i] + '\n'
        	outT.write(string)
	outT.close()
        outT = open('Script2_DLA64_sequences','w')
        for i in range(len(DLA64_seq)):
                string = '>DLA64_' + str(i+1) + '\n' + DLA64_seq[i] + '\n'
                outT.write(string)
        outT.close()
	nuc_freq_dic = summerize_site_nuc_frequency(result)
	positions = nuc_freq_dic.keys()
	positions.sort()
	outT = open('Script2_polymorphicSites_nucleotide_frequencies','w')
	header = '\tA\tC\tG\tT\n'
	outT.write(header)
	for i in range(len(positions)):
		pos = positions[i]
		rec = nuc_freq_dic[pos]
		string = str(pos) + '\t' + str(rec[0]) + '\t' + str(rec[1]) + '\t' + str(rec[2]) + '\t' + str(rec[3]) + '\n'
		outT.write(string)
	outT.close()
	read_pairwise_linkageSite_counts,pos_count = read_pairwise_linkage_counts(result)
	outT = open('Script2_pairwise_linkage','w')
	keys = read_pairwise_linkageSite_counts.keys()
	for i in range(len(keys)):
		key = keys[i]
		string = key + '\t' + str(read_pairwise_linkageSite_counts[key]) + '\n'
		outT.write(string)
	outT.close()
	outT = open('Script2_pairwise_position_counts','w')
	keys = pos_count.keys()
	for i in range(len(keys)):
        	key = keys[i]
        	string = key + '\t' + str(pos_count[key]) + '\n'
        	outT.write(string)
	outT.close()
	new_ContigGroupsAndOutliers,problem_contig,remaining_contig,problem_evidence = validate_contig_linkageGroup(read_pairwise_linkageSite_counts,contig_linkages,ContigGroupsAndOutliers,cutoff,pos_count)

        outT = open('Script2_problem_contig','w')
        for contig in problem_contig:
                evidence = problem_evidence[contig]
                string = contig + '\t' + ','.join(evidence) + '\n'
	        outT.write(string)
	outT.close()
	for i in range(len(new_ContigGroupsAndOutliers)):
		out2.write('>' + new_ContigGroupsAndOutliers[i])
	out2.close()
	remain_cont = remaining_contig.keys()
	for i in range(len(remain_cont)):
		contig = remain_cont[i]
		string = '>' + contig + '\n' + remaining_contig[contig] + '\n'
		out3.write(string)
	linkage_freq_dic = {}
	for i in range(len(new_ContigGroupsAndOutliers)):
		linkage_pattern  = new_ContigGroupsAndOutliers[i].split('\n')[0]
		if linkage_pattern != 'Ungrouped':
			linkage_freq_dic[linkage_pattern] = 0
			linkage = linkage_pattern.split('\t')
			for j in range(len(result)):
				if set(result[j]).issubset(set(linkage)):
					linkage_freq_dic[linkage_pattern] += 1
	outT = open('Script2_PolymorphicSites_Linkage','w')
	linkage_patterns = linkage_freq_dic.keys()
	for i in range(len(linkage_patterns)):
		linkage = linkage_patterns[i]
		linkage_lst = linkage.split('\t')
		linkage_lst.sort(key=natural_keys)
		linkage2 = '\t'.join(linkage_lst)
		string = str(linkage_freq_dic[linkage]) + '\t' + linkage + '\n'
		string2 = str(linkage_freq_dic[linkage]) + '\t' + linkage2 + '\n'
		out.write(string)
		outT.write(string2)
	out.close()
	outT.close()

