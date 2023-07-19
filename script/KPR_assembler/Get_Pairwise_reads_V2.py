#!/usr/bin/python

import sys

###########################
# build pools
def get_name(info,tag):
	lst = []
	if tag == 'fa':
		for i in range(len(info)):
			header = info[i].split('\n')[0]
			if header not in lst:
				lst.append(header)
	elif tag == 'sam':
		for i in range(len(info)):
			header = info[i].split('\t')[0]
			if header not in lst:
				lst.append(header)
	return (lst)

def intersection(lst1, lst2): 
	lst3 = [value for value in lst1 if value in lst2] 
	return lst3

def remove_values_from_list(the_list, val):
	return [value for value in the_list if value != val]

def check_pool(lst,inters,tag):
	out = []
	if tag == 'fa':
		for i in range(len(lst)):
			header = lst[i].split('\n')[0]
			if header in inters:
				inters = remove_values_from_list(inters,header)
				out.append('>'+lst[i])
	elif tag == 'sam':
		for i in range(len(lst)):
			header = lst[i].split('\t')[0]
			if header in inters:
				inters = remove_values_from_list(inters,header)
				out.append(lst[i] + '\n')
	return (out)


####################################
if __name__ == '__main__':
	forwd = sys.argv[1]
	rev = sys.argv[2]
	with open(forwd,'r') as f:
        	forward = f.read()
	with open(rev,'r') as f:
        	reverse = f.read()
	#data processing
	tag = 'sam'
	if forward[0] == '>' and reverse[0] == '>':
        	fd = forward.split('>')[1:]
        	rv = reverse.split('>')[1:]
        	tag = 'fa'
	elif forward[0] != '>' and reverse[0] != '>':
        	fd = forward.split('\n')[:-1]
        	rv = reverse.split('\n')[:-1]
	else:
        	fd = ''
        	rv = ''
	forward = get_name(fd,tag)
	reverse = get_name(rv,tag)
	overlap = intersection(forward,reverse)
	fd = check_pool(fd,overlap,tag)
	rv = check_pool(rv,overlap,tag)
	# output files
	f = open(forwd+'-pair','w')
	r = open(rev+'-pair','w')
	string = ''.join(fd)
	f.write(string)
	string = ''.join(rv)
	r.write(string)
	f.close()
	r.close()

