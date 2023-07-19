#!/usr/bin/python

import sys

###########################
# build pools
def get_name(info):
	lst = []
	for i in range(len(info)):
		header = info[i].split('\n')[0]
		if header not in lst:
			lst.append(header)
	return (lst)

def intersection(lst1, lst2): 
	lst3 = [value for value in lst1 if value in lst2] 
	return lst3

def remove_values_from_list(the_list, val):
	return [value for value in the_list if value != val]

def check_pool(lst,inters):
	out = []
	for i in range(len(lst)):
		header = lst[i].split('\n')[0]
		if header in inters:
			inters = remove_values_from_list(inters,header)
			out.append('>'+lst[i])
	return (out)


###############################
if __name__ == '__main__':
	forwd = sys.argv[1]
	rev = sys.argv[2]
	with open(forwd,'r') as f:
        	forward = f.read()
	with open(rev,'r') as f:
        	reverse = f.read()
	fd = forward.split('>')[1:]
	rv = reverse.split('>')[1:]
	# data processing
	forward = get_name(fd)
	reverse = get_name(rv)
	overlap = intersection(forward,reverse)
	fd = check_pool(fd,overlap)
	rv = check_pool(rv,overlap)
	# output files
	f = open(forwd+'-pair','w')
	r = open(rev+'-pair','w')
	string = ''.join(fd)
	f.write(string)
	string = ''.join(rv)
	r.write(string)
	f.close()
	r.close()





