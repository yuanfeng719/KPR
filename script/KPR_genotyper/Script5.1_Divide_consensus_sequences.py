#!/usr/bin/python

import sys

#################
# remvoe X from two end fo protein sequences
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


####################
if __name__ == '__main__':
	protein = sys.argv[1]
	result_path = sys.argv[2]
	with open(protein,'r') as f:
        	file=f.read()
	lst = file.split('>')[1:]
	if result_path[-1] != '/':
		result_path += '/'
	for i in range(len(lst)):
		header = lst[i].split('\n')[0]
		seq = ''.join(lst[i].split('\n')[1:])
		new_seq = remove_tail_char(seq,'X')
		out = open(result_path + header,'w')
		string = '>' + header + '\n' + new_seq + '\n'
		out.write(string)
		out.close()



