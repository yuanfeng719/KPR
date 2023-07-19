#!/usr/bin/python

import sys

def check_quality(rec,minimum):
        quality = rec.split('\t')[10]
        flag = True
        for i in range(len(quality)):
                if dic[quality[i]] < minimum:
                        flag = False
        return (flag)

if __name__ == '__main__':
	# https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
	dic = {'!':0,'"':1,'#':2,'$':3,'%':4,'&':5,"'":6,'(':7,')':8,'*':9,'+':10,',':11,'-':12,'.':13,'/':14,'0':15,'1':16,'2':17,'3':18,'4':19,'5':20,'6':21,'7':22,'8':23,'9':24,':':25,';':26,'<':27,'=':28,'>':29,'?':30,'@':31,'A':32,'B':33,'C':34,'D':35,'E':36,'F':37,'G':38,'H':39,'I':40,'J':41,'K':42}
	sam_file = sys.argv[1]
	minimum = int(sys.argv[2])
	with open(sam_file,'r') as f:
        	sam = f.read()
	lst = sam.split('\n')[:-1]
	out = open(sam_file + '_QualityFiltered','w')
	for i in range(len(lst)):
        	Flag = check_quality(lst[i],minimum)
        	if Flag == True:
                	string = lst[i] + '\n'
                	out.write(string)
	out.close()


