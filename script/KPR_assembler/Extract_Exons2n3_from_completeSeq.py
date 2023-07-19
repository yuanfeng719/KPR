#!/usr/bin/python

import sys

if __name__ == '__main__':
        file0 = sys.argv[1]
        outputPath = sys.argv[2]
        with open(file0,'r') as f:
                file = f.read()
        lst = file.split('>')[1:]

        out = open(outputPath + '/' + file0.split('/')[-1].split('.')[0] + '_E23.fa','w')
        for i in range(len(lst)):
                header = lst[i].split('\n')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                if 'DLA-88*5' in header:
                        length = 549
                else:
                        length = 546
                if 'GCTCCCACT' in seq:
                        start = seq.index('GCTCCCACT')
                elif 'DLA-12' in header:
                        start = 64
                else:
                        start = 73
                e23 = seq[start:start+length]
                string = '>' + header + '\n' + e23 + '\n'
                out.write(string)
        out.close()


