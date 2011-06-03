#!/usr/bin/python
#!/users/jl38/import/Python-2.6/python

import sys
import re

#argv[1]: you target directory for storing the input genome, each file is a chromosome
#read sysnin file from the stand input

lines=sys.stdin.readlines()
chrom={}
for l in lines:
#    print l.split()
    (a,b)=l.split()[0:2]
    if b not in chrom:
        chrom[b]=[]
    chrom[b].append(a)
for each in chrom:
    output=sys.argv[1]+'/'+each+'.txt'
    f=open(output,'w')
    for i in chrom[each]:
        f.write(i+' ')

