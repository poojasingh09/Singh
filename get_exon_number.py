#!/usr/local/bin/python



import sys
import os
# Count number of genes with one or more variants

import csv
#from collections import Counter

reflist = []
foundlist = []
uniquelist = []
genecounts = []

f1 = 'all.exons.exonCount'

f1 = open(f1, 'r')
theVARIANTs = f1.readlines()
for line in theVARIANTs:
	snpDescription = line.split('\t')[0]
	reflist.append(snpDescription)


theFiles = sys.argv[1:]
#print(theFiles[0])
for file in theFiles:
	f2 = open(file,'r')
	thelines = f2.readlines()
	f3 = open(file+'.uniq_DEU_Count','a')
	f3.write("gene" + "\t" + "DEU_Count" + "\n")
	for line in thelines:
		gene = line.split(":")[0]
		#finderStart = gene.find("on.gene.")
		#finderEnd = gene.find(":E")
		#geneID = gene[finderStart:finderEnd]
		foundlist.append(gene)
		#print geneID

	


#counts = 0
#k = 0
for name in foundlist:
	if name not in uniquelist:
		uniquelist.append(name)
		
		#k = k + 1
#Counter(foundlist)
for name in uniquelist:
	counts = 0
	for compName in foundlist:
		if compName == name:
			counts = counts + 1
	genecounts.append(counts)
	#print name, counts
#f3 = open(f3, 'a')
k = 0
for entry in uniquelist:
	f3.write(entry + "\t" +  str(genecounts[k]) + "\n")
	k = k + 1


f3.close()
