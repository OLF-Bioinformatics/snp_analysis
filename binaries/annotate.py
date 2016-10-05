#!/usr/bin/env python

'''
Created on Sep 13, 2016

@author: Marc-Olivier Duceppe
'''

from Bio import SeqIO
from sys import argv

# infile arg used to make compatible for both sorted and organized tables
script, gbk_file, position_file, annotation_table = argv

#GenBank file
gbk_handle = open(gbk_file, "rU")
gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_handle, "genbank"))
gbk_handle.close()

#Position file
pos_handle = open(position_file, "rU")
with pos_handle as f:
    posList = set([int(line.strip().split('-',2)[1]) for line in f]) #duplicates removed
    posListSorted = sorted(posList)
pos_handle.close()

#Output file (annotation_table)
table_handle = open(annotation_table, 'a')

#output file header
table_handle.write('chromosome' + '\t' + 'position' + '\t' + 'locus_tag' + '\t' + 'gene' +  '\t' + 'product' + '\n')

#iterate through chromosomes
chromList = list(gbk_dict.keys())
for chrom in chromList: #for each chromosome
    feat = list(gbk_dict[chrom].features)
    del feat[0] #delete first entry, which is the whole chromosome
    
    for f in feat: #for each feature
        if f.type == 'CDS':
            low = f.location.nofuzzy_start
            high = f.location.nofuzzy_end
            for pos in posListSorted: #for each SNP -> not optimal because it reads all the SNP for each feature
                if low <= pos <= high:
                    
                    mygene = '-'
                    if 'gene' in f.qualifiers:
                        mygene = str(', '.join(f.qualifiers['gene']))
                    
                    #SNPs not in as CDS are not in the output file
                    table_handle.write(chrom.split('.')[0] + '\t' + str(pos) + '\t' + str(', '.join(f.qualifiers['locus_tag'])) + '\t' + mygene + '\t' + str(', '.join(f.qualifiers['product'])) + '\n')
        #else:
        #    table_handle.write(chrom.split('.')[0] + '\t' + str(pos) + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\n') 

table_handle.close()
