#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 03/08/2014

@author: jacekrad
'''
from webservice import getGenes

yeast_file = open("yeast_transcriptome.txt")
protein_identifiers = []  # list of protein IDs

# total number of proteins read from file
protein_count = 0

for line in yeast_file:
    protein_count += 1
    fields = line.strip().split("\t")
    protein_id = fields[0]
    protein_identifiers.append(protein_id)
yeast_file.close()

# get both GO terms so we only need to make single call to UniProt
genes = getGenes(["GO:0005634", "GO:0003700"], "UniProtKB", "559292")

all_in_nucleus_proteins = genes.get("GO:0005634")
all_tf_proteins = genes.get("GO:0003700")

# filter for entries which are both in the result genes from UniProtKB
# as well as in out file (protein_identifiers)
yeast_in_nucleus_proteins = [prot_id for prot_id in all_in_nucleus_proteins \
                             if prot_id in protein_identifiers]

yeast_tf_proteins = [prot_id for prot_id in all_tf_proteins if prot_id \
                     in protein_identifiers]

yeast_in_nucleus_tf_proteins = [prot_id for prot_id in protein_identifiers \
                                if prot_id in all_in_nucleus_proteins \
                                and prot_id in all_tf_proteins]

# get number of entries for each list
in_nucleus_count = len(yeast_in_nucleus_proteins)
tf_count = len(yeast_tf_proteins) 
tf_in_nucleus_count = len(yeast_in_nucleus_tf_proteins)

# U - universal set (all proteins from file)
# T - set of transcription factor proteins
# N - set of in nucleus proteins
# T\N = T∩N' = T – T∩N
# N\T = N∩T' = N – T∩N
# T`∩N' = (T∪N)' = U - (T∪N) = U – T – N + T∩N
not_tf_nucleus_count = in_nucleus_count - tf_in_nucleus_count
tf_not_nucleus_count = tf_count - tf_in_nucleus_count
not_tf_not_nucleus_count = protein_count - tf_count - in_nucleus_count + tf_count
 
print "                           count probability"
print "    In nucleus and TF        :", tf_in_nucleus_count, \
    float(tf_in_nucleus_count) / protein_count
print "    Not in nucleus and TF    :", tf_not_nucleus_count, \
    float(tf_not_nucleus_count) / protein_count
print "    In nucleus and not TF    :", not_tf_nucleus_count, \
    float(not_tf_nucleus_count) / protein_count
print "    Not in nucleus and not TF:", not_tf_not_nucleus_count, \
    float(not_tf_not_nucleus_count) / protein_count

