#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 02/08/2014

@author: jacekrad
'''
from webservice import getGOTerms

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

# dictionary of go_terms with the gene ID as the key
go_terms_dict = getGOTerms(protein_identifiers)

in_nucleus_proteins = filter(lambda x: "GO:0005634" in x[1], \
                             go_terms_dict.iteritems())
tf_proteins = filter(lambda x: "GO:0003700" in x[1], \
                     go_terms_dict.iteritems())
tf_in_nucleus_proteins = filter(lambda x: "GO:0003700" in x[1] \
                                and "GO:0005634" in x[1], \
                                go_terms_dict.iteritems())

in_nucleus_count = len(in_nucleus_proteins)
tf_count = len(tf_proteins) 
tf_in_nucleus_count = len(tf_in_nucleus_proteins)

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
    
