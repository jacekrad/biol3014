'''
Created on 02/08/2014

@author: jacekrad
'''
from webservice import getGOTerms

yeast_file = open("yeast_transcriptome_100.txt")
protein_identifiers = [] # list of protein IDs

# total number of proteins read from file
protein_count = 0

# counters for TF and nucleus proteins
tf_nucleus_count = 0
not_tf_nucleus_count = 0
tf_not_nucleus_count = 0
not_tf_not_nucleus_count = 0

for line in yeast_file:
    protein_count += 1
    fields = line.strip().split("\t")
    protein_id = fields[0]
    protein_identifiers.append(protein_id) # ex 5.1
yeast_file.close()

# dictionary of go_terms with the gene ID as the key
go_terms_dict = getGOTerms(protein_identifiers)

for protein_id in go_terms_dict:
    go_terms = go_terms_dict.get(protein_id)
    if "GO:0005634" in go_terms: # check if nucleus
        if "GO:0003700" in go_terms: # check if TF
            tf_nucleus_count += 1
        else:
            not_tf_nucleus_count += 1
    else:
        if "GO:0003700" in go_terms: # check if TF
            tf_not_nucleus_count += 1
        else:
            not_tf_not_nucleus_count += 1



print "Exercise 5.1 - protein identifiers: ", protein_identifiers
print "Exercise 5.2 - GO Terms"
print go_terms_dict
print "Exercise 5.3 - nucleus and TF counts and probabilities:"
print "                           count probability"
print "    In nucleus and TF        :", tf_nucleus_count, \
    float(tf_nucleus_count) / protein_count
print "    Not in nucleus and TF    :", tf_not_nucleus_count, \
    float(tf_not_nucleus_count) / protein_count
print "    In nucleus and not TF    :", not_tf_nucleus_count, \
    float(not_tf_nucleus_count) / protein_count
print "    Not in nucleus and not TF:", not_tf_not_nucleus_count, \
    float(not_tf_not_nucleus_count) / protein_count
    
print "Exercise 5.4: probability of TF and in nucleus = ", \
    float(tf_nucleus_count) / protein_count
print "            probability of in nucleus given TF = P(N|TF) = P(N,TF) / P(TF)"
print "                                               = ", \
    (float(tf_nucleus_count) / protein_count) / (float(tf_nucleus_count + tf_not_nucleus_count) / protein_count)
