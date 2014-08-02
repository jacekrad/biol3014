'''
Created on 02/08/2014

@author: jacekrad
'''
from webservice import getGOTerms

yeast_file = open("yeast_transcriptome_10.txt")
protein_identifiers = [] # list of protein IDs
go_terms_dict = {} # dictionary of go_terms with the gene ID as the key

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
    protein_identifiers.append(protein_id)
    go_terms = getGOTerms(fields[0])
    go_terms_dict[protein_id] = go_terms
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
file.close()

print "Exercise 4.1 - protein identifiers: ", protein_identifiers
print "Exercise 4.2 - GO Terms"
for protein_id in go_terms_dict:
    print protein_id, go_terms_dict[protein_id]
print "Exercise 4.3 - nucleus and TF counts and probabilities:"
print "                           count probability"
print "    In nucleus and TF        :", tf_nucleus_count, \
    float(tf_nucleus_count) / protein_count
print "    Not in nucleus and TF    :", tf_not_nucleus_count, \
    float(tf_not_nucleus_count) / protein_count
print "    In nucleus and not TF    :", not_tf_nucleus_count, \
    float(not_tf_nucleus_count) / protein_count
print "    Not in nucleus and not TF:", not_tf_not_nucleus_count, \
    float(not_tf_not_nucleus_count) / protein_count
    
print "Exercise 4.4: probability of TF and in nucleus = ", \
    float(tf_nucleus_count) / protein_count
print "            probability of in nucleus given TF = P(N|TF) = P(N,TF) / P(TF)"
print "                                               = ", \
    (float(tf_nucleus_count) / protein_count) / (float(tf_nucleus_count + tf_not_nucleus_count) / protein_count)
