'''
Created on 19/08/2014

@author: jacekrad
'''
from sequence import *
from prob import *
from __builtin__ import exit
from sys import stderr
from spred import slidewin

prot = readFastaFile('prot2.fa', Protein_Alphabet)
sstr = readFastaFile('sstr3.fa', DSSP3_Alphabet)

if len(prot) != len(sstr):
    stderr.write("prot2.fa contains different number of sequences than sstr3.fa\n")
    exit(1)
    
# index
i = 0    
for protein in prot:
    secondary_structure = sstr[i]
    i += 1
    # make sure the name is the same for both protein
    # and its secondary structure entry
    if protein.name != secondary_structure.name:
        stderr.write(protein.name + " and " + secondary_structure.name \
                     + " are not the same where they should be\n")
    # check each amino acid and the corresponding secondary structure
    if len(protein.sequence) != len(secondary_structure.sequence):
        stderr.write(protein.name + " and it's secondary structure sequence " 
                     + " differ in length: " + str(len(protein.sequence)) \
                     + " vs " + str(len(secondary_structure.sequence)) + "\n")
    # checking for invalid characters in the secondary structure as specified in
    # the question, but this is not strictly required as the check is already done
    # by readFastaFile (check against specified alphabet)
    elif 'C' in secondary_structure or 'H' in secondary_structure or 'H' in secondary_structure:
        pass
    else:
        print "invalid character found in secondary structure of " + secondary_structure.name 
            
            
print "checked", i, "entries"

prot_trn = prot[0::2] # even-numbered indices
prot_tst = prot[1::2] # odd-numbered indices
sstr_trn = sstr[0::2] # even-numbered indices
sstr_tst = sstr[1::2] # odd-numbered indices

W = 5
slidewin(prot_trn[0], W)
nb = NaiveBayes([Protein_Alphabet for _ in range(W)], DSSP3_Alphabet)

for i in range(len(prot_trn)): # for each sequence
    subseqs = slidewin(prot_trn[i], W) # construct sub-seqs
    subtarg = sstr_trn[i][W/2:-W/2+1] # secondary structure elem
    for j in range(len(subseqs)):
        nb.observe(subseqs[j], subtarg[j]) # let NB count
