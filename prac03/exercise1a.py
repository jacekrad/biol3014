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

prot_trn = prot[0::2]  # even-numbered indices
prot_tst = prot[1::2]  # odd-numbered indices
sstr_trn = sstr[0::2]  # even-numbered indices
sstr_tst = sstr[1::2]  # odd-numbered indices

print len(prot_trn[0])

W = 5
print prot_trn[0]
print slidewin(prot_trn[0], W)
print len(slidewin(prot_trn[0], W))
nb = NaiveBayes([Protein_Alphabet for _ in range(W)], DSSP3_Alphabet)

# iterate over the sequence of proteins in the training set
#for i in range(len(prot_trn)):  
print "-------------------------------------------------------"
for i in range(1):  
    subseqs = slidewin(prot_trn[i], W)  # construct sub-seqs
    subtarg = sstr_trn[i][W / 2:-W / 2 + 1]  # secondary structure elem.  remove the overhang
    print subtarg
    print sstr_trn[i]
    #print subseqs, subtarg
    
    for j in range(len(subseqs)):
        nb.observe(subseqs[j], subtarg[j])  # let NB count
        print "observing", subseqs[j], subtarg[j]
        
# create an array of zeroes of length of the secondary structure alphabet
cm = numpy.zeros((len(DSSP3_Alphabet), len(DSSP3_Alphabet)))

# iterate over the proteins in the test set
for i in range(len(prot_tst)):
    subseqs = slidewin(prot_tst[i], W)
    subtarg = sstr_tst[i][W/2:-W/2+1]
    for j in range(len(subseqs)):
        out = nb[subseqs[j]]
        c_targ = DSSP3_Alphabet.index(subtarg[j])  
        c_pred = out.prob().index(max(out.prob()))
        cm[c_targ, c_pred] += 1
        
        