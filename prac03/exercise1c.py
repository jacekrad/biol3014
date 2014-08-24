'''
Created on 20/08/2014

@author: jacekrad
'''
from sequence import *
from prob import *
from __builtin__ import exit
from sys import stderr
from spred import slidewin
from ml import Qk

prot = readFastaFile('prot2.fa', Protein_Alphabet)
sstr = readFastaFile('sstr3.fa', DSSP3_Alphabet)

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
# for i in range(len(prot_trn)):  
print "-------------------------------------------------------"
for i in range(len(prot_trn)):  
    subseqs = slidewin(prot_trn[i], W)  # construct sub-seqs
    subtarg = sstr_trn[i][W / 2:-W / 2 + 1]  # secondary structure elem.  remove the overhang
    print len(subtarg)
    print len(sstr_trn[i])
    # print subseqs, subtarg
    
    for j in range(len(subseqs)):
        nb.observe(subseqs[j], subtarg[j])  # let NB count
        #print "observing", subseqs[j], subtarg[j]
        
# create an array of zeroes of length of the secondary structure alphabet
cm = numpy.zeros((len(DSSP3_Alphabet), len(DSSP3_Alphabet)))

# iterate over the proteins in the test set
for i in range(len(prot_tst)):
    subseqs = slidewin(prot_tst[i], W)
    subtarg = sstr_tst[i][W / 2:-W / 2 + 1]
    for j in range(len(subseqs)):
        out = nb[subseqs[j]]
        c_targ = DSSP3_Alphabet.index(subtarg[j])  
        c_pred = out.prob().index(max(out.prob()))
        #print out
        #print c_targ, c_pred
        cm[c_targ, c_pred] += 1
        
print DSSP3_Alphabet        
print cm
 
 
Q = Qk(cm, DSSP3_Alphabet)

       
#z = numpy.zeros((2, 4))
#print z


