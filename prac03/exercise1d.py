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

print "W, Q3, C, E, H"

for W in range(3, 50):
    nb = NaiveBayes([Protein_Alphabet for _ in range(W)], DSSP3_Alphabet)

    # iterate over the sequence of proteins in the training set
    # for i in range(len(prot_trn)):  
    for i in range(len(prot_trn)):  
        subseqs = slidewin(prot_trn[i], W)  # construct sub-seqs
        subtarg = sstr_trn[i][W / 2:-W / 2 + 1]  # secondary structure elem.  remove the overhang
        for j in range(len(subseqs)):
            nb.observe(subseqs[j], subtarg[j])  # let NB count
            # print "observing", subseqs[j], subtarg[j]
        
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
            # print out
            # print c_targ, c_pred
            cm[c_targ, c_pred] += 1        
    Q3 = Qk(cm, DSSP3_Alphabet)
    print W, ',', Q3[0], ',', Q3[1].get('C'), ',', Q3[1].get('E'), ',', Q3[1].get('H')

