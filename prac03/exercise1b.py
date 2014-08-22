'''
Created on 20/08/2014

@author: jacekrad
'''
from sequence import *
from prob import *
from __builtin__ import exit
from sys import stderr
from spred import slidewin

prot = readFastaFile('prot2.fa', Protein_Alphabet)
sstr = readFastaFile('sstr3.fa', DSSP3_Alphabet)

prot_trn = prot[0::2]  # even-numbered indices
prot_tst = prot[1::2]  # odd-numbered indices
sstr_trn = sstr[0::2]  # even-numbered indices
sstr_tst = sstr[1::2]  # odd-numbered indices

W = 5
print "Sequence length  :", len(prot_trn[0])
print "Number of windows:", len(slidewin(prot_trn[0], W))

