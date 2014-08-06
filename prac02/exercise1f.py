'''
Created on 06/08/2014

@author: s4361277
'''
from sequence import *
from gibbs import *
from gibbs import *
seqs = readFastaFile('hth_40.fa', Protein_Alphabet)
W = 24 # the width of the motif sought
g = GibbsMotif(seqs, W)
q = g.discover()
p = g.getBackground()
a = getAlignment(seqs, q, p)
k = 0
for seq in seqs:
    print "%s \t%s" % (seq.name, seq[a[k]:a[k]+W])
    k += 1

foreground = None
background = readDistrib("blosum62.distrib")

pwm = PWM(foreground, background)

ids = []
seqs = readFastaFile('HPRD.fa', Protein_Alphabet)
for seq in seqs:
    hits = pwm.search(seq)
    if len(hits) > 0:
        print "%s \t%d \t%s \t%5.3f" % (seq.name, hits[0][0], hits[0][1], hits[0][2])
        ids.append(seq.name)
report = getGOReport(ids) 
for row in report:
    print "%s \t%d \t%s" % row