'''
Created on 06/08/2014

@author: s4361277
'''
from sequence import *

foreground = readDistribs("ex1e-W30-iteration2-p-max.distrib")
background = readDistribs("ex1e-W30-iteration2-q-max.distrib")

print type(foreground), type(background)

pwm = PWM(foreground, background[0])

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