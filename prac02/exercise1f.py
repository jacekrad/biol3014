'''
Created on 06/08/2014

@author: s4361277
'''
from sequence import *
from reportlab.lib.set_ops import intersect

foreground = readDistribs("ex1e-W30-iteration2-p-max.distrib")
q1e_background = readDistribs("ex1e-W30-iteration2-q-max.distrib")
hprd_background = readDistribs("HPRD.distrib")

print type(foreground), type(q1e_background), type(hprd_background)

# create two Position Wighted Matrices.  First used the background
# saved in question 1e and the second PWM is created using the 
# background from question 1b, ie from the HPRD sequence data
pwm1 = PWM(foreground, q1e_background[0])
pwm2 = PWM(foreground, hprd_background[0])

print pwm1, pwm2
ids1 = []
ids2 = []
hprd_sequences = readFastaFile('HPRD.fa', Protein_Alphabet)
hth_sequences = readFastaFile('hth_40.fa', Protein_Alphabet)
hth_ids = [seq.name for seq in hth_sequences]

# these are the sequences we are going to search.  These are the
# HPRD sequences that do not have an entry in our hth_40.fa FASTA file
search_sequences = [seq for seq in hprd_sequences if not(seq.name in hth_ids)]

print "will search", len(search_sequences), "sequences out of a total of", \
    len(hprd_sequences), "entries found."

# because each sequence is a match we'll count the total number of hits
# for the two different sequences

pwm1_hit_count = 0
pwm2_hit_count = 0
pwm1_theshhold = 4.5
count = 0 
for sequence in search_sequences:
    hits1 = pwm1.search(sequence) # search using first PWM
    hits2 = pwm2.search(sequence) # search using second PWM
    pwm1_hit_count += len(hits1)
    pwm2_hit_count += len(hits2)
    if len(hits1) > 0:
        #print "%s \t%d \t%s \t%5.3f" % (sequence.name, hits1[0][0], hits1[0][1], hits1[0][2])
        ids1.append(sequence.name)
    #if len(hits2) > 0:
    #    ids2.append(sequence.name)
        
#print "number of sequences matched by pwm1 and pwm2", len(ids1), len(ids2), len(intersect(ids1, ids2))
#print "total hits by pwm1 and pwm2", pwm1_hit_count, pwm2_hit_count
print ids1
print "getting report"
report = getGOReport(ids1) 
print "report rows"
for row in report:
    print "%s \t%d \t%s" % row
