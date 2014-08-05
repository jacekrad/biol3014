'''
Created on 06/08/2014

@author: s4361277
'''
from sequence import *
from reportlab.lib.set_ops import intersect

foreground = readDistribs("ex1e-W30-iteration2-q-max.distrib")
q1e_background = readDistrib("ex1e-W30-iteration2-p-max.distrib")
hprd_background = readDistrib("HPRD.distrib")

# create two Position Wighted Matrices.  First used the background
# saved in question 1e and the second PWM is created using the 
# background from question 1b, ie from the HPRD sequence data
pwm1 = PWM(foreground, q1e_background)
pwm2 = PWM(foreground, hprd_background)

print "PWM1\n", pwm1, "\nPWM2\n", pwm2
ids1 = [] # ids with hits from PWM1 search
ids2 = [] # ids with hits from PWM2 search

hprd_sequences = readFastaFile('HPRD.fa', Protein_Alphabet)
hth_sequences = readFastaFile('hth_40.fa', Protein_Alphabet)
hth_ids = [seq.name for seq in hth_sequences]

# these are the sequences we are going to search.  These are the
# HPRD sequences that do not have an entry in our hth_40.fa FASTA 
# file, which was our training set
search_sequences = [seq for seq in hprd_sequences if not(seq.name in hth_ids)]

print "will search", len(search_sequences), "sequences out of a total of", \
    len(hprd_sequences), "entries found."

# because each sequence is a match we'll count the total number of hits
# for the two different sequences.  hits is not the same as the number of
# sequences matched as a sequence can be matched in more than one place
pwm1_hit_count = 0
pwm2_hit_count = 0

for sequence in search_sequences:
    hits1 = pwm1.search(sequence) # search using first PWM
    hits2 = pwm2.search(sequence) # search using second PWM
    pwm1_hit_count += len(hits1)
    pwm2_hit_count += len(hits2)
    # for both pwm1 and pwm2 search hits we only print the first of the hits
    # for each sequence as we were not asked to print them all.  There will
    # be sequences where the motif provides hits in more than one location
    if len(hits1) > 0:
        print "PWM1 hit: %s \t%d \t%s \t%5.3f" % (sequence.name, \
                                                  hits1[0][0], hits1[0][1], \
                                                  hits1[0][2])
        ids1.append(sequence.name)
    if len(hits2) > 0:
        print "PWM2 hit: %s \t%d \t%s \t%5.3f" % (sequence.name, \
                                                  hits2[0][0], hits2[0][1], \
                                                  hits2[0][2])
        ids2.append(sequence.name)
        
print "number of sequences matched by pwm1 and pwm2, and total unique", \
    len(ids1), len(ids2), len(intersect(ids1, ids2))
print "total hits by pwm1 and pwm2", pwm1_hit_count, pwm2_hit_count
print ids1
print "getting report"
report = getGOReport(ids1) 
print "report rows"
for row in report:
    print "%s \t%d \t%s" % row
