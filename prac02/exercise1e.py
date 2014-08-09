'''
Created on 07/08/2014

@author: jacekrad
'''
from gibbs import *
from sequence import *
import sys

seqs = readFastaFile("hth_40.fa", Protein_Alphabet)
width_values = [10, 24, 30]

# list of the files where we dump the results
# these will go to stderr for logo post processing
alignment_filenames = []

# maximum values for saving
maxLL = 0.0
max_p = None
max_q = None
p_filename = None
q_filename = None

for W in width_values:
    for i in range(1,4): # create 3 sets of results
        g = GibbsMotif(seqs, W)
        q = g.discover()
        p = g.getBackground()
        a = getAlignment(seqs, q, p)
        k = 0
        results_filename = "ex1e-W" + str(W) + "-iteration" + str(i) + ".aln"
        if g.maxLL > maxLL:
            maxLL = g.maxLL
            max_p = p
            max_q = q
            p_filename = "ex1e-W" + str(W) + "-iteration" + str(i) + "-p-max.distrib"
            q_filename = "ex1e-W" + str(W) + "-iteration" + str(i) + "-q-max.distrib"
            print "New maxLL distribution is ", q_filename
        sys.stderr.write(results_filename + "\n")
        results_file = open(results_filename, 'w')
        for seq in seqs:
            results_file.write("%s \t%s\n" % (seq.name, seq[a[k]:a[k]+W]))
            k += 1
        results_file.close()

# save distributions with highest log odds
print "Writing best distributions to ", p_filename, " and ", q_filename
max_p.writeDistrib(p_filename)
writeDistribs(max_q, q_filename)
