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

for W in width_values:
    for i in range(1,4): # create 3 sets of results
        g = GibbsMotif(seqs, W)
        q = g.discover()
        p = g.getBackground()
        a = getAlignment(seqs, q, p)
        k = 0
        results_filename = "ex1e-W" + str(W) + "-iteration" + str(i) + ".aln"
        p_filename = "ex1e-W" + str(W) + "-iteration" + str(i) + "-p.aln"
        q_filename = "ex1e-W" + str(W) + "-iteration" + str(i) + "-q.aln"
        p.writeDistrib(p_filename)
        writeDistribs(q, q_filename)
        sys.stderr.write(results_filename + "\n")
        results_file = open(results_filename, 'w')
        for seq in seqs:
            results_file.write("%s \t%s\n" % (seq.name, seq[a[k]:a[k]+W]))
            k += 1
        results_file.close()