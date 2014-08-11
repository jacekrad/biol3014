'''
Created on 06/08/2014

@author: jacekrad
'''
from gibbs import *
from sequence import *

seqs = readFastaFile("hth_40.fa", Protein_Alphabet)

W = 10 # the width of the motif sought

# create a GibbsMotif object from a list of sequences
# and of length W    
g = GibbsMotif(seqs, W)
    
# execute the core Gibbs Sampling algorithm to discover 
# the motif          
q = g.discover()
   
# get the probability distribution for the background used
# in the discovery calculated above
p = g.getBackground()

# getAlignments is called and alignment for sequences seq
# is calculated from the foreground q and background p
# the resulting alignment is assigned to a
a = getAlignment(seqs, q, p)

k = 0
for seq in seqs:
    print "%s \t%s" % (seq.name, seq[a[k]:a[k]+W])
    k += 1
