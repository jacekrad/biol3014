'''
Created on 30/07/2014

@author: s4361277
'''
myHydrophobics = "FLIMVPAWG"
mySequence = ["M","H","K","L"]
for myaa in mySequence:
    print "The amino acids is",myaa
    if myaa in myHydrophobics:
        print " It is hydrophobic"
        print "End of the program"