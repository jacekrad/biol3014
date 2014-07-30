'''
Created on 30/07/2014

@author: s4361277
'''
myHydrophobics = "FLIMVPAWG"
mySequence = "MHKL"
mySequenceHydro = ""
for myaa in mySequence:
    print "The amino acids is",myaa
    if myaa in myHydrophobics:
        mySequenceHydro += "*"
        #print " It is hydrophobic"
    else:
        mySequenceHydro += " "
print mySequence
print mySequenceHydro
#print "End of the program"