'''
Created on 05/08/2014

@author: jacekrad
'''
from sequence import * 
from prob import * 
from symbol import * 
from webservice import *
from sys import stderr

blosum62_background = readDistrib("blosum62.distrib")


def get_distrib_from_fasta(fasta_filename):
    '''
    This function creates a distribution from a set sequences
    read from a FASTA file.  All sequences are read and every
    residue is used to construct the distribution.  This
    function assumes that the FASTA file contains protein alphabet
    sequence.  The function returns a Distrib object containing the
    generated distribution.
    '''
    protein_counts = {}
    # initilise the dictionary with zeros
    for residue in Protein_Alphabet.symbols:
        protein_counts[residue] = 0
        
    sequences = readFastaFile(fasta_filename)
    for sequence in sequences:
        for letter in sequence:
            protein_counts[letter] += 1
            
    return Distrib(Protein_Alphabet, protein_counts)
        
hprd_background = get_distrib_from_fasta("HPRD.fa")

hprd_background.writeDistrib("HPRD.distrib")

print hprd_background
print blosum62_background

# compare the two distributions by calculating deltas for each amino acid
sys.stderr.write("AA,HPRD,BLOSUM62,DELTA\n")
for amino_acid in hprd_background.alpha:
    delta = hprd_background[amino_acid] - blosum62_background[amino_acid]
    sys.stderr.write(amino_acid + "," + str(hprd_background[amino_acid]) \
                     + "," + str(blosum62_background[amino_acid]) + "," \
                     + str(delta) + "\n")
    