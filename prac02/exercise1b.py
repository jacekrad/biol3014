'''
Created on 05/08/2014

@author: jacekrad
'''
from sequence import * 
from prob import * 
from symbol import * 
from webservice import *

blosum62_background = readDistrib("blosum62.distrib")


def get_distrib_from_fasta(fasta_filename):
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