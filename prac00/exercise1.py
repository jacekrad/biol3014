'''
Created on 30/07/2014

@author: s4361277
'''

my_aa = "Z"
my_aacodes = {'Y':'TYR', 'F':'PHE', 'C':'CYS', 'V':'VAL', 'H':'HIS', 'L':'LEU', 'D':'ASP','A':'ALA', 'S':'SER', 'E':'GLU', 'W':'TRP', 'P':'PRO', 'I':'ILE', 'R':'ARG', 'K':'LYS','N':'ASN', 'M':'MET', 'T':'THR', 'G':'GLY', 'Q':'GLN'}


print "The provided letter is ", my_aa
if my_aa in my_aacodes:
    print "This letter is in the list, it corresponds to the amino acid", my_aacodes.get(my_aa)
else:
    print "This letter is not in the list and therefore does not correspond to an amino acid"