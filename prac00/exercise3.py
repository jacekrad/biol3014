'''
Created on 30/07/2014

@author: s4361277
'''
def myIsDNA (mySeq):
    '''
    This function returns True iff mySeq consists entirely of DNA letters or
    more specifically mySeq does not contain any of the following letters:
    QWERYUIOPSDFHJKLZXVBNM.  Function returns False otherwise.
    '''
    # list of non-DNA letters
    myNonDNA = "QWERYUIOPSDFHJKLZXVBNM" 
    
    # iterate over the non-DNA letters
    for myAA in myNonDNA:
        
        # check if the non-DNA letter is contained in our sequence
        # and if it is then return False
        if myAA in mySeq:
            return False
        
    # finished iterating over the non-DNA letter and have not returned False yet so return True as mySeq IS DNA
    return True

# call the function
print myIsDNA("GTTCGACCA")