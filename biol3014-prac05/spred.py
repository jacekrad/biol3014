'''
A module to enable experimentation with various methods for predicting properties
assigned to sequence elements, e.g. secondary structure of proteins.
A neural net wrapper class is provided.
A couple of example applications are found at the end of this module.
'''
import numpy
import symbol
import prob
import sequence
import ml

def slidewin(seq, winsize):
    """ Produce a list of sub-sequences of a given length from a complete sequence """
    subseqs = []
    for i in range(len(seq) - winsize + 1):
        subseqs.append(seq[i : i + winsize])
    return subseqs

def _onehotIndex(alpha, sym):
    """ Create array with "one-hot" bit codes (only adding "ones" to an all-"zero" array) """
    symlen = len(sym)
    alphalen = len(alpha)
    indices = [ alpha.index(sym[i]) + (i * alphalen) for i in range(symlen) ]
    return indices

class SeqNN():
    """ A neural net wrapper for multinomial classification of sequence input """
    
    def __init__(self, inp_len, inp_alpha, outp_alpha, nhidden, cascade = 0):
        """ Construct a neural net with numeric inputs and outputs 
            depending on alphabets used for inputs and outputs.
            inp_len: number of symbols to use as input
            inp_alpha: input alphabet
            outp_alpha: output alphabet (defines number of classes)
            nhidden: number of "hidden" nodes in the net
            cascade: if non-zero, number of positions to feed into a cascaded structure-to-structure NN (also the number of hidden nodes of this NN)
        """
        self.nn1 = ml.NN(inp_len * len(inp_alpha), nhidden, len(outp_alpha))  # neural net
        self.nn2 = None
        self.cascade = cascade
        if cascade > 0:
            self.nn2 = ml.NN(cascade * len(outp_alpha), cascade, len(outp_alpha))  # cascaded neural net
        self.inp_len  = inp_len
        self.inp_alpha  = inp_alpha
        self.outp_alpha = outp_alpha
    
    def _encodeseq(self, seqs, targets = None):
        """ Convert a list of sequences into numeric input suitable as input to NN. """
        try:
            len(seqs[0]) # if this does not throw error, it is a multi-input already
        except TypeError:
            seqs = [ seqs ]
            targets = [ targets ]
        totlen = 0
        alpha = None
        for seq in seqs:
            if not alpha:
                alpha = seq.alphabet
            totlen += len(seq) - self.inp_len + 1
        im = numpy.zeros((totlen, self.inp_len * len(alpha)))
        if targets:
            om = numpy.zeros((totlen, len(self.outp_alpha)))
        row = 0 
        for i in range(len(seqs)):
            subseqs = slidewin(seqs[i], self.inp_len)
            if targets:
                # Note how we remove the targets at the ends of the sequence
                subtarg = targets[i][self.inp_len/2:-self.inp_len/2+1]
            for k in range(len(subseqs)):
                im[row, _onehotIndex(alpha,  subseqs[k])] = 1
                if targets: om[row, self.outp_alpha.index(subtarg[k])] = 1
                row += 1
        print "There are", row, "entries in data set"
        if targets:
            return im, om
        else:
            return im, None
    
    def observeAll(self, seqs, targets, eta = 0.1, niter = 1):
        """ Train a classifier to map from all possible windows to the target symbols.
            Decompose each sequence to all full-width sub-sequences. Map each sub-sequence
            to the target symbol for the symbol in the centre of the sub-sequence. """
        assert len(seqs) == len(targets), "Number of input sequences need to match the number of target sequences"
        im, om = self._encodeseq(seqs, targets)
        for i in range(niter):  # train first NN
            rmse = self.nn1.train(im, om, eta = eta, niter = 1)
            print i, ":", rmse
        if not self.cascade:    # if there's no cascaded NN, finish here
            return rmse
        nn1seqs = []            # a list of new SS sequences ...
        for seq in seqs:        # ... based on AA sequences
            nn1seq = self.predict(seq, useCascade = False) # construct a new sequence which consists of SS predictions
            nn1seqs.append(nn1seq)
        im, om = self._encodeseq(nn1seqs, targets)  # construct input/output patterns from SS sequences
        for i in range(niter):  # train cascaded NN
            rmse = self.nn2.train(im, om, eta = eta, niter = 1)
            print i, ":", rmse
        return rmse
    
    def testAll(self, seqs, targets):
        """ Test the neural network on the specified sequences and target sequences.
            Returns a confusion matrix with the predictions. """
        assert len(seqs) == len(targets), "Number of input sequences needs to match the number of target sequences"
        if not self.cascade:
            im, om = self._encodeseq(seqs, targets)
            cm = self.nn1.test(im, om)
            return cm
        else:
            nn1seqs = []
            for seq in seqs:
                nn1seq = self.predict(seq, useCascade = False)
                nn1seqs.append(nn1seq)
            im, om = self._encodeseq(nn1seqs, targets)
            cm = self.nn2.test(im, om)
            return cm
    
    def predict(self, inpseq, useCascade = True):
        """ Classify each symbol in a sequence.
            Return the predictions as a list of symbols. """
        W = self.nn1.ninput / len(self.inp_alpha)
        if useCascade and self.cascade:
            nn1seq = self.predict(inpseq, useCascade = False)
            subseqs = slidewin(nn1seq, self.cascade)
            predsyms = ['C' for _ in range(len(inpseq))] # use coil for positions in flanking regions 
            for i in range(len(subseqs)):    # for each input sub-sequence of the primary NN
                input = numpy.zeros(self.cascade * len(self.outp_alpha))
                input[_onehotIndex(self.outp_alpha, subseqs[i])] = 1
                outvec = self.nn2.feedforward(input)
                d = prob.Distrib(self.outp_alpha)
                for k in range(len(outvec)):
                    d.observe(self.outp_alpha[k], outvec[k])
                predsyms[i + self.cascade / 2] = d.getmax()    # use the symbol with the highest probability
            return sequence.Sequence(predsyms, self.outp_alpha)
        else: # only predict using the first NN
            subseqs = slidewin(inpseq, W)
            predsyms = ['C' for _ in range(len(inpseq))] # use coil for positions in flanking regions 
            for i in range(len(subseqs)):    # for each input sub-sequence of the primary NN
                input = numpy.zeros(self.inp_len * len(self.inp_alpha))
                input[_onehotIndex(self.inp_alpha, subseqs[i])] = 1
                outvec = self.nn1.feedforward(input)
                d = prob.Distrib(self.outp_alpha)
                for k in range(len(outvec)):
                    d.observe(self.outp_alpha[k], outvec[k])
                predsyms[i + W / 2] = d.getmax()    # use the symbol with the highest probability
            return sequence.Sequence(predsyms, self.outp_alpha)

##################################################################################################################
#      Example applications of ML methods including NNs and Naive Bayes for secondary structure prediction.      #
##################################################################################################################

if __name__=='__main__': # examples to run unless this module is merely "imported"
    import os, time
    os.chdir('/Users/mikael/workspace/binf/data')  # Note you will need to change this to find your directory of choice
    prot = sequence.readFastaFile('prot2.fa', symbol.Protein_Alphabet)  # proteins
    sstr = sequence.readFastaFile('sstr3.fa', symbol.DSSP3_Alphabet)    # secondary structure of prot
    # separate training and test data
    prot_trn = prot[0::2] # even-numbered indices
    prot_tst = prot[1::2] # odd-numbered indices
    sstr_trn = sstr[0::2] # even-numbered indices
    sstr_tst = sstr[1::2] # odd-numbered indices
    W = 15

if __name__=='__main__':   # NN (should read "__main__" for it to be executed on "Run")  
    nHid = 30
    nn = SeqNN(W, symbol.Protein_Alphabet, symbol.DSSP3_Alphabet, nHid, cascade = W)
    #nn.nn = ml.readNNFile('sstr3.nn')
    #print "Successfully loaded network"
    start = time.time()
    print nn.observeAll(prot_trn, sstr_trn, eta = 0.01, niter = 20)
    print "It took", time.time() - start, "seconds to train NN."
    cm = nn.testAll(prot_tst, sstr_tst)
    print cm
    Qk, Q = ml.Qk(cm, symbol.DSSP3_Alphabet)
    print Q
    print Qk
    nn.nn1.writeFile('sstr3.nn1')
    nn.nn2.writeFile('sstr3.nn2')
    print "Successfully saved network/s"

if __name__=='__main__REMOVE_ME':    # NB (should read "__main__" for it to be executed on "Run")
    nb = prob.NaiveBayes([ symbol.Protein_Alphabet for _ in range(W) ], symbol.DSSP3_Alphabet)
    start = time.time()
    for i in range(len(prot_trn)):
        subseqs = slidewin(prot_trn[i], W)
        # Note how we remove the targets at the ends of the sequence
        subtarg = sstr_trn[i][W/2:-W/2+1]
        for j in range(len(subseqs)):
            nb.observe(subseqs[j], subtarg[j])
    print "It took", time.time() - start, "seconds to train NB."
    cm = numpy.zeros((len(symbol.DSSP3_Alphabet), len(symbol.DSSP3_Alphabet)))
    for i in range(len(prot_tst)):
        subseqs = slidewin(prot_tst[i], W)
        # Note how we remove the targets at the ends of the sequence
        subtarg = sstr_tst[i][W/2:-W/2+1]
        for j in range(len(subseqs)):
            out = nb[subseqs[j]]
            c_targ = symbol.DSSP3_Alphabet.index(subtarg[j])  
            c_pred = out.prob().index(max(out.prob()))
            cm[c_targ, c_pred] += 1
    print cm
    Qk, Q = ml.Qk(cm, symbol.DSSP3_Alphabet)
    print Q
    print Qk
        