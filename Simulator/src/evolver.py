import numpy as np
import dendropy as dp
from scipy import linalg
import random as rn
import misc
import time
        
class Evolver(object):
    def __init__(self, **kwargs):
        
        #Provided by user
        self.parts    = kwargs.get("partitions", None)   # this should be a list of tuples. Each tuple is (length, model). 
        self.numparts = len(self.parts)
        assert (len(self.numparts) >= 1), "You have nothing to evolve. Partitions, please!"
        self.tree     = kwargs.get("tree", None)
        assert (self.tree is not None), "You have no tree. Treeeeeee!"
        self.outfile  = kwargs.get("outfile", "seqs_"+time.strftime("%m.%d.%H.%M.%S")+".fasta")        
        self.ratefile = kwargs.get("ratefile", None)
        
        self.seqlen   = 0
        for i in range(self.numparts):
            self.seqlen += self.parts[i][0]
                
        # Internals
        self.alndict = {}
        self.zero = 1e-8
        self.molecules = misc.Genetics()
        
        dimensions = self.parts[0][1].Q.shape[0]
        if dimensions == 4:
            self.code = self.molecules.nucleotides
        elif dimsensions == 20:
            self.code = self.molecules.amino_acids
        elif dimensions == 61:
            self.code = self.molecules.codons
        else:
            raise AssertionError("Matrix is totally weird.")
        
        
        

    def seq2int(self, entry):
        ''' Take a DNA/PROT value and return its integer index 0-60 '''
        return self.code.index(entry)
    
    
    def int2seq(self, index):
        ''' Take a DNA/PROT index (0-60) and return its corresponding letter/codon ''' 
        return self.code[index]
    
    
    def intseq_to_string(self, intseq):
        ''' Take a sequence coded as ints and turn to actual codon string '''
        stringseq = ''
        for i in intseq:
            single_seq = self.int2seq(i)
            stringseq += single_seq
        return stringseq
    
        
    def generateRootSeq(self):
        ''' Select starting sequence based on state frequencies, for each partition. '''
        rootSeq = np.empty(self.seqlen, dtype=int)
        index=0
        for i in range(self.numparts):
            seqlen = self.parts[i][0]
            freqs  = self.parts[i][1].params["stateFreqs"]
            for j in range(seqlen):
                rootSeq[index] = self.generateSeq(freqs)
                index += 1
        return rootSeq    


    def sim_sub_tree(self, tree, baseSeq = None):
        ''' Traverse the tree and simulate. '''
        
        # We are at the base and must generate root sequence
        if (baseSeq is None):
            tree.seq = self.generateRootSeq()        
        else:
            self.evolve_branch(tree, baseSeq)
                
        # We are at an internal node. Keep evolving
        if len(tree.children)>0:
            for node in tree.children:
                self.sim_sub_tree(node, tree.seq)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[tree.name]=tree.seq
    
    
    def generateSeq(self, probArray):
        ''' Sample a sequence letter (nuc,aa,or codon). probArray can be any list/numpy array of probabilities that sum to 1.'''
        #### CHECKED FXN ON 2/6/14. WORKS AS INTENDED #####
        # Assertion is overkill but who cares
        assert ( abs(np.sum(probArray) - 1.) < self.zero), "Probabilities do not sum to 1. Cannot generate a molecule."
        r = rn.uniform(0,1)
        i=0
        sum=probArray[i]
        while sum < r:
            i+=1
            sum+=probArray[i]
        return i        

    def checkParentBranch(self, node, baseSeq):
        ''' Check that the baseSeq exists and the branch length is reasonable ''' 
        assert (baseSeq != None), "There is no parent sequence."
        bl = float( node.branch )
        assert (bl >= 0), "Branch length is negative. Must be >= 0."
        return bl            


    def writeSequences(self):
        ''' Write resulting sequences to a file'''
        print "Writing sequences to file, in fasta format."
        out_handle=open(self.outfile, 'w')
        for entry in self.alndict:
            seq = self.intseq_to_string(self.alndict[entry])
            out_handle.write(">"+entry+"\n"+seq+"\n")
        out_handle.close()    
        

    ##################### THIS FUNCTION IS TOTALLY DEPRECATED ############################
    def printRates(self):
        ''' Print dN/dS to file for each site. ASSUMES THAT SITES ARE NOT SHUFFLED. ''' 
        
        f = open(outfile, 'w')
        index = 0
        for part in partitions:
            partLen = part[0]
            partRate = str( part[1].params["omega"] )
            for i in range(partLen):
                f.write(str(index)+'\t'+partRate+'\n')
                index+=1
        f.close()    
    ######################################################################################
    
        
    def evolve_branch(self, node, baseSeq):
        '''Base class function. Not implemented.'''
        return 0
    


class StaticEvolver(Evolver):
    ''' Evolve according to a static landscape (no temporal variation). Site heterogeneity ok, but no temporal/branch heterogeneity. ''' 
    
    def __init__(self, **kwargs):
        super(StaticEvolver, self).__init__(**kwargs)
    
    def evolve_branch(self, node, baseSeq):
        
        bl = self.checkParentBranch(node, baseSeq)
        
        # If there is no branch length then there is nothing to evolve. Attach baseSeq to node
        if bl <= self.zero:
            print bl, "branch length of 0 detected"
            node.seq = baseSeq
        
        else:
        
            # Evolve each partition, i, and add onto newSeq as we go
            newSeq = np.empty(self.seqlen, dtype=int)
            index = 0
            for i in range(self.numparts):
            
                # set the length and the instantaneous rate matrix for this partition
                seqlen  = self.parts[i][0]
                instMat = self.parts[i][1].Q
                
                # Generate probability matrix for evolution along this branch and assert correct
                Qt = np.multiply(instMat, bl)
                probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
                for i in range(len(self.code)):
                    assert( abs(np.sum(probMatrix[i]) - 1.) < self.zero ), "Row in P(t) matrix does not sum to 1."
    
                # Move along baseSeq and evolve. 
                for j in range(seqlen):
                    newSeq[index] = self.generateSeq( probMatrix[baseSeq[index]] )
                    index+=1
                    
            # Attach final sequence to node
            node.seq = newSeq 
            
            
#### CREATED BUT NOT DIFFERENT AT THIS TIME FROM STATICEVOLVER CLASS ####
class ShiftingEvolver(Evolver):
    ''' Evolve according to a changing landscape (temporal variation) ''' 
    def __init__(self, **kwargs):
        super(ShiftingEvolver, self).__init__(**kwargs)

    def evolve_branch(self, node, baseSeq):
        
        bl = self.checkParentBranch(node, baseSeq)
        
        # If there is no branch length then there is nothing to evolve. Attach baseSeq to node
        if bl < self.zero:
            print bl, "branch length of 0 detected"
            node.seq = baseSeq
        
        else:
            ## Evolve for each partition and then join together
            newSeq = np.empty(self.seqlen, dtype=int)
            index = 0
            for i in range(self.numparts):
            
                # set the length and the instantaneous rate matrix for this partition
                seqlen  = self.parts[i][0]
                instMat = self.parts[i][1].Q
                
                # Generate probability matrix for evolution along this branch and assert correct
                Qt = np.multiply(instMat, bl) # Matrix has already been scaled properly.
                probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
                for i in range(61):
                    assert( abs(np.sum(probMatrix[i]) - 1.) < self.zero ), "Row in P(t) matrix does not sum to 1."
    
                # Move along baseSeq and evolve. 
                for j in range(seqlen):
                    newSeq[index] = self.generateSeq( probMatrix[baseSeq[index]] )
                    index+=1
                    
            # Attach final sequence to node
            node.seq = newSeq 







