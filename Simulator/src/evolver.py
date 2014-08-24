import misc
import numpy as np
from scipy import linalg
import random as rn
from time import strftime
        
        
class Evolver(object):
    def __init__(self, partitions, rootModel):
        
        self.zero = 1e-8
        self.molecules = misc.Genetics()
        
        self.parts    = partitions   # this should be a list of tuples. Each tuple is (length, {flag:model, ...}). 
        self.numparts = len(self.parts)
        assert ( self.numparts >= 1), "You have nothing to evolve. Partitions, please!"
        self.rootModel = rootModel # some checking should probably be done here.
        self.alndict = {}
        self.seqlen = 0
        for i in range(self.numparts):
            self.seqlen += self.parts[i][0]       
        self.setCode()

        
    def setCode(self):
        ''' Determine whether we have nuc, amino, or codon alphabet and assign a genetic code. '''   
        # An incredible hack! Go python, go!!
        dim = self.parts[0][1].values()[0].Q.shape[0]
        
        if dim == 4:
            self.code = self.molecules.nucleotides
        elif dim == 20:
            self.code = self.molecules.amino_acids
        elif dim == 61:
            self.code = self.molecules.codons
        else:
            raise AssertionError("Matrix is totally weird.")
            
            
    def seq2int(self, entry):
        ''' Take a DNA/PROT value and return its integer index '''
        return self.code.index(entry)
    
    
    def int2seq(self, index):
        ''' Take a DNA/PROT index and return its corresponding molecule ''' 
        return self.code[index]
    
    
    def intseq_to_string(self, intseq):
        ''' Take a sequence coded as ints and turn to actual molecule string '''
        stringseq = ''
        for i in intseq:
            stringseq += self.int2seq(i)
        return stringseq   

    
    def generateSeq(self, probArray):
        ''' Sample a sequence letter (nuc,aa,or codon). probArray can be any list/numpy array of probabilities that sum to 1.'''
        r = rn.uniform(0,1)
        i=0
        sum=probArray[i]
        while sum < r:
            i+=1
            sum+=probArray[i]
        return i        

    def generateRootSeq(self):
        ''' Select starting sequence based on state frequencies, for each partition, and return full root sequence. '''
        
        rootSeq = np.zeros(self.seqlen, dtype=int)
        rootSeq[rootSeq == 0] = 100
        completed = 0
        for i in range(self.numparts):
            partlen = self.parts[i][0]
            freqs  = self.parts[i][1][self.rootModel].params['stateFreqs']
            root_iter = np.nditer(rootSeq[completed:completed+partlen], op_flags=['writeonly'], flags = ['c_index'])
            for x in root_iter:
                x[...] = self.generateSeq(freqs)
            completed += partlen
        assert( np.sum(rootSeq == 100) == 0), "Root sequence improperly evolved. womp womp."
        return rootSeq 
        

    def checkParentBranch(self, node, parentSeq, parentModel):
        ''' Check that the parent sequence exists, branch length is reasonable, and assign a model. ''' 
        assert (parentSeq != None), "There is no parent sequence."

        branchModel = node.modelFlag
        if branchModel is None:
            node.modelFlag = parentModel

        branchLength = float( node.branch )
        assert (branchLength >= 0.), "Branch length is negative. Must be >= 0."

        return branchLength, node.modelFlag           




    def writeSequences(self, **kwargs):
        ''' Write resulting sequences to a file, currently only in fasta format.'''
        
        outfile  = kwargs.get("outfile", "seqs_"+strftime("%m.%d.%H.%M.%S")+".fasta")  
        out_handle=open(outfile, 'w')
        for entry in self.alndict:
            seq = self.intseq_to_string(self.alndict[entry])
            out_handle.write(">"+entry+"\n"+seq+"\n")
        out_handle.close()    
        
        
        
        
    def evolveBranch(self, node, parentNode):
        ''' Crux function to evolve sequences along a branch.'''
    
        # Ensure brank length ok, parent sequence exists, and model is assigned.
        parentSeq = parentNode.seq
        parentModel = parentNode.modelFlag
        branchLength, branchModel = self.checkParentBranch(node, parentSeq, parentModel)

        # Evolve only if branch length is greater than 0.
        if branchLength <= self.zero:
            newSeq = parentSeq
        else:
            completed = 0
            for i in range(self.numparts):
            
                # set the length and the instantaneous rate matrix for this partition at this node
                partlen  = self.parts[i][0]
                instMat = self.parts[i][1][branchModel].Q
                
                # Generate probability matrix for evolution along this branch and assert correct
                Qt = np.multiply(instMat, branchLength)
                probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
                assert( -1.*self.zero < np.sum(probMatrix, axis=1).all() - 1. < self.zero )
                
                # Now evolve
                parent_mirror = parentSeq[completed:completed+partlen]
                for x in seq_iter: 
                   x[...] = self.generateSeq(probMatrix[ parent_mirror[seq_iter.index] ] ) 
                completed += partlen  
                             
        # Attach final sequence to node
        node.seq = parentSeq 



    def simulate(self, currentNode, parentNode = None):
        ''' Traverse the tree and simulate. '''

        # We are at the base and must generate root sequence
        if (parentNode is None):
            currentNode.seq = self.generateRootSeq() 
            currentNode.modelFlag = self.rootModel 
        else:
            self.evolveBranch(currentNode, parentNode)
            
            
        # We are at an internal node. Keep evolving
        if len(currentNode.children)>0:
            for childNode in currentNode.children:
                self.simulate(childNode, currentNode)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[currentNode.name]=currentNode.seq
            