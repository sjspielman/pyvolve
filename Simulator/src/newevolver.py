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
        
        #####################################################################
        # Note that these 3 lines only make sense when there are no indels. #
        self.seqlen = 0
        for i in range(self.numparts):
            self.seqlen += self.parts[i][0]  
        #####################################################################     
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
            raise AssertionError("wut.")

    def writeSequences(self, **kwargs):
        ''' Write resulting sequences to a file, currently only in fasta format. '''
        outfile  = kwargs.get("outfile", "seqs_"+strftime("%m.%d.%H.%M.%S")+".fasta")  
        out_handle=open(outfile, 'w')
        for entry in self.alndict:
            seq = self.intseq_to_string(self.alndict[entry])
            out_handle.write(">"+entry+"\n"+seq+"\n")
        out_handle.close()               
            
    def seq2int(self, entry):
        ''' Take a DNA/PROT value and return its integer index '''
        return self.code.index(entry)
    
    
    def int2seq(self, index):
        ''' Take a DNA/PROT index and return its corresponding molecule ''' 
        return self.code[index]
    
    ######### FUNCTION WILL LIKELY NEED OVERHAUL TO ACCOUNT FOR SITE CLASS ###############
    def intseq_to_string(self, intseq):
        ''' Take a sequence coded as ints and turn to actual molecule string '''
        stringseq = ''
        for i in intseq:
            stringseq += self.int2seq(i)
        return stringseq   
        

    def checkParentBranch(self, node, parentSeq, parentModel):
        ''' Check that the parent sequence exists, branch length is reasonable, and assign a model. ''' 
        assert (parentSeq != None), "There is no parent sequence."

        branchModel = node.modelFlag
        if branchModel is None:
            node.modelFlag = parentModel

        branchLength = float( node.branch )
        assert (branchLength >= 0.), "Branch length is negative. Must be >= 0."

        return branchLength, node.modelFlag           



    def generateSite(self, nodeID, state):
        ''' Generate a new Site object when generating the root (state 0) and insertion (state 1) sites. '''
        site        = misc.Site()
        site.origin = nodeID
        site.state  = state
        return site
           
           

    def generateRootSeq(self):
        ''' Select starting sequence based on state frequencies, for each partition, and return full root sequence. '''
        
        currentNode.seq = []
        for n in range(self.numparts):
            partlen = self.parts[n][0]
            partSites = []
            freqs = self.parts[n][1][self.rootModel].substParams['stateFreqs']
            index = 0
            for j in range(partlen):
                newsite = generateSite(rootID, 0)
                newsite.intSeq = self.generateSeq(freqs)
                partSites[index] = newsite
                index += 1
            currentNode.seq.append(partSites)



    def generateSeq(self, probArray):
        ''' Sample a sequence letter (nuc,aa,or codon). probArray can be any list/numpy array of probabilities that sum to 1.'''
        
        assert ( abs(np.sum(probArray) - 1.) < self.zero), "Probabilities do not sum to 1. Cannot generate a new sequence."
        r = rn.uniform(0,1)
        i=0
        sum=probArray[i]
        while sum < r:
            i+=1
            sum+=probArray[i]
        return i        
        

    def simulate(self, currentNode, parentNode = None):
        ''' Traverse the tree and simulate. '''

        # We are at the base and must generate root sequence
        if (parentNode is None):
            self.generateRootSeq( currentNode )
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
        
            newSeq = [] # Since indels, must be dynamic list.
            seqIndex = 0
            for n in range(self.numparts):
            
                # Create partition length and model variables, for clarity
                partlen = len(self.parts[n][0])
                model = self.parts[n][1][branchModel]
                
                # Simulate substitutions. This will also create the partSites list for this partition.
                partSites, mean_sub_rate = self.simSubst(partlen, model, parentSeq, branchLength)
            
                # TO DO: Simulate indels. Probably will be something like this:
                # partSites = self.simIndel(partlen, model, partSites, branchLength)
                
                newSeq.append(partSites)
                
        # Attach final newSeq list of evolved partitions to node
        node.seq = newSeq 
            
    def simIndel(self, length, model, mean_sub_rate, partSites):
        ''' go, indels, go go go. Gillespie. '''
        
        # Scale insertion and deletion rates by this branch's mean_sub_rate
        # model.indelParams = the dictionary of indel parameters. insRate, delRate, insDist, delDist, etc...
        
        # 
        
        
             
            
    def simSubst(self, length, model, parentSeq, branchLength):  
        ''' Simulate substitution process along a branch for a single partition.
            Here, we will only be tweaking the Site.intSeq attributes (no others!!).
         '''
                                    
        # Generate probability matrix for evolution along this branch and assert correct
        Q = model.Q
        Qt = np.multiply(Q, branchLength)
        probMatrix = linalg.expm( Qt )
        for i in range(len(self.code)):
            assert( abs(np.sum(probMatrix[i]) - 1.) < self.zero ), "Row in P(t) matrix does not sum to 1."
        
        # Update as we evolve. This value will be subsequently used for indels.
        mean_sub_rate = 0.
    
        # March along parentSeq and substitute, while also saving cumulative subprob
        partSites = []
        for j in range( length ):
            
            # Assign parent position to current node.
            partSites[j] = parentSeq[j]
            
            # If position is not a gap, evolve it
            else:
                partSites[j].intSeq = self.generateSeq( probMatrix[ partSites[j].intSeq ] )     
                
                # Keep building this up so can get a final mean
                mean_sub_rate += np.sum(probMatrix[ parentSeq[j].intSeq ]) - probMatrix[intSeq]
        
        return partSites, mean_sub_rate /= length
            
            
            
            
            
            
            
            
            
            
            