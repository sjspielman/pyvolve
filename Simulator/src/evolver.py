import misc
from indel import *
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
    ######################################################################################   

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
                newsite = generateSite(currentNode.name, 0)
                newsite.intSeq = self.generateUnifProb(freqs)
                partSites[index] = newsite
                index += 1
            currentNode.seq.append(partSites)



    def generateUnifProb(self, probArray):
        ''' Sample a sequence letter (nuc,aa,or codon). probArray can be any list/numpy array of probabilities that sum to 1.'''
        assert ( abs(np.sum(probArray) - 1.) < self.zero), "Probabilities do not sum to 1. Cannot generate a new sequence."
        r = rn.uniform(0,1)
        i=0
        sum=probArray[i]
        while sum < r:
            i+=1
            sum+=probArray[i]
        return i        
        

    def simulate(self, currentNode, parentNode = None, waitTime = None):
        ''' Traverse the tree and simulate. '''

        # We are at the base and must generate root sequence
        if (parentNode is None):
            self.generateRootSeq( currentNode )
            currentNode.modelFlag = self.rootModel 
        else:
            waitTime = self.evolveBranch(currentNode, parentNode, waitTime)
            
        # We are at an internal node. Keep evolving
        if len(currentNode.children)>0:
            for childNode in currentNode.children:
                # only give waitTime to one child.
                if childNode is currentNode.children[0]:            
                    self.simulate(childNode, currentNode, waitTime)
                else:
                    self.simulate(childNode, currentNode)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[currentNode.name]=currentNode.seq
            
          
             
            
    def simSubst(self, parentSeq, model, insertedSites, branchLength):
        ''' Simulate substitution process along a branch for a single partition.
            Here, we will only be tweaking the Site.intSeq attributes (no others!!).
         '''
        ###################### NEEDS OVERHAUL TO EVOLVE NEWLY INSERTED BITS SEPARATELY FROM PREVIOUSLY PRESENT BITS #############################
                                    
        # Generate probability matrix for evolution along this branch and assert correct
        Q = model.Q
        Qt = np.multiply(Q, branchLength)
        probMatrix = linalg.expm( Qt )
        for i in range(len(self.code)):
            assert( abs(np.sum(probMatrix[i]) - 1.) < self.zero ), "Row in P(t) matrix does not sum to 1."
     
        
        # March along parentSeq and substitute
        partSites = []
        for j in range( parentSeq.numsites ):
            
            # Assign parent position to current node.
            partSites[j] = parentSeq[j]
            
            # If position is not a gap, evolve it
                partSites[j].intSeq = self.generateUnifProb( probMatrix[ partSites[j].intSeq ] )     

        return partSites   
        
        
    def evolveBranch(self, node, parentNode, waitTime = None):
        ''' Crux function to evolve sequences along a branch. '''
    
        # Ensure brank length ok, parent sequence exists, and model is assigned.
        parentSeq = parentNode.seq
        parentModel = parentNode.modelFlag
        branchLength, branchModel = self.checkParentBranch(node, parentSeq, parentModel)

        # Evolve only if branch length is greater than 0.
        if branchLength <= self.zero:
            newSeq = parentSeq
        else:
        
            newSeq = [] # Since indels, must be dynamic list.
            for n in range(self.numparts):
                model = self.parts[n][1][branchModel] # variable for clarity
                
                # Simulate indels and create the partSites for this partition
                partSites, insertedSites, nextWaitTime = self.simIndel(parentSeq, node.name, model, branchLength, waitTime)
                
                # Simulate substitutions.
                partSites = self.simSubst(parentSeq, model, insertedSites, branchLength)
            
                
                newSeq.append(partSites)
                
        # Attach final newSeq list of evolved partitions to node
        node.seq = newSeq 
        
        return nextWaitTime
        
        
        
    ############################### INDEL-SPECIFIC CODE ##########################################
            
    def calcProbIndelEvent(insSpace, delSpace, insRate, delRate, meanDelSize):
        ''' calculate probability of indel event occuring, which is a fxn of the current sequence length.
            NOTE: as length will stay the same from a parent to child branch, it's ok that things get carried over.
            # current scheme, as of 7/7/14: P_event = rate_ins*(length+1) + rate_del*(length) . <- iSG uses this.
            # ME: only real sequence bits can be deleted, as in cannot delete a gap. Therefore, the probability of deletions depends not on a naive length, but number of sequence sites that there are. Otherwise will overestimate deletion rate.
        '''
        insSpace = float(insSpace)
        delSpace = float(delSpace)
        p_event = insRate*(insSpace+1.) + delRate*(delSpace + meanDelSize + 1.) 
        p_ins = (insRate*(insSpace+1.)) / p_event
        p_del = delRate*delSpace / p_event
        
        return p_event, p_ins, p_del
        
            
    def shouldIinsert(self, p_ins, p_del):
        ''' random number draw to determine if insertion or deletion will occur 
        '''
        r = rn.uniform(0,1)
        if r <= p_ins:
            return True
        else:
            return False
    
          
    def generateInsertion(self, model, numsites, nodeID):
        ''' 1. Get a length for the insertion
            2. Generate the inserted sequence using base frequencies
            3. Decide on a location for the inserted sequence
        '''
        # Generate a length
        length = self.generateUnifProb(model.indelParams['insDist']) + 1 # +1 since lengths aren't indexed from zero, and indel length distributions begin at 1.

        # Generated new inserted sites and hold in list called insertion
        freqs = model.substParams['stateFreqs']
        insertion = []
        for j in range(length):
            newsite = generateSite(nodeID, 1)
            newsite.intSeq = self.generateUnifProb(freqs)
            insertion.append(newsite)
            
        # Generate a location
        location = rn.uniform(0, numsites)
        
        return location, insertion
        
        
        
        
    
    def generateDeletion(self, model, deletableSites):
        ''' 1. generate length for deletion
            2. generate starting position for deletion, conditioned on length
        '''
         # Generate a length
        length = self.generateUnifProb(model.indelParams['delDist']) + 1 # +1 since lengths aren't indexed from zero, and indel length distributions begin at 1.

        # Generate a location
        location = rn.randint(0, len(deletableSites) - length)
        return location, length
        
        
    def updateDeletable(self, deletableSites, seq):
        ''' edit list which contains sites that are ok to delete. this needs to be updated as indels are simulated.
        '''
        for n in range(len(seq)):
            if seq[n].intSeq != -1
                deletableSites.append(n)
        return deleteableSites
        
            
    def simIndel(self, parentSeq, nodeID, model, bl, waitTime)
        ''' ABOVE ARGUMENTS NEED OVERHAUL.
            Go, indels, go go go!! Welcome to Gillespieeeeeeeeee!!!! Wooahhhooahhh hey indels! Yeah, yeah, yeah! YEAH YEAH YEAH!     
        '''
        # Tracking variables. Not for whole simulation, just for this branch.
        insertedSites = [] # Will contain list of tuples of len=2. In each tuple, first entry are positions that were inserted for a given insertion event, and second entry is their remaining time (to later use for p=e^(Qt) for those positions).
        deletableSites = [] # Make sure that we don't delete gaps
        deletableSites = self.updateDeletable(deletableSites, parentSeq.intSeq)
    
        # Scale insertion and deletion by branchLength
        insRate = model.indelParams['insRate'] * bl
        delRate = model.indelParams['delRate'] * bl

  
        # Generate event probabilities and waiting time (if the latter doesn't already exist from parent branch) 
        p_event, p_ins, p_del = self.calcProbIndelEvent(seqlen, numsites, insRate, delRate, model.indelParams['meanDelLen'])
        if firstWaitTime is None:
            waitTime = np.random.exponential(scale = 1./p_event), p_event            
        else:
            waitTime = firstWaitTime
        
        ## SIMULATE ##
        remainingTime = branchLength - waitTime
        while remainingTime >= self.zero:
        
            # Perform insertion. track it and add it to sequence.
            if self.shouldIinsert(p_ins, p_del):
                insLocation, insertion = self.generateInsertion( model, numsites, nodeID )
                insertedSites.append( (range(insLocation, len(insertion)), remainingTime) ) # We can always update this in case the inserted sequences are deleted, but I suspect we don't have to because the simSubst function will skip these sites.
                for entry in insertion:
                    partSites.insert(insLocation, entry)
                partSites.insert(insLocation, insertion) 
              
            # Perform deletion and edit relevant site states
            else:
                delLength, delLocation = self.generateDeletion( model, deletableSites )
                deleteme = deletableSites[delLocation: delLocation + delLength] # indices of partSites that we need to switch to deletions
                for x in deleteme:
                    # Change intSeq attr to None (gap). Change state attr from 0->2 (core to deleted core) or from 1->3 (insertion to deleted insertion)
                    assert(partSites[x].intSeq is not None), "You've just deleted a gap! WOAH NOW!"
                    assert(partSites[x].state == 1 or partSites[x].state == 3), "intSeq is not a deletion, but its state is a deletion. WOAH NOW AGAIN!"
                    partSites[x].intSeq = None
                    partSites[x].state += 2
            
            # update deletableSites list
            deletableSites = self.updateDeletable(deletableSites, partSites)          
         
            # next leap
            waitTime = self.generateWaitTime( len(partSites), insRate, delRate)
            remainingTime -= waitTime

    return partSites, insertedSites, abs(remainingTime)
            
            
       