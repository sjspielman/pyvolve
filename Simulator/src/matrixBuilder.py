import numpy as np
from misc import Genetics, Model


class MatrixBuilder(object):
    def __init__(self, model):
        
        # Need to be provided by user
        self.params = model.params
        self.molecules = Genetics()
        self.zero  = 1e-10
        
        ### COMMENTING OUT FOR NOW, BUT MUST RETURN TO THIS!!!!!!!! #####
        ##### THIS CONDITION ACTUALLY IS OK - INVARIANTS ARE PERMITTED, BUT THEY DON'T REALLY HAVE A MATRIX. THEREFORE THEY NEED TO BE SOMEHOW CODED AS INVARIANT. ########
        ## IN CONCLUSION, I WILL NEED TO WORK THIS CONDITION IN MORE CAREFULLY. 
        # Double check that stateFreqs is ok (more than 1 character). 
        #for entry in self.stateFreqs:
        #    assert (1. - entry > self.zero), "You must permit evolution to occur!! Can't only allow one character at a site."    
    
        
    def isTI(self, source, target):
        ''' Returns True for transition, False for transversion.
            Used for nucleotide and codon models.
        '''
        tiPyrim = source in self.molecules.pyrims and target in self.molecules.pyrims
        tiPurine = source in self.molecules.purines and target in self.molecules.purines    
        if tiPyrim or tiPurine:
            return True
        else:
            return False
    
    def getCodonFreq( self, codon):
        ''' Get the frequency for a given codon. 
            Used for codon and mutation-selection models.
        '''
        assert (len(self.params['stateFreqs']) == 61), "This function should ONLY be used when codon frequencies are of interest."
        return self.params['stateFreqs'][self.molecules.codons.index(codon)]
    
    
    def orderNucleotidePair(self, nuc1, nuc2):
        ''' Alphabetize a pair of nucleotides to easily determine reversible mutation rate. 
            The string AG should remain AG, but the string GA should become AG, etc.
            Used for nucleotide, codon, and mutation-selection models.
        '''
        return ''.join(sorted(nuc1+nuc2))
    
    def getNucleotideDiff(self, sourceCodon, targetCodon):
        ''' Get the the nucleotide difference between two codons.
            Returns a string giving sourcenucleotide, targetnucleotide. 
            If this string has 2 characters, then only a single change separates the codons. 
        '''
        nucDiff = ''
        for i in range(3):
            if sourceCodon[i] == targetCodon[i]:    
                continue
            else:    
                nucDiff+=sourceCodon[i]+targetCodon[i]
        if len(nucDiff) != 2:
            return False
        else:
            return nucDiff
        
        
    def buildQ(self):
        ''' Builds instantaneous matrix, Q. 
            For nucleotides, self.size = 4. Amino acids, self.size = 20. Codons, self.size = 61.
        '''    
        self.instMatrix = np.ones( [self.size, self.size] )
        for s in range(self.size):
            for t in range(self.size):
                rate = self.calcInstProb( s, t )                
                self.instMatrix[s][t] = rate
                
            # Fill in the diagonal position so the row sums to 0.
            if np.sum(self.instMatrix[s]) > self.zero: # This check ensures that there are no -0 values in the matrix.
                self.instMatrix[s][s]= -1 * np.sum( self.instMatrix[s] )
            assert ( np.sum(self.instMatrix[s]) < self.zero ), "Row in matrix does not sum to 0."
        self.scaleMatrix()
        return self.instMatrix
        
        
        
    def scaleMatrix(self):
        ''' Scale the instantaneous matrix Q so -Sum(pi_iQ_ii)=1. Ensures branch lengths meaningful for evolving. '''
        scaleFactor = 0
        for i in range(self.size):
            scaleFactor += ( self.instMatrix[i][i] * self.params['stateFreqs'][i] ) ##### IS THIS OK FOR EMPIRICAL MODELS? CHECK THIS!!!
        scaleFactor*=-1.
        self.instMatrix = np.divide( self.instMatrix, scaleFactor )
        ######## CHECK THAT THE SCALING WORKED OUT ##############
        sum=0.
        for i in range(self.size):
            sum += ( self.instMatrix[i][i] * self.params['stateFreqs'][i] )
        assert( abs(sum + 1.) <  self.zero ), "Matrix scaling was a bust."
    
    
    
    def calcInstProb(self, source, target):
        ''' BASE CLASS FUNCTION. NOT IMPLEMENTED.
            Children classes primarily use this function to calculate an entry for the instantaneous rate matrix.
        '''



class aminoAcid_MatrixBuilder(MatrixBuilder):
    ''' This class implements functions relevant to constructing amino acid model instantaneous matrices (Q).
        Deals with empirical matrices, which are coded in empiricalMatrices.py.
    '''        
    def __init__(self, model):
        super(aminoAcid_MatrixBuilder, self).__init__(model)
        self.size = 20
        self.code = self.molecules.amino_acids
        self.initEmpiricalMatrix()
        
    def initEmpiricalMatrix(self):
        import empiricalMatrices as em
        try:
            aaModel = self.params['aaModel'].lower() # I have everything coded in lower case
        except KeyError:
            print "Need an empirical model specification, please"
        try:
            self.empMat = eval("em."+aaModel+"_matrix")
        except:
            print "Couldn't figure out your empirical matrix specification. Note that we currently only support the JTT, WAG, or LG empirical amino acid models."
            
    def calcInstProb(self, source, target):
        ''' Simply return s_ij * p_j'''
        return self.empMat[source][target] * self.params['stateFreqs'][target]        
        
        

class nucleotide_MatrixBuilder(MatrixBuilder):
    ''' This class implements functions relevant to constructing nucleotide model instantaneous matrices (Q).
        All models are essentially nested versions of GTR.
    '''        
    def __init__(self, model):
        super(nucleotide_MatrixBuilder, self).__init__(model)
        self.size = 4
        self.code = self.molecules.nucleotides

    def calcInstProb(self, source, target):
        ''' Calculate instantaneous probability for nucleotide substitutions. '''
        sourceNuc = self.code[source]
        targetNuc = self.code[target]
        substProb = self.params['stateFreqs'][target] * self.params['mu'][sourceNuc+targetNuc]
        return substProb


class codon_MatrixBuilder(MatrixBuilder):    
    ''' This class implements functions relevant to constructing codon model instantaneous matrices (Q).
        Note that the GY and MG models are essentially nested versions of MGREV (which could also be MGHKY, really anything), so can include them all here. Nested versions will merely have extraneous variables fixed to 1 (done elsewhere).
        Model citations:
            GY94:      Yang Z. 1998.
            MG94:      Muse SV, Gaut BS. 1994.
            MG94(REV): Kosakovsky Pond SL, Muse SV. 2005.
            ECM:       Kosiol and Goldman, 2007.
    '''        
    
    
    ############ REQUIRES MASSIVE OVERHAUL TO ACCOMODATE THE FACT THAT MG94 USES NUC FREQS
    ############ ADDITIONALLY, NEED TO CODE EMPIRICAL CODON MODELS.
    ######## THIS WILL LIKELY REQUIRE TWO SUBCLASSES: MXNCODON AND EMPCODON ##############
    
    
    def __init__(self, model):
        super(codon_MatrixBuilder, self).__init__(model)
        self.size = 61
        self.code = self.molecules.codons
        
        # We'll need the nucleotide frequencies for this model, not the codon frequencies....
        #if self.params['codonType'] == 'MG94':
            



    def isSyn(self, sourceCodon, targetCodon):
        ''' Returns True for synonymous codon change, False for nonsynonymous codon change.'''
        ''' Input arguments source and target are each three-letter codons. '''
        if ( self.molecules.codon_dict[sourceCodon] == self.molecules.codon_dict[targetCodon] ):
            return True
        else:
            return False
    
    def calcSynProb( self, target, sourceNuc, targetNuc ):
        ''' Calculate the probability of synonymous change to the targetCodon.'''
        nucPair = self.orderNucleotidePair( sourceNuc, targetNuc)

        if len(target) == 3:
            targetFreq = self.getCodonFreq(target)
        elif len(target) == 1:
            targetFreq = self.getNucFreq(target)
            
        return ( self.getCodonFreq( targetCodon ) * self.params['alpha'] * self.params['mu'][nucPair] )
    
    
    def calcNonsynProb( self, targetCodon, sourceNuc, targetNuc ):
        ''' Calculate the probability of nonsynonymous change to the targetCodon. Note that TI/TV issues are taken care of via the mutation parameters. '''
        nucPair = self.orderNucleotidePair(sourceNuc, targetNuc)
        return ( self.getCodonFreq( targetCodon ) * self.params['beta'] * self.params['mu'][nucPair] )


    def calcInstProb(self, source, target):
        ''' Calculate instantaneous probabilities for codon model matrices.    ''' 
        sourceCodon = self.code[source]
        targetCodon = self.code[target]
        nucDiff = self.getNucleotideDiff(sourceCodon, targetCodon)
        if not nucDiff:
            return 0
        else:
            if self.isSyn(sourceCodon, targetCodon):
                return self.calcSynProb(targetCodon, nucDiff[0], nucDiff[1])
            else:
                return self.calcNonsynProb(targetCodon, nucDiff[0], nucDiff[1])



class mutSel_MatrixBuilder(MatrixBuilder):    
    ''' Implements functions relevant to constructing mutation-selection balance model instantaneous matrices (Q).
        Currently, Halpern and Bruno.
    '''
    def __init__(self, model):
        super(mutSel_MatrixBuilder, self).__init__(model)
        self.size = 61
        self.code = self.molecules.codons
        # PARAMETERS: mu (BH model has non-reversible mutation rates, which I think might be a violation, but will code for now), amino acid frequencies/propensities.
        # Kappa can be included, would be incoporated into mu's before reaching here though.    

    def calcSubstitutionProb(self, sourceFreq, targetFreq, mu_forward, mu_backward):
        ''' Given pi(i) and pi(j) and nucleotide mutation rates, where pi() is the equilibrium frequency/propensity of a given codon, return substitution probability.
            Substitution probability = prob(fixation) * forward_mutation_rate.
        '''
        assert (sourceFreq > 0. and targetFreq > 0. and sourceFreq != targetFreq), "calcSubstitutionProb called when should not have been!" 
        numerator = np.log( (targetFreq*mu_forward)/(sourceFreq*mu_backward) )
        denominator = 1 - ( (sourceFreq*mu_backward)/(targetFreq*mu_forward) )    
        fixProb = numerator/denominator
        substProb = fixProb * mu_forward
        return substProb
        
    def calcInstProb(self, source, target):
        ''' Calculate instantaneous probability for source -> target substitution. ''' 
        sourceCodon = self.code[source]
        targetCodon = self.code[target]
        nucDiff = self.getNucleotideDiff(sourceCodon, targetCodon)
        if nucDiff:
            sourceFreq = self.getCodonFreq(sourceCodon)
            targetFreq = self.getCodonFreq(targetCodon)
            if sourceFreq == 0 or targetFreq == 0:
                return 0
            else:            
                nucPair_forward = nucDiff # mu for source -> target
                mu_forward = self.params["mu"][nucPair_forward]
                if sourceFreq == targetFreq:
                    return mu_forward
                else:
                    nucPair_backward = nucPair_forward[1] + nucPair_forward[0]  # mu for target -> source
                    mu_backward = self.params["mu"][nucPair_backward]
                    return self.calcSubstitutionProb(sourceFreq, targetFreq, mu_forward, mu_backward) 
        else:
            return 0























    
    
            