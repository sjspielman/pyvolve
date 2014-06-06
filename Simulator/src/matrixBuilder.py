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
            Source and target are NUCLEOTIDES, NOT INDICES.
        '''
        tiPyrim = source in self.molecules.pyrims and target in self.molecules.pyrims
        tiPurine = source in self.molecules.purines and target in self.molecules.purines    
        if tiPyrim or tiPurine:
            return True
        else:
            return False
            
    def isSyn(self, source, target):
        ''' Returns True for synonymous codon change, False for nonsynonymous codon change.
            Input arguments source and target are codon indices.
        '''
        sourceCodon = self.molecules.codons[source]
        targetCodon = self.molecules.codons[target]
        if ( self.molecules.codon_dict[sourceCodon] == self.molecules.codon_dict[targetCodon] ):
            return True
        else:
            return False           
    
    def orderNucleotidePair(self, nuc1, nuc2):
        ''' Alphabetize a pair of nucleotides to easily determine reversible mutation rate. 
            The string AG should remain AG, but the string GA should become AG, etc.
            Used for nucleotide, codon, and mutation-selection models.
        '''
        return ''.join(sorted(nuc1+nuc2))
    
    
    def getNucleotideDiff(self, source, target, positionFlag = False, multiple = False):
        ''' Get the the nucleotide difference between two codons.
            Input arguments are codon indices.
            multiple - should we check for multiple changes or not? Nearly always False, but not for ECM.
            Returns a string giving sourcenucleotide, targetnucleotide. 
            If this string has 2 characters, then only a single change separates the codons. 
        '''
        sourceCodon = self.molecules.codons[source]
        targetCodon = self.molecules.codons[target]
        position = None # which position in the codon has changed.
        nucDiff = ''
        for i in range(3):
            if sourceCodon[i] == targetCodon[i]:    
                continue
            else:
                position = i    
                nucDiff+=sourceCodon[i]+targetCodon[i]
        
        # Lots of dynamic returns.
        if multiple: # NEVER needs to be positional.
            return nucDiff
        
        if positionFlag: # NEVER needs to be multiple.
            if len(nucDiff) != 2:
                return False, None
            else:
                return nucDiff, position
        # Default behavior.
        else:
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
    def __init__(self, *args):
        super(aminoAcid_MatrixBuilder, self).__init__(*args)
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
            raise AssertionError("\n\nCouldn't figure out your empirical matrix specification. Note that we currently only support the JTT, WAG, or LG empirical amino acid models.")
            
    def calcInstProb(self, source, target):
        ''' Simply return s_ij * p_j'''
        return self.empMat[source][target] * self.params['stateFreqs'][target]        
        
        





class nucleotide_MatrixBuilder(MatrixBuilder):
    ''' This class implements functions relevant to constructing nucleotide model instantaneous matrices (Q).
        All models are essentially nested versions of GTR.
    '''        
    def __init__(self, *args):
        super(nucleotide_MatrixBuilder, self).__init__(*args)
        self.size = 4
        self.code = self.molecules.nucleotides

    def calcInstProb(self, source, target):
        ''' Calculate instantaneous probability for nucleotide substitutions. '''
        sourceNuc = self.code[source]
        targetNuc = self.code[target]
        return self.params['stateFreqs'][target] * self.params['mu'][sourceNuc+targetNuc]







class empCodon_MatrixBuilder(MatrixBuilder):
    ''' This child class implements functions relevant to constructing *empirical* codon model instantaneous matrices (Q).
        Note that these matrices can also include parameters for kappa (k_ti, k_tv, as described in Kosiol2007) and omega (we will use beta, alpha to allow for dS variation). 
        Currently supporting only ECM (6/5/14). 
    ''' 
    def __init__(self, *args):
        super(empCodon_MatrixBuilder, self).__init__(*args)
        self.size = 61
        self.code = self.molecules.codons
        self.initEmpiricalMatrix() # defines attributes self.restricted (bool), self.empMat
        
        
        
    def initEmpiricalMatrix(self):
        ''' Brings in the empirical rate matrix. Similar, but not identical, to function in amino acid child class.
            Here, we can bring in either the restricted (1 inst change only) or unrestricted (1-3 inst changes) matrix.
        '''
        import empiricalMatrices as em
        try:
            self.restricted = self.params['restricted']
            assert(type(self.restricted) is bool), ("\n\nNeed to specify True or False for restricted.")
            if self.restricted:
                self.empMat = em.ecmrest_matrix
            else:
                self.empMat = em.ecmunrest_matrix
        except KeyError:
            raise AssertionError("\n\nMust specify if you want restricted or unrestricted ECM.")



    def setKappaParam(self, nucDiff):
        ''' Calculations for the "kappa" parameter(s) for ECM. See paper for details.'''
        num_ti = 0
        num_tv = 0
        for i in range(0, len(nucDiff), 2):
            sourceNuc = nucDiff[i]
            targetNuc = nucDiff[i+1]
            if self.isTI(sourceNuc, targetNuc):
                num_ti += 1
            else:
                num_tv += 1
        return self.params['k_ti']**num_ti * self.params['k_tv']**num_tv



    def calcInstProb(self, source, target):
        '''look, a description!'''     
        
        nucDiff = self.getNucleotideDiff(source, target, False, True) # no positional needed, but account for multiple
        if len(nucDiff) == 0  or (self.restricted and len(nucDiff) != 2):
            return 0 # no change, or restricted is True and there is more than 1 change.
        else:
            
            kappaParam = self.setKappaParam(nucDiff) # which ones changed? set kappa stuff accordingly
            if self.isSyn(source, target):
                return self.empMat[source][target] * self.params['stateFreqs'][target] * self.params['alpha'] * kappaParam
            else:
                return self.empMat[source][target] * self.params['stateFreqs'][target] * self.params['beta'] * kappaParam















class mechCodon_MatrixBuilder(MatrixBuilder):    
    ''' This child class implements functions relevant to constructing *mechanistic* codon model instantaneous matrices (Q).
        Model citations:
            GY94:      Yang Z. 1994,1998.
            MG94:      Muse SV, Gaut BS. 1994.
            MG94(REV): Kosakovsky Pond SL, Muse SV. 2005.
    '''        
 
 
    def __init__(self, *args):
        super(mechCodon_MatrixBuilder, self).__init__(*args)
        self.size = 61
        self.code = self.molecules.codons
        # Assign self.modelClass to GY94 or MG94 based on state frequencies.
        if self.params['stateFreqs'].shape == (61,):
            self.modelClass = 'GY94'
        else:
            if self.params['stateFreqs'].shape == (3,4):
                self.modelClass = 'MG94'
            # Convert to positional nucleotide frequencies.
            elif self.params['stateFreqs'].shape == (4,):
                self.modelClass = 'MG94'
                self.params['stateFreqs'] = np.array( [self.params['stateFreqs'], self.params['stateFreqs'], self.params['stateFreqs'] ])
            else:
                raise AssertionError("\n\nIf you'd like to use a mechanistic codon model, you must work with either codon state frequencies (GY94 model) or nucleotide state frequencies (MG94 model). \nYou're doing neither, and that's unfortunate.\nI'm quitting now.")


    def getTargetFreq(self, target, position):
        ''' Function to return target frequency for mechanistic codon model.
            Note that the argument target is a **codon** index.
            Note that the arguement position will only come into play for the MG94 model class.
        '''
        if self.modelClass == 'GY94':
            return self.params['stateFreqs'][target]
        else:
            return self.params['stateFreqs'][position, target]
  

    def calcSynProb(self, targetFreq, nucPair):
        ''' Calculate instantaneous probability of synonymous change for mechanistic codon model.'''
        return targetFreq * self.params['alpha'] * self.params['mu'][nucPair]
    
    
    def calcNonsynProb(self, targetFreq, nucPair):
        ''' Calculate instantaneous probability of nonsynonymous change for mechanistic codon model.'''
        return targetFreq * self.params['beta'] * self.params['mu'][nucPair]


    def calcInstProb(self, source, target):
        ''' Calculate instantaneous probabilities for mechanistic codon model matrices.
            Arguments:
                source,target are indices.
        ''' 
        nucDiff, position = self.getNucleotideDiff(source, target, True)
        if not nucDiff:
            return 0
        else:
            nucPair = self.orderNucleotidePair( nucDiff[0], nucDiff[1] )
            if self.modelClass == 'GY94':
                targetFreq = self.getTargetFreq(target, position)
            else:
                targetNuc = self.molecules.nucleotides.index(nucDiff[1])
                targetFreq = self.getTargetFreq(targetNuc, position)
            if self.isSyn(source, target):
                return self.calcSynProb(targetFreq, nucPair)
            else:
                return self.calcNonsynProb(targetFreq, nucPair)









##################    TO DO     ###########################################
#####  1. ALLOW FOR CODON BIAS. DO THIS BY PROVIDING W/ CODON FREQUENCIES INSTEAD. SHOULD BE DYNAMIC BASED ON TYPE OF FREQ PROVIDED
#####  2. IMPLEMENT TAMURI, SELLAHIRSH, RODRIGUE FLAVORS.
#####  3. IMPLEMENT A NUCLEOTIDE-LEVEL MUTATION-SELECTION MODEL
############################################################################

class mutSel_MatrixBuilder(MatrixBuilder):    
    ''' Implements functions relevant to constructing mutation-selection balance model instantaneous matrices (Q).
        Currently, Halpern and Bruno.
    '''
    def __init__(self, *args):
        super(mutSel_MatrixBuilder, self).__init__(*args)
        self.size = 61
        self.code = self.molecules.codons
        # PARAMETERS: mu, amino acid frequencies/propensities.
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
        nucDiff = self.getNucleotideDiff(source, target) # no positional aspects and only single changes allowed.
        if nucDiff:
            sourceFreq = self.params['stateFreqs'][source]
            targetFreq = self.params['stateFreqs'][target]
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























    
    
            