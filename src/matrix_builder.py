#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

### TO DO: INCORPORATE MG INTO MECHCODON_MATRIX CHILD CLASS ###


'''
Generate the instantaneous rate matrix for Markov chain.
'''



import numpy as np
from misc import ZERO, Genetics
MOLECULES = Genetics()



class MatrixBuilder(object):
    '''
        Parent class for model instantaneous matrix generation.
        
        Child class include the following:
        1. *aminoAcid_Matrix*       : Specific functionality for empirical amino acid models (currently, JTT, WAG, and LG and SCG05).
        2. *nucleotide_Matrix*      : Specific functionality for nucleotide models (GTR and nested).
        3. *mechCodon_Matrix*       : Specific functionality for so-called mechanistic codon models, which include GY-style and MG-style models (dN/dS models)
        5. *ECM_Matrix*             : Specific functionality for the ECM (Kosiol2007) empirical codon model
        6. *mutSel_Matrix*          : Specific functionality for the mutation-selection-balance model (Halpern and Bruno 1998). Extended to work for either codon or nucleotides.


        TODO:
            Incorpprate SCG05 empirical codon model. 
    '''
    
    def __init__(self, model):
        
        self.params = model.params
     


    def _is_TI(self, source, target):
        ''' 
            Determine if a given nucleotide change is a transition or a tranversion.  Used in child classes nucleotide_Matrix, mechCodon_Matrix, ECM_Matrix, mutSel_Matrix .
            Returns True for transition, False for transversion.
            
            Arguments "source" and "target" are the actual nucleotides (not indices).
        '''
        
        ti_pyrim = source in MOLECULES.pyrims and target in MOLECULES.pyrims
        ti_purine = source in MOLECULES.purines and target in MOLECULES.purines    
        if ti_pyrim or ti_purine:
            return True
        else:
            return False
            
    def _is_syn(self, source, target):
        '''
            Determine if a given codon change is synonymous or nonsynonymous. Used in child classes mechCodon_Matrix, ECM_Matrix .
            Returns True for synonymous, False for nonsynonymous.
             
            Arguments arguments "source" and "target" are codon indices (0-60, alphabetical).
        '''
        
        source_codon = MOLECULES.codons[source]
        target_codon = MOLECULES.codons[target]
        if ( MOLECULES.codon_dict[source_codon] == MOLECULES.codon_dict[target_codon] ):
            return True
        else:
            return False           
    
    
    
    def _get_nucleotide_diff(self, source, target):
        ''' 
            Get the nucleotide difference(s) between two codons. Used in child classes ECM_Matrix, mechCodon_Matrix, mutSel_Matrix .
            Returns a string representing the nucleotide differences between source and target codon.
            For instance, if source is AAA and target is ATA, the string AT would be returned. If source is AAA and target is ACC, then ACAC would be returned.
            
            Input arguments source and target are codon indices (0-60, alphabetical).
        '''        
        
        source_codon = MOLECULES.codons[source]
        target_codon = MOLECULES.codons[target]
        return "".join( [source_codon[i]+target_codon[i] for i in range(len(source_codon)) if source_codon[i] != target_codon[i]] )




    def buildQ(self):
        ''' 
            General function to generate the instantaneous matrix, Q. Used with all child classes.
        '''    
        
        self.inst_matrix = np.zeros( [self._size, self._size] ) # For nucleotides, self._size = 4; amino acids, self._size = 20; codons, self._size = 61.
        for s in range(self._size):
            for t in range(self._size):
                # Non-diagonal
                rate = self._calc_instantaneous_prob( s, t )                
                self.inst_matrix[s][t] = rate
                
            # Fill in the diagonal position so the row sums to 0, but ensure it doesn't became -0
            self.inst_matrix[s][s]= -1. * np.sum( self.inst_matrix[s] )
            if self.inst_matrix[s][s] == -0.:
                self.inst_matrix[s][s] = 0.
            assert ( abs(np.sum(self.inst_matrix[s])) < ZERO ), "Row in instantaneous matrix does not sum to 0."
        self._scale_matrix()
        return self.inst_matrix
        
        
    def _scale_matrix(self):
        ''' 
            Scale the instantaneous matrix Q so -Sum(pi_iq_ii)=1, where q_ii is diagonal matrix element in row/column i and p_i is the state frequency for i.
            This scaling ensures branch lengths meaningful for evolving. 
        '''
        
        scaling_factor = 0.
        for i in range(self._size):
            scaling_factor += ( self.inst_matrix[i][i] * self.params['state_freqs'][i] )
        self.inst_matrix /= -1.*scaling_factor
        ## Double check that scaling worked
        sum = 0.
        for i in range(self._size):
            assert ( -1.* ZERO < np.sum(self.inst_matrix[i]) < ZERO ), "After scaling, row in matrix does not sum to 0."
            sum += ( self.inst_matrix[i][i] * self.params['state_freqs'][i] )
        assert( abs(sum + 1.) <  ZERO ), "Matrix scaling was a bust."
    
    
    
    def _calc_instantaneous_prob(self, source, target):
        ''' 
            BASE CLASS FUNCTION. NOT IMPLEMENTED.
            Child classes use this function to calculate a given element in the instantaneous rate matrix.
            Returns the substitution probability from source to target, for a given model.
            
            Arguments "source" and "target" are indices for the relevant aminos (0-19) /nucs (0-3) /codons (0-60). 
        '''
        return 0



class aminoAcid_Matrix(MatrixBuilder):
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing amino acid model instantaneous matrices.
        Note that all empirical amino acid replacement matrices are in the file empirical_matrices.py.
    '''        
    
    def __init__(self, *args):
        super(aminoAcid_Matrix, self).__init__(*args)
        self._size = 20
        self._code = MOLECULES.amino_acids
        self._init_empirical_matrix()
        
        
    def _init_empirical_matrix(self):
        '''
            Function to load the appropriate replacement matrix from empirical_matrices.py 
        '''
        import empirical_matrices as em
        try:
            aa_model = self.params['aa_model'].lower() # I have everything coded in lower case
        except KeyError:
            print "Need an empirical model specification, please"
        try:
            self.emp_matrix = eval("em."+aa_model+"_matrix")
        except:
            raise AssertionError("\n\nCouldn't figure out your empirical matrix specification. Note that we currently only support the JTT, WAG, or LG empirical amino acid models.")
            
    def _calc_instantaneous_prob(self, source, target):
        ''' 
            Returns the substitution probability (s_ij * p_j, where s_ij is replacement matrix entry and p_j is target amino frequency) from source to target for amino acid empirical models.
            Arguments "source" and "target" are indices for the relevant aminos (0-19).
        '''
        return self.emp_matrix[source][target] * self.params['state_freqs'][target]        






class nucleotide_Matrix(MatrixBuilder):
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing nucleotide model instantaneous matrices (Q).
        All models computed here are essentially nested versions of GTR.
    '''        
    
    def __init__(self, *args):
        super(nucleotide_Matrix, self).__init__(*args)
        self._size = 4
        self._code = MOLECULES.nucleotides

    def _calc_instantaneous_prob(self, source, target):
        ''' 
            Returns the substitution probability (\mu_ij * p_j, where \mu_ij are nucleotide mutation rates and p_j is target nucleotide frequency) from source to target for nucleotide models.
            Arguments "source" and "target" are indices for the relevant nucleotide (0-3).
        '''
        source_nuc = self._code[source]
        target_nuc = self._code[target]
        if source_nuc == target_nuc:
            return 0.
        else:
            return self.params['state_freqs'][target] * self.params['mu']["".join(sorted(source_nuc + target_nuc))]







class ECM_Matrix(MatrixBuilder):
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing a matrix specifically for the ECM (described in Kosiol2007) model.
        We support both restricted (instaneous single changes only) and unrestricted (instantaneous single, double, or triple) versions of this model (see paper for details).
        
        !!! NOTE: The ECM model supports omega (dN/dS) and kappa (TI/TV) ratios in their calculations, and therefore I have included these parameters here. HOWEVER, I do NOT recommend their use.
    
    ''' 
    
    def __init__(self, *args):
        super(ECM_Matrix, self).__init__(*args)
        self._size = 61
        self._code = MOLECULES.codons
        self._init_empirical_matrix() # defines attributes self.restricted (bool), self.empMat
        
        
        
    def _init_empirical_matrix(self):
        '''
            Function to load the appropriate replacement matrix from empirical_matrices.py 
        '''
        import empirical_matrices as em
        try:
            self.restricted = self.params['restricted']
            assert(type(self.restricted) is bool), ("\n\nNeed to specify True or False for restricted.")
            if self.restricted:
                self.emp_matrix = em.ecmrest_matrix
            else:
                self.emp_matrix = em.ecmunrest_matrix
        except KeyError:
            raise AssertionError("\n\nMust specify if you want restricted or unrestricted ECM.")



    def _set_kappa_param(self, nuc_diff):
        ''' 
            Calculations for the "kappa" parameter(s) for the ECM model. See the 2007 paper for details.
        '''
        num_ti = 0
        num_tv = 0
        for i in range(0, len(nuc_diff), 2):
            source_nuc = nuc_diff[i]
            target_nuc = nuc_diff[i+1]
            if self._is_TI(source_nuc, target_nuc):
                num_ti += 1
            else:
                num_tv += 1
        return self.params['k_ti']**num_ti * self.params['k_tv']**num_tv




    def _calc_instantaneous_prob(self, source, target):
        ''' 
            Returns the substitution probability from source to target for ECM models.
            Arguments "source" and "target" are indices for the relevant codons (0-60).
        '''  
        
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) == 0  or (self.restricted and len(nuc_diff) != 2):
            return 0.
        else:
            kappa_param = self._set_kappa_param(nuc_diff)
            if self._is_syn(source, target):
                return self.emp_matrix[source][target] * self.params['state_freqs'][target] * self.params['alpha'] * kappa_param
            else:
                return self.emp_matrix[source][target] * self.params['state_freqs'][target] * self.params['beta'] * kappa_param










class mechCodon_Matrix(MatrixBuilder):    
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to "mechanistic" (dN/dS) codon models.
        Both GY-style (target codon frequencies) and MG-style (target nucleotide frequencies) are accepted.
        Both dS and dN variation are allowed, as are GTR mutational parameters (not strictly HKY85).
  
        TO DO : ACTUALLY INCORPORATE MG-STYLE MODELS. WE NEED TO GET THE NUCLEOTIDE FREQUENCIES FROM THE CODON FREQUENCIES AND USE THEM ACCORDINGLY IN THIS CIRCUMSTANCE !!!
  
    '''        
 
 
    def __init__(self, *args):
        super(mechCodon_Matrix, self).__init__(*args)
        self._size = 61
        self._code = MOLECULES.codons
        
        # These next lines of code will eventually be incorporated elsewhere when model parameter defaults are set up.
        if 'omega' in self.params.keys():
            self.params['beta'] = self.params['omega']
            self.params['alpha'] = 1.0
        
#        self._model_type = self.params['model_type'] # This can be either GY or MG. If MG, we'll need to do some extra computations as target nucleotide are used in the matrix.
#        assert(self._model_type == 'GY' or self._model_type == 'MG'), "\n\nFor mechanistic codon models, you must specify a model_type as GY (uses target *codon* frequencies) or MG (uses target *nucleotide* frequencies.) I RECOMMEND MG!!"
                
              

    def _calc_syn_prob(self, target, nuc_pair):
        ''' 
            Calculate instantaneous probability of synonymous change for mechanistic codon models.
        '''
        return self.params['state_freqs'][target] * self.params['alpha'] * self.params['mu'][nuc_pair]
    
    
    def _calc_nonsyn_prob(self, target, nuc_pair):
        ''' 
            Calculate instantaneous probability of nonsynonymous change for mechanistic codon models.
        '''
        return self.params['state_freqs'][target] * self.params['beta'] * self.params['mu'][nuc_pair]


    def _calc_instantaneous_prob(self, source, target):
        ''' 
            Returns the substitution probability from source to target for mechanistic codon models.
            Arguments "source" and "target" are indices for the relevant codons (0-60).
        ''' 
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) != 2:
            return 0.
        else:
            nuc_pair = "".join(sorted(nuc_diff[0] + nuc_diff[1]))
            if self._is_syn(source, target):
                return self._calc_syn_prob(target, nuc_pair)
            else:
                return self._calc_nonsyn_prob(target, nuc_pair)




class mutSel_Matrix(MatrixBuilder):    
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing mutation-selection balance model instantaneous matrices, according to the HalpernBruno 1998 model.
        Here, this model is extended such that it can be used for either nucleotide or codon. This class will automatically detect which one you want based on your state frequencies.

    '''
    
    def __init__(self, *args):
        super(mutSel_Matrix, self).__init__(*args)
        
         # Assign self._model_class to codon or nucleotide based on state frequencies.
        if self.params['state_freqs'].shape == (61,):
            self._model_class = 'codon'
            self._size = 61
            self._code = MOLECULES.codons
        elif self.params['state_freqs'].shape == (4,):
            self._model_class = 'nuc'
            self._size = 4
            self._code = MOLECULES.nucleotides
        else:
            raise AssertionError("\n\nMutSel models need either codon or nucleotide frequencies.")

           

    def _calc_instantaneous_prob(self, source, target):
        ''' 
            Calculate the substitution probability from source to target for mutation-selection-balance models.
            Arguments "source" and "target" are indices for the relevant codons (0-60) or nucleotide (0-3).
        '''        
        
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) != 2:
            return 0.
        else:
            pi_i  = self.params['state_freqs'][source]           # source frequency
            pi_j  = self.params['state_freqs'][target]           # target frequency 
            mu_ij = self.params["mu"][nuc_diff]                  # source -> target mutation rate
            mu_ji = self.params["mu"][nuc_diff[1] + nuc_diff[0]] # target -> source mutation rate

            if pi_i <= ZERO or pi_j <= ZERO: 
                inst_prob = 0.
            elif abs(pi_i - pi_j) <= ZERO:
                inst_prob = mu_ij
            else:
                pi_mu = (pi_j*mu_ji)/(pi_i*mu_ij)
                inst_prob =  np.log(pi_mu)/(1. - 1./pi_mu) * mu_ij
            return inst_prob
            
            
            
            
            
            
            
