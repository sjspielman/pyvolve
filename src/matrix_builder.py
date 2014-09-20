import numpy as np
from misc import ZERO, Genetics
MOLECULES = Genetics()



class Matrix_Builder(object):
    def __init__(self, model):
        
        self.params = model.params
     


    def _is_TI(self, source, target):
        ''' Returns True for transition, False for transversion.
            Used for nucleotide and codon models.
            Source and target are NUCLEOTIDES, NOT INDICES.
        '''
        ti_pyrim = source in MOLECULES.pyrims and target in MOLECULES.pyrims
        ti_purine = source in MOLECULES.purines and target in MOLECULES.purines    
        if ti_pyrim or ti_purine:
            return True
        else:
            return False
            
    def _is_syn(self, source, target):
        ''' Returns True for synonymous codon change, False for nonsynonymous codon change.
            Input arguments source and target are codon indices.
        '''
        source_codon = MOLECULES.codons[source]
        target_codon = MOLECULES.codons[target]
        if ( MOLECULES.codon_dict[source_codon] == MOLECULES.codon_dict[target_codon] ):
            return True
        else:
            return False           
    
    def _order_nucleotide_pair(self, nuc1, nuc2):
        ''' Alphabetize a pair of nucleotides to easily determine reversible mutation rate. 
            The string AG should remain AG, but the string GA should become AG, etc.
            Used for nucleotide, codon, and mutation-selection models.
        '''
        return ''.join(sorted(nuc1+nuc2))
    
    
    def _get_nucleotide_diff(self, source, target):
        ''' Get the the nucleotide difference between two codons.
            Input arguments are codon indices.
            multiple - should we check for multiple changes or not? Nearly always False, but not for ECM.
            Returns a string giving sourcenucleotide, targetnucleotide. 
            If this string has 2 characters, then only a single change separates the codons. 
        '''
        source_codon = MOLECULES.codons[source]
        target_codon = MOLECULES.codons[target]
        nuc_diff = ''
        for i in range(3):
            if source_codon[i] == target_codon[i]:    
                continue
            else:
                nuc_diff += source_codon[i]+target_codon[i]
        return nuc_diff




    def buildQ(self):
        ''' Builds instantaneous matrix, Q. 
            For nucleotides, self.size = 4. Amino acids, self.size = 20. Codons, self.size = 61.
        '''    
        self.instMatrix = np.zeros( [self.size, self.size] )
        for s in range(self.size):
            for t in range(self.size):
                # Non-diagonal
                rate = self.calcInstProb( s, t )                
                self.instMatrix[s][t] = rate
                
            # Fill in the diagonal position so the row sums to 0, but ensure it doesn't became -0
            self.inst_matrix[s][s]= -1. * np.sum( self.inst_matrix[s] )
            if self.inst_matrix[s][s] == -0.:
                self.inst_matrix[s][s] = 0.
            assert ( abs(1. - np.sum(self.inst_matrix[s])) < ZERO ), "Row in instantaneous matrix does not sum to 0."
        self._scale_matrix()
        return self.inst_matrix
        
        
    def _scale_matrix(self):
        ''' Scale the instantaneous matrix Q so -Sum(pi_iQ_ii)=1. Ensures branch lengths meaningful for evolving. '''
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
        ''' BASE CLASS FUNCTION. NOT IMPLEMENTED.
            Child classes primarily use this function to calculate an entry for the instantaneous rate matrix.
        '''



class aminoAcid_Matrix(Matrix_Builder):
    ''' This class implements functions relevant to constructing amino acid model instantaneous matrices (Q).
        Deals with empirical matrices, which are coded in empiricalMatrices.py.
    '''        
    def __init__(self, *args):
        super(aminoAcid_Matrix, self).__init__(*args)
        self._size = 20
        self._code = MOLECULES.amino_acids
        self._init_empirical_matrix()
        
    def _init_empirical_matrix(self):
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
        ''' Simply return s_ij * p_j'''
        return self.emp_matrix[source][target] * self.params['state_freqs'][target]        






class nucleotide_Matrix(Matrix_Builder):
    ''' This class implements functions relevant to constructing nucleotide model instantaneous matrices (Q).
        All models are essentially nested versions of GTR.
    '''        
    def __init__(self, *args):
        super(nucleotide_Matrix, self).__init__(*args)
        self._size = 4
        self._code = MOLECULES.nucleotides

    def _calc_instantaneous_prob(self, source, target):
        ''' Calculate instantaneous probability for nucleotide substitutions. '''
        source_nuc = self.code[source]
        target_nuc = self.code[target]
        if source_nuc == target_nuc:
            return 0.
        else:
            return self.params['state_freqs'][target] * self.params['mu'][self._order_nucleotide_pair(source_nuc,target_nuc)]







class empiricalCodon_Matrix(Matrix_Builder):
    ''' This child class implements functions relevant to constructing *empirical* codon model instantaneous matrices (Q).
        Note that these matrices can also include parameters for kappa (k_ti, k_tv, as described in Kosiol2007) and omega (we will use beta, alpha to allow for dS variation). 
        Currently supporting only ECM (6/5/14). 
    ''' 
    def __init__(self, *args):
        super(empiricalCodon_Matrix, self).__init__(*args)
        self._size = 61
        self._code = MOLECULES.codons
        self._init_empirical_matrix() # defines attributes self.restricted (bool), self.empMat
        
        
        
   def _init_empirical_matrix(self):
        ''' Brings in the empirical rate matrix. Similar, but not identical, to function in amino acid child class.
            Here, we can bring in either the restricted (1 inst change only) or unrestricted (1-3 inst changes) matrix.
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
        ''' Calculations for the "kappa" parameter(s) for ECM. See paper for details.'''
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
        '''look, a description!'''     
        
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) == 0  or (self.restricted and len(nuc_diff) != 2):
            return 0.
        else:
            
            kappa_param = self._set_kappa_param(nuc_diff) # which ones changed? set kappa stuff accordingly
            if self._is_syn(source, target):
                return self.emp_matrix[source][target] * self.params['state_freqs'][target] * self.params['alpha'] * kappa_param
            else:
                return self.emp_matrix[source][target] * self.params['state_freqs'][target] * self.params['beta'] * kappa_param










class codonGY_Matrix(Matrix_Builder):    
    ''' This child class implements functions relevant to constructing GY94-style codon model instantaneous matrices (Q).
        By GY94-style, I mean target codon frequencies are used in the matrix. Note that this class can still accomodate distinct nonsyn and syn rates and arbitrary mutational parameters.
    '''        
 
 
    def __init__(self, *args):
        super(codonGY_Matrix, self).__init__(*args)
        self._size = 61
        self._code = MOLECULES.codons
       

    def _calc_syn_prob(self, target, nuc_pair):
        ''' Calculate instantaneous probability of synonymous change for mechanistic codon model.'''
        return self.params['state_freqs'][target] * self.params['alpha'] * self.params['mu'][nuc_pair]
    
    
    def _calc_nonsyn_prob(self, target, nuc_pair):
        ''' Calculate instantaneous probability of nonsynonymous change for mechanistic codon model.'''
        return self.params['state_freqs'][target] * self.params['beta'] * self.params['mu'][nuc_pair]


    def _calc_instantaneous_prob(self, source, target):
        ''' Calculate instantaneous probabilities for mechanistic codon model matrices.
            source,target are codon indices.
        ''' 
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) != 2:
            return 0.
        else:
            nuc_pair = self._order_nucleotide_pair( nuc_diff[0], nuc_diff[1] )
            if self._is_syn(source, target):
                return self._calc_syn_prob(target, nuc_pair)
            else:
                return self._calc_nonsyn_prob(target, nuc_pair)




class mutSel_Matrix(Matrix_Builder):    
    ''' Implements functions relevant to constructing mutation-selection balance model instantaneous matrices (Q).
        The codon mutation-selection model implemented is that described in Halpern Bruno 1998.
    '''
    
    def __init__(self, *args):
        super(mutSel_Matrix, self).__init__(*args)
        
         # Assign self._model_class to codon or nucleotide based on state frequencies.
        if self.subst_params['state_freqs'].shape == (61,):
            self._model_class = 'codon'
            self._size = 61
            self._code = MOLECULES.codons
        elif self.subst_params['state_freqs'].shape == (4,):
            self._model_class = 'nuc'
            self._size = 4
            self._code = MOLECULES.nucleotides
        else:
            raise AssertionError("\n\nMutSel models need either codon or nucleotide frequencies.")

           

    def _calc_substitution_prob(self, pi_i, pi_j, mu_ij, mu_ji):
        ''' Given pi(i) and pi(j) and nucleotide mutation rates, where pi() is the equilibrium frequency of a given codon, return substitution probability.
            Substitution probability = prob(fixation) * forward_mutation_rate.
        '''

        if pi_i <= ZERO or pi_j <= ZERO: 
            return 0.
        elif abs(pi_i - pi_j) <= ZERO:
            return mu_ij
        else:
            pi_mu = (pi_j*mu_ji)/(pi_i*mu_ij)
            return np.log(pi_mu)/(1. - 1./pi_mu) * mu_ij
    


    def _calc_instantaneous_prob(self, source, target):
        ''' Calculate instantaneous probability for source -> target substitution. ''' 
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) != 2:
            return 0.
        else:
            source_freq = self.params['state_freqs'][source]
            target_freq = self.params['state_freqs'][target]
            mu_forward = self.params["mu"][nuc_diff]
            mu_backward = self.params["mu"][nuc_diff[1] + nuc_diff[0]]
            
            return self._calc_substitution_prob(source_freq, target_freq, mu_forward, mu_backward) 

            
            
            
            
            
            
            
            
            
