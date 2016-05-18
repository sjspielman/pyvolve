#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################




'''
    This module will generate the Markov chain's instantaneous rate matrix, Q.
'''



import numpy as np
from scipy import linalg
from copy import deepcopy
from .genetics import *
from .state_freqs import *
ZERO      = 1e-8
MOLECULES = Genetics()


class MatrixBuilder(object):
    '''
        Parent class for model instantaneous matrix creation.
        
        Child class include the following:
            1. *AminoAcid_Matrix* 
                - Empirical amino acid models
            2. *Nucleotide_Matrix* 
                - Nucleotide models (GTR and nested)
            3. *MechCodon_Matrix*
                - So-called mechanistic codon models, which include GY-style and MG-style models (dN/dS models)
            4. *ECM_Matrix*
                - ECM (Kosiol2007) empirical codon model
            5. *MutSel_Matrix* 
                - Mutation-selection model (Halpern and Bruno 1998), extended for either codon or nucleotides
        
    '''
    
    def __init__(self, model_type, parameters):
        '''
            Requires two positional argument:
                1. **model_type**, the type of model which will be built
                2. **parameters**, a dictionary containing parameters about the substitution process which will be checked.
        '''
        self.model_type = model_type.lower()
        self.params = parameters  
         

    def _build_matrix( self, parameters = None):
        ''' 
            Generate an instantaneous rate matrix.
        '''    
        if parameters is None:
            parameters = self.params
        matrix = np.zeros( [self._size, self._size] ) # For nucleotides, self._size = 4; amino acids, self._size = 20; codons, self._size = 61.
        for s in range(self._size):
            for t in range(self._size):
                # Non-diagonal
                rate = self._calc_instantaneous_prob( s, t, parameters )                
                matrix[s][t] = rate
                
            # Fill in the diagonal position so the row sums to 0, but ensure it doesn't became -0
            matrix[s][s]= -1. * np.sum( matrix[s] )
            if matrix[s][s] == -0.:
                matrix[s][s] = 0.
            assert ( abs(np.sum(matrix[s])) < ZERO ), "\n\nRow in instantaneous matrix does not sum to 0."
        return matrix




    def _scaling_factor_from_matrix(self, frequencies, matrix):
        '''
            Determine a scaling factor from a given frequency distribution and corresponding un-normalized rate matrix.
        '''
        scaling_factor = 0.
        for i in range(self._size):
            scaling_factor += ( matrix[i][i] * frequencies[i] )
        return scaling_factor

          

    def _obtain_scaling_factor(self):
        '''
            Calculate the appropriate matrix scaling factor.
        '''
        if self.scale_matrix == "persite":
            factor = self._scaling_factor_from_matrix(self.params["state_freqs"], self.inst_matrix)
        else:
            frequencies, matrix = self._build_scaling_matrix()
            factor = self._scaling_factor_from_matrix(frequencies, matrix)
        
        return factor            




    def _build_scaling_matrix(self):  
        '''
            Build the matrix used to determine a scaling factor for either average or neutral scaling.
        '''
        params = self._build_scaling_params()
        matrix = self._build_matrix(parameters = params)
        return params["state_freqs"], matrix


 
    def __call__(self):
        ''' 
            Generate, scale, and return instantaneous rate matrix.
        '''    
        
        # Construct matrix
        self.inst_matrix = self._build_matrix()

        # Scale matrix
        scaling_factor = self._obtain_scaling_factor()
        self.inst_matrix /= -1.*scaling_factor

        return self.inst_matrix



    def _init_empirical_matrix(self, name):
        '''
            Function to load the appropriate replacement matrix from empirical_matrices.py 
        '''
        from . import empirical_matrices as em
        try:
            self.emp_matrix = eval("em."+name+"_matrix")
        except:
            raise ValueError("\n\nCouldn't figure out your empirical matrix specification.")




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
        
        
    
    def _calc_instantaneous_prob(self, source, target, params = None):
        ''' 
            Calculate a given element in the instantaneous rate matrix.
            Returns the substitution probability from source to target, for a given model.
            Arguments "source" and "target" are *indices* for the relevant aminos (0-19) /nucs (0-3) /codons (0-60). 
        '''
        print("Parent class function. Not called.")










class AminoAcid_Matrix(MatrixBuilder):
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing amino acid model instantaneous matrices.
        Note that all empirical amino acid replacement matrices are in the file empirical_matrices.py.
    '''        
    
    def __init__(self, *args):
        super(AminoAcid_Matrix, self).__init__(*args)
        self._size = 20
        self.code = MOLECULES.amino_acids
        self.scale_matrix = "persite"
        self._init_empirical_matrix(self.model_type)
   
            
            
    def _calc_instantaneous_prob( self, source, target, parameters = None ):
        ''' 
            Returns the substitution probability (s_ij * p_j, where s_ij is replacement matrix entry and p_j is target amino frequency) from source to target for amino acid empirical models.
            Arguments "source" and "target" are indices for the relevant aminos (0-19).
        '''
        if parameters is None:
            parameters = self.params
        return self.emp_matrix[source][target] * parameters['state_freqs'][target]        







class Nucleotide_Matrix(MatrixBuilder):
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing nucleotide model instantaneous matrices.
        All models computed here are essentially nested versions of GTR.
    '''        
    
    def __init__(self, *args):
        super(Nucleotide_Matrix, self).__init__(*args)
        self._size = 4
        self._code = MOLECULES.nucleotides
        self.scale_matrix = "persite"


    def _calc_instantaneous_prob(self, source, target, parameters = None):
        ''' 
            Returns the substitution probability (\mu_ij * p_j, where \mu_ij are nucleotide mutation rates and p_j is target nucleotide frequency) from source to target for nucleotide models.
            Arguments "source" and "target" are indices for the relevant nucleotide (0-3).
        '''
        if parameters is None:
            parameters = self.params
        source_nuc = self._code[source]
        target_nuc = self._code[target]
        if source_nuc == target_nuc:
            return 0.
        else:
            return parameters['state_freqs'][target] * parameters['mu']["".join(sorted(source_nuc + target_nuc))]









class MechCodon_Matrix(MatrixBuilder):    
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to "mechanistic" (dN/dS) codon models.
        Models include both GY-style or MG-style varieties, although users should *always specify codon frequencies* to class instance!
        Both dS and dN variation are allowed, as are GTR mutational parameters (not strictly HKY85).
    
    '''        

    def __init__(self, *args):
        super(MechCodon_Matrix, self).__init__(*args)
        self._size = 61
        self._code = MOLECULES.codons
        if "neutral_scaling" not in self.params:
            self.params["neutral_scaling"] = False    
       
        # First check neutral scaling
        if self.params["neutral_scaling"] is True:
            self.scale_matrix = "neutral"
        else:
            # This key will be in the dictionary if it is a heterogeneous codon model, for which we want average scaling. Without heterogeneity, persite is equivalent.
            if "hetmodel_mean_dnds" in self.params:
                self.scale_matrix = "average"
            else:
                self.scale_matrix = "persite"


    def _calc_prob(self, target_codon, target_nuc, nuc_pair, subrate):
        ''' 
            Calculate instantaneous probability of (non)synonymous change for mechanistic codon models.
            Argument *subrate* is either dN or dS.
        '''

        prob =  self.params['mu'][nuc_pair] * subrate 
        if self.model_type == 'gy':
            prob *= self.params['state_freqs'][target_codon]
        else:
            prob *= self.params["nuc_freqs"][ MOLECULES.nucleotides.index(target_nuc) ]        
        return prob
    


    def _calc_instantaneous_prob(self, source, target, parameters = None):
        ''' 
            Returns the substitution probability from source to target for mechanistic codon models.
            Arguments "source" and "target" are indices for the relevant codons (0-60).
        
            Third argument can be specified as non-self when we are computing neutral scaling factor.
        ''' 
        if parameters is None:
            parameters = self.params
            
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) != 2:
            return 0.
        else:
            nuc_pair = "".join(sorted(nuc_diff[0] + nuc_diff[1]))
            if self._is_syn(source, target):
                return self._calc_prob(target, nuc_diff[1], nuc_pair, parameters['alpha'])
            else:
                return self._calc_prob(target, nuc_diff[1], nuc_pair, parameters['beta'])


    def _build_scaling_params(self):
        '''
            Build scaling parameters for a dN/dS model.
        '''
        new_parameters = {"alpha": 1., "state_freqs": self.params["state_freqs"], "mu": self.params["mu"]} 
        if self.scale_matrix == "average":
            try:
                new_parameters["beta"] = self.params["hetmodel_mean_dnds"]
            except:
                raise KeyError("\n\nCannot compute the matrix scaling factor for heterogeneous codon models without a specified average dN/dS value.")
                
        elif self.scale_matrix == "neutral":
            new_parameters["beta"] = 1. 
                   
        else:
            raise AssertionError("\n\nCan only construct a dN/dS model scaling factor under 'average' and 'neutral' scaling.")
        
        return new_parameters






class MutSel_Matrix(MatrixBuilder):    
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing mutation-selection balance model instantaneous matrices, according to the HalpernBruno 1998 model.
        Here, this model is extended such that it can be used for either nucleotide or codon. This class will automatically detect which one you want based on the provided state frequencies or fitness values.

    '''

    def __init__(self, *args):
        super(MutSel_Matrix, self).__init__(*args)
        try:
            self._size = len(self.params["state_freqs"])
        except:
            self._size = len(self.params["fitness"])
        self.scale_matrix = "neutral"
        if self._size == 4:
            self._code = MOLECULES.nucleotides
        elif self._size == 61:
            self._code = MOLECULES.codons
        else:
            raise ValueError("\n\nMutSel model matrices must be of dimensions 4x4 or 61x61.")



        
    def _calc_fixrate_state_freqs(self, source, target, nuc_diff, parameters):
        ''' 
            Calculate fixation probability using state frequencies and mutation rates.
        '''
            
        pi_i  = parameters['state_freqs'][source]           # source frequency
        pi_j  = parameters['state_freqs'][target]           # target frequency 
        mu_ij = parameters["mu"][nuc_diff]                  # source -> target mutation rate
        mu_ji = parameters["mu"][nuc_diff[1] + nuc_diff[0]] # target -> source mutation rate
        
        # If either frequency is equal to 0, then the rate is 0.
        if abs(pi_i) <= ZERO or abs(pi_j) <= ZERO:
            fixation_rate = 0.
        
        # Otherwise, compute scaled selection coefficient as np.log( pi_mu ) = np.log( (mu_ji*pi_j)/(mu_ij*pi_i) )
        else:
            pi_mu = (mu_ji*pi_j)/(mu_ij*pi_i)
            # If pi_mu == 1, L'Hopitals gives fixation rate of 1 (substitution probability is the forward mutation rate) 
            if abs(1. - pi_mu) <= ZERO:
                fixation_rate = 1. 
            else:
                fixation_rate =  np.log(pi_mu)/(1. - 1./pi_mu)
        return fixation_rate




    def _calc_fixrate_fitness(self, source, target, parameters):
        ''' 
            Calculate fixation probability using fitness values and mutation rates.
        '''            
        sij = parameters['fitness'][target] - parameters['fitness'][source]  
        if abs(sij) <= ZERO:
            fixation_rate = 1. 
        else:
            fixation_rate = (sij)/(1 - np.exp(-1.*sij))
        return fixation_rate           
        
        
                   

    def _calc_instantaneous_prob(self, source, target, parameters = None):
        ''' 
            Calculate the substitution probability from source to target for mutation-selection-balance models.
            Arguments "source" and "target" are indices for the relevant codons (0-60) or nucleotide (0-3).
        
            Third argument can be specified as non-self when we are computing neutral scaling factor.

        '''        
        if parameters is None:
            parameters = self.params
        
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) != 2:
            return 0.
        else:
            if self.params["calc_by_freqs"]:
                fixation_rate = self._calc_fixrate_state_freqs(source, target, nuc_diff, parameters)

            else:
                fixation_rate = self._calc_fixrate_fitness(source, target, parameters)
          
            return fixation_rate * parameters['mu'][nuc_diff]
            
 
 
    def _build_scaling_params(self):
        '''
            Build scaling parameters for a mutation-selection model. Note that only neutral scaling is allowed.
        '''
        assert(self.scale_matrix == "neutral"), "\n\nCan only construct a mutation-selection model scaling factor under 'neutral' scaling."

        new_parameters = {"mu": self.params["mu"]}         
        new_parameters["state_freqs"] = np.repeat(1./self._size, self._size)
        new_parameters["fitness"] = np.ones(self._size)
            
        return new_parameters





class ECM_Matrix(MatrixBuilder):
    ''' 
        Child class of MatrixBuilder. This class implements functions relevant to constructing a matrix specifically for the ECM (described in Kosiol2007) model.
        We support both restricted (instantaneous single changes only) and unrestricted (instantaneous single, double, or triple) versions of this model (see paper for details).
        
        !!! NOTE: The ECM model supports omega (dN/dS) and kappa (TI/TV) ratios in their calculations, and therefore I have included these parameters here. HOWEVER, I do NOT recommend their use.
    
    ''' 
    
    def __init__(self, *args):      
        super(ECM_Matrix, self).__init__(*args)
        if self.model_type == 'ecmrest':
            self.restricted = True
        elif self.model_type == 'ecmunrest':

            self.restricted = False
        else:
            raise ValueError("\n\nECM model must be specified as REST or UNREST, for restricted or unrestricted, respectively.")

        self._size = 61
        self._code = MOLECULES.codons
        self.scale_matrix = "persite" # It's completely unclear how these models should work, so stick with this.
        self._init_empirical_matrix(self.model_type)
        
 


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




    def _calc_instantaneous_prob(self, source, target, parameters = None):
        ''' 
            Returns the substitution probability from source to target for ECM models.
            Arguments "source" and "target" are indices for the relevant codons (0-60).
        '''  
        if parameters is None:
            parameters = self.params
        nuc_diff = self._get_nucleotide_diff(source, target)
        if len(nuc_diff) == 0  or (self.restricted and len(nuc_diff) != 2):
            return 0.
        else:
            kappa_param = self._set_kappa_param(nuc_diff)
            if self._is_syn(source, target):
                return self.emp_matrix[source][target] * parameters['state_freqs'][target] * parameters['alpha'] * kappa_param
            else:
                return self.emp_matrix[source][target] * parameters['state_freqs'][target] * parameters['beta'] * kappa_param






           
