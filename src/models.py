#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################


'''
    EvoModel() classes: Model() and CodonModel().
'''

import numpy as np
from copy import deepcopy
from matrix_builder import *


class EvoModels(object):
    ''' 
        Parent class for Model(), CodonModel().
        REQUIRED POSITIONAL ARGUMENTS:
            *params*      = dictionary of parameters pertaining to substitution process. For all models, this includes a vector of stationary frequencies. Each individual evolutionary model will have its own additional parameters.
            *model_type*  = type of model (matrix) we'll use. Options include the following:
                1. nucleotide
                2. amino_acid
                3. codon/GY94/MG94 (if codon specified, result will be GY94 matrix)
                4. ECM             (Kosiol et al. 2007)
                5. mutsel          (mutation-selection model, for nucleotide and codons)  
    '''
    
    def __init__(self, params, model_type, **kwargs):
        self.params       = params
        self.model_type   = model_type
        assert( type(self.params) is dict ), "params argument must be a dictionary."
        assert( self.model_type == 'nucleotide' or self.model_type == 'amino_acid' or self.model_type == 'codon' or self.model_type == 'GY94' or self.model_type == 'MG94' or self.model_type == 'ECM' or self.model_type == 'mutsel' ), "Inappropriate model type specified."
        if self.model_type == 'codon':
            self.model_type = 'GY94'
        self.name = None
          

    def construct_model(self):
        '''
            Setup matrix/ces, rate probabilities.
        '''
        print "Parent class method. Not executed."

  
    def assign_name(self, name):
        '''
            Assign name to a model. 
            Names *must* be used in cases of branch heterogeneity in order to map to model flags in the phylogeny. Otherwise, no need.
        '''
        self.name = name
 
        
    def _assign_matrix(self):
        '''
            Compute and assign Q matrix/ces. 
        '''
        print "Parent class method. Not executed."

        
         
    def _assign_rate_probs(self, category_variable):
        '''
            Determine rate heterogeneity class/category probabilities, either from provided argument or, as default, set all categories to equal probabilities.
            Argument *category_variable* is either self.rate_factors (if Model()) or self.matrices (if CodonModel())
        '''
        
     
        # Nothing provided, default.
        if self.rate_probs is None:
            self.rate_probs = np.repeat( 1./len(category_variable), len(category_variable) )  

        # Assign according to user-provided argument.        
        else:
            # Ensure numpy array
            if type( self.rate_probs ) is list:
                self.rate_probs = np.array( self.rate_probs )
            
            # Size sanity check.
            assert( len(self.rate_probs) == len(category_variable) ), "Different numbers of probabilities and matrices..."
            
            # Check sum to 1. If not, renormalize.
            if abs( 1. - np.sum(self.rate_probs)) > ZERO:
                self.rate_probs /= np.sum(self.rate_probs)


    def num_classes(self):
        ''' 
            How many rate classes? 
        '''
        return len(self.rate_probs)   
            
            
    def codon_model(self):
        '''
            Is this a codon model?"
        '''
        if isinstance(self, CodonModel):
            return True
        else:
            return False
            
        
        

    

      
    
    
    
    
    
class Model(EvoModels):
    '''
        Defines a Model() object. Used for models for which heterogeneity is determined by a scalar factor (all but dN/dS models).
    '''
    
    def __init__(self, *args, **kwargs):
        super(Model, self).__init__(*args, **kwargs)
        self.rate_factors = [1.]  # Default Rate heterogeneity factors (default is site homogeneity).
        
        

    def construct_model(self, **kwargs):
        '''
            Setup matrix/ces, rate probabilities.
            Optional arguments:
                *rate_factors*    = list/np.array of scalars for rate heterogeneity
                *rate_probs*      = list/np.array of probabilities (summing to 1!) for each rate category
                *alpha*           = alpha shape parameter if rates should be drawn from a gamma
                *num_categories*  = number of rates. Use in conjunction with alpha to draw that many factors!
                
        '''
        self.rate_factors = kwargs.get('rate_factors', np.array([1.]))    
        self.rate_probs   = kwargs.get('rate_probs', None )
        alpha = kwargs.get('alpha', None)
        k     = kwargs.get('num_categories', None)

        if alpha is not None:
            if k is None:
                if self.rate_probs is not None:
                    k = len(self.rate_probs)
                else:
                    raise AssertionError("You must specify the number of categories (argument num_categories=...) when constructing model if you want gamma rates.") 
            self._assign_gamma_rates(alpha, k)
        
        self._assign_matrix()
        self._assign_rate_probs(self.rate_factors)
        self._sanity_rate_factors()
        
        
        
        
    def _assign_gamma_rates(self, alpha, k):
        '''
            Draw *k* rates from a gamma distribution with shape parameter *alpha*.
        '''    
        self.rate_factors = np.random.gamma(alpha, scale = alpha, size = k) 
        

    def _assign_matrix(self):
        '''
            Construct Q model matrix. 
        '''
        if self.model_type == 'nucleotide':
            self.matrix = nucleotide_Matrix(self.params)()
        
        elif self.model_type == 'amino_acid':
            self.matrix = aminoAcid_Matrix(self.params)()
        
        elif self.model_type == 'GY94' or self.model_type == 'MG94':
            self.matrix = mechCodon_Matrix(self.params, self.model_type)()
        
        elif self.model_type == 'ECM':
            self.matrix = ECM_Matrix(self.params)()
        
        elif self.model_type == 'mutsel':
            self.matrix = mutSel_Matrix(self.params)()
        
        else:
            raise AssertionError("WHAT ARE WE DOING HERE?! Please contact stephanie.spielman@gmail.com .")
            
    
 
    def _sanity_rate_factors(self):
        '''
            Sanity checks on site heterogeneity set-up.
            1. Ensure that rate_factors is type np.array
            2. Ensure that rates are properly normalized with probabilities
        '''
        
        if type( self.rate_factors ) is list:
            self.rate_factors = np.array( self.rate_factors )
        if abs( 1. - np.sum(self.rate_probs * self.rate_factors)) > ZERO:
            self.rate_factors /= np.sum(self.rate_factors * self.rate_probs)









class CodonModel(EvoModels):
    '''
        Defines a CodonModel() object. 
        Class reserved for cases of *codon model heterogeneity* where dN/dS (or other) varies, and hence matrices vary.
    '''               
                
        
    def __init__(self, *args, **kwargs):
        ''' 
            For CodonModel, argument *params* should be a dictionary containing state_freqs, mutational parameters, a list of betas, a list of alphas.
                Example dictionary, {'state_freqs':f, 'kappa':2.75, 'beta':[1, 2.5, 0.5], 'alpha':[1, 1, 1]}
        '''
        super(CodonModel, self).__init__(*args, **kwargs)
        assert( self.model_type == 'GY94' or self.model_type == 'MG94' or self.model_type == 'ECM' ), "CodonModels supported include only GY94, MG94, and ECM."

   
    def construct_model(self, **kwargs):
        '''
            Setup matrices, rate probabilities.
            Optional arguments:
                *rate_probs*      = list/np.array of probabilities (summing to 1!) for each rate category
                *alpha*           = alpha shape parameter if rates should be drawn from a gamma
                *num_categories*  = number of rates. Use in conjunction with alpha to draw that many factors!

        '''
        self._assign_matrix()
        self.rate_probs = kwargs.get('rate_probs', None)
        self._assign_rate_probs( self.matrices )            
    
    
    def _assign_matrix(self):
        '''
            Construct list of Q matrices.
        '''
        self.matrices = []
        assert( len(self.params['beta']) == len(self.params['alpha']) ), "num dn is not same as num ds"
        for i in range(len( self.params['beta'] )):
            temp_params = deepcopy(self.params)
            temp_params['beta'] = self.params['beta'][i]
            temp_params['alpha'] = self.params['alpha'][i]
            if self.model_type == 'GY94' or self.model_type == 'MG94':
                self.matrices.append( mechCodon_Matrix(temp_params, self.model_type)() )
            else:
                self.matrices.append( ECM_Matrix(temp_params)() )     
        assert( len(self.matrices) > 0), "You have no matrices for your CodonModel :("


 
 
 
 