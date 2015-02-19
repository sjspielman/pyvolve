#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################


'''
    This module defines evolutionary model objects, EvoModel() and its child classes Model() and CodonModel().
    All evolutionary models contain information about the substitution process (rate matrix) and information about rate heterogeneity.
    The Model() class uses a single rate matrix, and heterogeneity is modeled using discrete scaling factors and associated probabilities.
    The CodonModel() class is used specifically in cases of codon model (dN/dS or omega) rate heterogeneity. Rate heterogeneity is implemented using a set of matrices with distinct dN/dS values, and each matrix has an associated probability. 
'''

import numpy as np
from copy import deepcopy
from matrix_builder import *


class EvoModels(object):
    ''' 
        Parent class for child classes Model() and CodonModel(). 
    '''
    
    def __init__(self, model_type, **kwargs):
        '''
            The EvoModel class will construct an evolutionary model object which will be used to evolve sequence data.            
            Instantiation requires a single positional argument:
                
                1. **model_type** is type of model (matrix) that is being used. These matrices are described explicitly in the matrix_builder module. Options include the following:
                                
                                    
                    +------------+-----------------------------------------------------------+
                    | model_type |                         Notes                             | 
                    +============+===========================================================+
                    | nucleotide | Arbitrary GTR                                             | 
                    +------------+-----------------------------------------------------------+
                    | JTT        | Jones, Taylor, and Thornton 1994 (amino acids)            |
                    +------------+-----------------------------------------------------------+
                    | WAG        | Whelan and Goldman 2002      (amino acids)                |
                    +------------+-----------------------------------------------------------+
                    | LG         | Le and Gascuel 2008        (amino acids)                  |
                    +------------+-----------------------------------------------------------+
                    | amino_acid | Defaults to LG                                            |
                    +------------+-----------------------------------------------------------+
                    | GY94       | Goldman and Yang 1994, Nielsen and Yang 1998              | 
                    +------------+-----------------------------------------------------------+
                    | MG94       | Muse and Gaut 1994                                        | 
                    +------------+-----------------------------------------------------------+
                    | codon      | Defaults to GY94                                          | 
                    +------------+-----------------------------------------------------------+
                    | ECM        | Kosiol et al. 2007                                        |   
                    +------------+-----------------------------------------------------------+
                    | mutsel     | Halpern and Bruno 2008 (may also be used for nucleotides) |  
                    +------------+-----------------------------------------------------------+
                    
            
            Optional keyword arguments include, 
                1. **params** is a dictionary of parameters pertaining to substitution process. For all models, this includes a vector of stationary frequencies. Each individual evolutionary model will have its own additional parameters. Note that if this argument is not provided, default parameters for your selected model will be assigned.
                
                2. **scale_matrix** = <'yang', 'neutral', 'False/None'>. This argument determines how rate matrices should be scaled. By default, all matrices are scaled according to Ziheng Yang's approach, in which the mean substitution rate is equal to 1. However, for codon models (GY94, MG94), this scaling approach effectively causes sites under purifying selection to evolve at the same rate as sites under positive selection, which may not be desired. Thus, the 'neutral' scaling option will allow for codon matrices to be scaled such that the mean rate of *neutral* subsitution is 1. You may also opt out of scaling by providing either False or None to this argument, although this is not recommended.
                
       '''
    
        
        self.model_type   = model_type.upper()
        self.params       = kwargs.get('params', {})
        self.scale_matrix = kwargs.get('scale_matrix', 'yang') # 'Yang', 'neutral', or False/None

        accepted_models = ['NUCLEOTIDE', 'AMINO_ACID', 'JTT', 'WAG', 'LG', 'CODON', 'GY94', 'MG94', 'MUTSEL']
        assert( self.model_type in accepted_models), "Inappropriate model type specified."
        assert( type(self.params) is dict ), "params argument must be a dictionary."
        
        # Default codon, aa models
        if self.model_type == 'CODON':
            self.model_type = 'GY94'
        if self.model_type == 'AMINO_ACID':
            self.model_type = 'LG'
        self.name = None
          

    def construct_model(self):
        '''
            Construct EvoModel. Setup substitution matrix(ces) and rate heterogeneity probabilities.
            Calls _assign_matrix and _assign_rate_probs, as needed.
            
            Parent class method. Not executed.
        '''
        print "Parent class method. Not executed."

  
    def assign_name(self, name):
        '''
            Assign name to an EvoModel instance. 
            In cases of branch/temporal homogeneity, names are unneeded.
            However, in cases of **branch heterogeneity, each model must be named**. Names are used to map to model flags given in the phylogeny.

        '''
        self.name = name
 
        
    def _assign_matrix(self):
        '''
            Compute and assign substitution matrix(ces).
            
            Parent class method. Not executed.
        '''
        print "Parent class method. Not executed."

        
         
    def _assign_rate_probs(self, category_variable):
        '''
            Compute probabilities for rate heterogeneity.
            By default, equal probabilities will be assigned to all rate categories (either scaling factors for Model() or dN/dS values and corresponding matrices for CodonModel()).
            
            Alternatively, a single optional argument giving user-specified rate categories may be provided.
            If Model(), the argument should be a list of rate factors
            If CodonModel(), the argument should be a list of matrices.
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
            Return the number of rate classes associated with a given model.
        '''
        return len(self.rate_probs)   
            
            
    def codon_model(self):
        '''
            Return True if the model is a CodonModel(), and return False otherwise.
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
        '''
            Model() construction requires arguments as described under the EvoModel() documentation. 
        '''

        super(Model, self).__init__(*args, **kwargs)
        self.rate_factors = [1.]  # Default Rate heterogeneity factors (default is site homogeneity).
        
        

    def construct_model(self, **kwargs):
        '''
            Construct Model by building the substitution matrix and defining rate heterogeneity probabilities.
            
            Optional keyword arguments include, 
            
                1. **rate_factors**, a list/numpy array of scalar factors for rate heterogeneity. Default: rate homogeneity.
                2. **rate_probs**, a list/numpy array of probabilities (which sum to 1!) for each rate category. Default: equal.
                3. **alpha**, the alpha shape parameter which should be used to draw rates from a discrete gamma distribution. Supply this argument to have gamma-distribtued rates.
                4. **num_categories**, the number of rate categories to create. Supply this argument to draw a certain number of rates from a gamma distribution.               
                
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
            Draw **k** rates from a discrete gamma distribution with shape parameter **alpha**.
        '''    
        self.rate_factors = np.random.gamma(alpha, scale = alpha, size = k) 
        

    def _assign_matrix(self):
        '''
            Construct the model rate matrix, Q, based on model_type by calling the matrix_builder module. 
        '''
        if self.model_type == 'NUCLEOTIDE':
            self.matrix = nucleotide_Matrix(self.params, self.scale_matrix)()
        
        elif self.model_type in ['JTT', 'WAG', 'LG']:
            self.params["aa_model"] = self.model_type
            self.matrix = aminoAcid_Matrix(self.params, self.scale_matrix)()
        
        elif self.model_type == 'GY94' or self.model_type == 'MG94':
            self.matrix = mechCodon_Matrix(self.params, self.model_type, self.scale_matrix)()
        
        elif self.model_type == 'ECM':
            self.matrix = ECM_Matrix(self.params, self.scale_matrix)()
        
        elif self.model_type == 'MUTSEL':
            self.matrix = mutSel_Matrix(self.params, self.scale_matrix)()
        
        else:
            raise AssertionError("WHAT ARE WE DOING HERE?! Please contact stephanie.spielman@gmail.com .")
            
    
 
    def _sanity_rate_factors(self):
        '''
            Perform sanity checks on rate heterogeneity set-up:
                1. Ensure that rate_factors is type np.array
                2. Ensure that rates are properly normalized with probabilities
        '''
        
        if type( self.rate_factors ) is list:
            self.rate_factors = np.array( self.rate_factors )
        if abs( 1. - np.sum(self.rate_probs * self.rate_factors)) > ZERO:
            self.rate_factors /= np.sum(self.rate_factors * self.rate_probs)









class CodonModel(EvoModels):

    '''
        Defines a CodonModel() object. This class is reserved specifically for cases of **codon model heterogeneity**, in which dN/dS (omega) varies across sites, and hence matrices must vary.
    '''               
                
        
    def __init__(self, *args, **kwargs):
        
        '''
            CodonModel() instantiation requires arguments as described under the EvoModel() documentation. Importantly, the first positional argument, the **params** dictionary, must contain state_freqs, mutational parameters, a list of betas (dN), and an associated list of alphas (dS). For example, model with 3 categories of dN/dS heterogeneity might look like, ```{'state_freqs':f, 'kappa':2.75, 'beta':[1, 2.5, 0.5], 'alpha':[1, 1, 1]}```
        '''
        
        super(CodonModel, self).__init__(*args, **kwargs)
        assert( self.model_type == 'GY94' or self.model_type == 'MG94' or self.model_type == 'ECM' ), "CodonModels supported include only GY94, MG94, and ECM."

   
    def construct_model(self, **kwargs):
        '''
            Construct CodonModel by building substitution matrices and defining rate heterogeneity probabilities.
            
            Optional keyword arguments include, 
            
                1. **rate_probs**, a list/numpy array of probabilities (which sum to 1!) for each dN/dS category. Default: equal.

        '''
        self._assign_matrix()
        self.rate_probs = kwargs.get('rate_probs', None)
        self._assign_rate_probs( self.matrices )            
    
        
    
    
    def _assign_matrix(self):
        '''
            Construct each model rate matrix, Q, to create a list of codon-model matrices.
        '''
        self.matrices = []
        assert( len(self.params['beta']) == len(self.params['alpha']) ), "num dn is not same as num ds"
        for i in range(len( self.params['beta'] )):
            temp_params = deepcopy(self.params)
            temp_params['beta'] = self.params['beta'][i]
            temp_params['alpha'] = self.params['alpha'][i]
            if self.model_type == 'GY94' or self.model_type == 'MG94':
                self.matrices.append( mechCodon_Matrix(temp_params, self.model_type, self.scale_matrix)() )
            else:
                self.matrices.append( ECM_Matrix(temp_params, self.scale_matrix)() )     
        assert( len(self.matrices) > 0), "You have no matrices for your CodonModel :("


 
 
 
 