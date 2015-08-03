#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################


'''
    Define evolutionary model objects.
'''

import numpy as np
from copy import deepcopy
from matrix_builder import *
from genetics import *
MOLECULES = Genetics()

class Model():
    ''' 
        This class defines evolutionary model objects.
        All evolutionary models contain information about the substitution process (rate matrix) and information about rate heterogeneity.
        Note that, in cases of rate heterogeneity, non-dN/dS models use a single rate matrix and model heterogeneity using discrete scaling factors and associated probabilities.
        Alternatively, rate heterogeneity in dN/dS models is implemented using a set of matrices with distinct dN/dS values, and each matrix has an associated probability. 
    '''
    
    def __init__(self, model_type, parameters=None, **kwargs):
        '''
            The Model class will construct an evolutionary model object which will be used to evolve sequence data.            
            Instantiation requires a single positional argument (but a second one is recommended, read on!):
                
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
                    | MTMAM      | Yang, Nielsen, and Hasagawa 1998   (amino acids)          |
                    +------------+-----------------------------------------------------------+
                    | MTREV24    | Adachi and Hasegawa  1996 (amino acids)                   |
                    +------------+-----------------------------------------------------------+
                    | DAYHOFF    | Dayhoff, Schwartz, and Orcutt  1978 (amino acids)         |
                    +------------+-----------------------------------------------------------+
                    | AB         | Mirsky, Kazandjian, and Anisimova 2015    (amino acids)   |
                    +------------+-----------------------------------------------------------+
                    | GY         | Goldman and Yang 1994 (modified), Nielsen and Yang 1998   | 
                    +------------+-----------------------------------------------------------+
                    | MG         | Muse and Gaut 1994                                        | 
                    +------------+-----------------------------------------------------------+
                    | codon      | Defaults to GY-style model                                | 
                    +------------+-----------------------------------------------------------+
                    | ECM        | Kosiol et al. 2007                                        |   
                    +------------+-----------------------------------------------------------+
                    | mutsel     | Halpern and Bruno 2008 (may also be used for nucleotides) |  
                    +------------+-----------------------------------------------------------+

            To use your own rate matrix (which you must create on your own), enter "custom" for the model_type argument, and provide the custom matrix (numpy array or list of lists) in the **parameters** dictionary with the key "matrix". Please note that pyvolve stores nucleotides, amino acids, and codons in alphabetical order of their abbreviations:
            *  Nucleotides: A, C, G, T
            *  Amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
            *  Codons:      AAA, AAC, AAG, AAT, ACA, ... TTG, TTT [note that stop codons should be *excluded*]
            
            If you wish to evolve *custom states* (neither nucleotide, amino acids, nor codons), for instance to evolve characters, also include the key "code" in the parameters dictionary. The associated value should be a list of strings, e.g. ["0", "1", "2"], and the length of this list should be the same as a dimension of the square custom matrix provided. Note that this argument is not required if wish to evolve nucleotides, amino-acids, and/or codons. 

            Please be careful here - Pyvolve takes your matrix (mostly) at face-value (provided it has proper dimensions and rows sum to 0). In particular, the matrix will not be scaled!!! 


            A second positional argument, **parameters** may additionally be specified. This argument should be a dictionary of parameters pertaining to substitution process. Each individual evolutionary model will have its own parameters. Note that if this argument is not provided, default parameters for your selected model will be assigned. Note that this argument is **required** for mechanistic codon (dN/dS) models, as this rate ratio must be assigned!
            
                            
            Optional keyword arguments include, 
                1. **scale_matrix** = <'yang', 'neutral'>. This argument determines how rate matrices should be scaled. By default, all matrices are scaled according to Ziheng Yang's approach, in which the mean substitution rate is equal to 1. However, for codon models (GY-style and MG-style), this scaling approach effectively causes sites under purifying selection to evolve at the same rate as sites under positive selection, which may not be desired. Thus, the 'neutral' scaling option will allow for codon matrices to be scaled such that the mean rate of *neutral* subsitution is 1.            
                2. **name**, the name for a Model object. Names are not needed in cases of branch homogeneity, but when there is **branch heterogeneity**, names are required to map the model to the model flags provided in the phylogeny.
                3. **rate_factors**, for specifying rate heterogeneity in nucleotide or amino acid models. This argument should be a list/numpy array of scalar factors for rate heterogeneity. Default: rate homogeneity.
                4. **rate_probs**, for specifying rate heterogeneity probabilities in nucleotide, amino acid, or codon models. This argument should be a list/numpy array of probabilities (which sum to 1!) for each rate category. Default: equal.
                5. **alpha**, for specifying rate heterogeneity in nucleotide or amino acid models if gamma-distributed heterogeneity is desired. The alpha shape parameter which should be used to draw rates from a discrete gamma distribution.
                6. **num_categories**, for specifying the number of gamma categories to draw for rate heterogeneity in nucleotide or amino acid models. Should be used in conjunction with the "alpha" parameter. Default: 4.               
                7. **pinv**, for specifying a proportion of invariant sites when gamma heterogeneity is used. When specifying custom rate heterogeneity, a proportion of invariant sites can be specified simply with a rate factor of 0.
                8. **save_custom_frequencies**, for specifying a file name in which to save the state frequencies from a custom matrix. Pyvolve automatically computes the proper frequencies and will save them to a file named "custom_matrix_frequencies.txt", and you can use this argument to change the file name. Note that this argument is really only relevant for custom models.

       '''
    
        
        
        self.model_type   = model_type.upper()
        if parameters == None:
            self.params = {}
        else:
            self.params = parameters
        
        self.scale_matrix = kwargs.get('scale_matrix', 'yang')     # 'yang' or 'neutral'
        self.name         = kwargs.get('name', None)               # Can be overwritten through .assign_name()
        self.rate_probs   = kwargs.get('rate_probs', None)         # Default is no rate hetereogeneity
        self.rate_factors = kwargs.get('rate_factors', np.ones(1)) # Default is no rate hetereogeneity
        self.alpha        = kwargs.get('alpha', None)
        self.k_gamma      = kwargs.get('num_categories', None)
        self.pinv         = kwargs.get('pinv', 0.)                 # If > 0, this will be the last entry in self.rate_probs, and 0 will be the last entry in self.rate_factors
        self._save_custom_matrix_freqs = kwargs.get('save_custom_frequencies', "custom_matrix_frequencies.txt")
        self.code         = None
        
        self._sanity_model()
        self._check_codon_model()
        self._construct_model()    

        
        
 
    def _assign_code(self):
        ''' 
            Assign genetic code or custom code provided specifically for a custom matrix.
        '''    
        if "code" in self.params:
            self.code = self.params["code"]
        else:
            dim = len(self.params['state_freqs']) 
            if dim == 4:
                self.code = MOLECULES.nucleotides
            elif dim == 20:
                self.code = MOLECULES.amino_acids
            elif dim == 61:
                self.code = MOLECULES.codons
            else:
                raise AssertionError("\n\nUnknown genetic code.")



    def _sanity_model(self):
        '''
            A series of brief sanity checks on Model init.
        '''
        
        accepted_models = ['NUCLEOTIDE', 'JTT', 'WAG', 'LG', 'AB', 'MTMAM', 'MTREV24', 'DAYHOFF', 'CODON', 'GY', 'MG', 'MUTSEL', 'ECM', 'ECMREST', 'ECMUNREST', 'CUSTOM']
        assert( self.model_type in accepted_models), "\n\nInappropriate model type specified."
        assert( type(self.params) is dict ), "\n\nThe parameters argument must be a dictionary."
        
        # Default codon, ecm models
        if self.model_type == 'CODON':
            self.model_type = 'GY'
            print("Using default codon model, GY-style.")
        if self.model_type == 'ECM':
            self.model_type = 'ECMREST'
            print("Using restricted ECM model.")
        if self.model_type == 'CUSTOM':
            assert("matrix" in self.params), "\n\nTo use a custom model, you must provide a matrix in your params dictionary under the key 'matrix'. Your matrix must be symmetric and rows must sum to 1 (note that pyvolve will normalize the matrix as needed). Also note that pyvolve orders nucleotides, amino acids, and codons alphabetically by their abbreviations (e.g. amino acids are ordered A, C, D, ... Y)."          
            if "state_freqs" in self.params:
                print("Since you have specified a custom matrix, your provided state frequencies will be *ignored*. Pyvolve will calculate them for you, from the provided matrix. These frequencies will be saved, for your convenience, to a file",self._save_custom_matrix_freqs,".")
        


    def _check_codon_model(self):
        '''
            Determine if this is a heterogenous codon model and assign self.codon_model accordingly.
            We also remove any omega keys are replace right away with beta.
        '''
        self.codon_model = False
        
        if "omega" in self.params:
            self.params["beta"] = self.params["omega"]
            self.params.pop("omega")
        
        if "beta" in self.params:
            # Iterable?
            try:
                (x for x in self.params["beta"])
                self.codon_model = True
            except:
                self.codon_model = False
                try:
                    self.params["beta"] = float(self.params["beta"])
                except:
                    raise AssertionError("\n\nTo specify a dN/dS value, provide either an integer or float (for rate homogeneity) or a list/numpy array of values (for rate heterogeneity) using the key 'omega' (or keys 'alpha' and 'beta' for dS and dN, respectively, separately).")


    def _construct_model(self):
        '''
            Construct the evolutionary model.
        '''
        
        # Assign matrix (or matrices, for a heterogeneous codon model)
        self._assign_matrix()
       
        # Assign rate heterogeneity
        if self.codon_model:
            self._assign_rate_probs(self.matrix)
        else:
            self._assign_rates()
        
        # Assign code
        self._assign_code()



    def _assign_matrix(self):
        '''
            Construct the model rate matrix, Q, based on model_type by calling the matrix_builder module. Alternatively, call the method self._assign_codon_model_matrices() if we have a heterogenous codon model.
        '''
        
        if self.model_type == 'NUCLEOTIDE':
            self.matrix = nucleotide_Matrix(self.params, self.scale_matrix)()
        
        elif self.model_type in ['JTT', 'WAG', 'LG', 'AB', 'MTMAM', 'MTREV24', 'DAYHOFF']:
            self.params["aa_model"] = self.model_type
            self.matrix = aminoAcid_Matrix(self.params, self.scale_matrix)()
        
        elif self.model_type == 'GY' or self.model_type == 'MG':
            if self.is_codon_model():
                self._assign_codon_model_matrices()
            else:
                self.matrix = mechCodon_Matrix(self.params, self.model_type, self.scale_matrix)()
        
        elif 'ECM' in self.model_type:
            self.params["rest_type"] = self.model_type.split("ECM")[1]
            self.matrix = ECM_Matrix(self.params, self.scale_matrix)()
        
        elif self.model_type == 'MUTSEL':
            temp_mat = mutSel_Matrix(self.params, self.scale_matrix)
            self.matrix = temp_mat()
            if "state_freqs" not in self.params:
                self.params["state_freqs"] = temp_mat.extract_state_freqs(self.matrix)
        
        elif self.model_type == 'CUSTOM':
            self._assign_custom_matrix()
        
        else:
            raise AssertionError("You have reached this in error! Please file a bug report, with this error, at https://github.com/sjspielman/pyvolve/issues .")



    def _assign_custom_matrix(self):
        '''
            Create rate matrix from user-provided. Also, extract the equilibrium frequencies and save to file.
            We must check that dimensions are square and that rows sum to 1. If they sum to 1 with a tolerance of 1e-4, we accept it and "re-tolerize". If not, return an error.
            Further, if a custom code was specified, check it!
        '''
        
        custom_matrix = np.array( self.params['matrix'] )
        
        # Check shape and code. Assigns code attribute, as well.
        if "code" in self.params:
            assert(type(self.params["code"]) is list), "\n When providing a custom code for your custom matrix, provide a list of *strings*. Each item in this list is a state (so states can be arbitrarily named!), and therefore the length of this list should equal a dimension of your square matrix!"
            for item in self.params["code"]:
                assert(type(item) is str), "\n When providing a custom code for your custom matrix, provide a list of *strings*. Each item in this list is a state (so states can be arbitrarily named!), and therefore the length of this list should equal a dimension of your square matrix!"
            dim = len(self.params["code"])
            assert( custom_matrix.shape == (dim, dim) ), "\n The dimensions for your custom matrix must be the same as your custom code!" 
        else:
            assert( custom_matrix.shape == (4,4) or custom_matrix.shape == (20,20) or custom_matrix.shape == (61,61) ), "\n Custom transition matrix must be symmetric with dimensions 4x4 (nucleotides), 20x20 (amino-acids), or codons (61x61). If you wish to use a custom code which does not have these states, then specify this code with the argument custom_code."
            dim = custom_matrix.shape[0]
         
        # Check that sums to zero with a relatively permissive tolerance
        assert ( np.allclose( np.zeros(dim), np.sum(custom_matrix, 1) , rtol=1e-4) ), "Rows in custom transition matrix do not sum to 0."
        
        # "Re-normalize" matrix with better tolerance and confirm
        for s in range(dim):
            temp_sum = np.sum(custom_matrix[s]) - np.sum(custom_matrix[s][s])
            custom_matrix[s][s] = -1. * temp_sum
            assert ( abs(np.sum(custom_matrix[s])) < ZERO ), "Re-normalized row in custom transition matrix does not sum to 0."

        # Now, calculate state frequencies from this matrix.. hacky. leave me alone.
        temp = mutSel_Matrix( {"state_freqs": np.repeat(0.25, 4) })
        state_freqs = temp.extract_state_freqs(custom_matrix, size = dim)
        del temp

        self.matrix = custom_matrix
        self.params["state_freqs"] = state_freqs
        np.savetxt(self._save_custom_matrix_freqs, state_freqs) 
      
      
      
      
    def _assign_codon_model_matrices(self):
        '''
            Construct each model rate matrix, Q, to create a list of codon-model matrices. Also, perform some sanity checks.
        '''
        
        # Sanity checks. Note that any 'omega' keys have already been replaced by 'beta', in the self._check_codon_model method here.
        assert("beta" in self.params), "You must provide dN values (using either the key 'beta' or 'omega') in params dictionary to run this model!"
        
        # alpha should be the same length as beta
        if "alpha" in self.params:
            assert( len(self.params['beta']) == len(self.params['alpha']) ), "To specify both dN and dS heterogeneity, provide lists (or numpy arrays), of the *same lengths*, for keys 'alpha' and 'beta'."
        else:
            self.params['alpha'] = np.repeat(1., len(self.params["beta"]))  
        
        # We need to add state_freqs if missing, as well, since the actual params dictionary does not get passed to matrixBuilder
        if "state_freqs" not in self.params:
            self.params["state_freqs"] = np.repeat(1./61, 61)
            

        # Construct matrices
        self.matrix = []
        for i in range(len( self.params['beta'] )):
            temp_params = deepcopy(self.params)
            temp_params['beta'] = self.params['beta'][i]
            temp_params['alpha'] = self.params['alpha'][i]
            self.matrix.append( mechCodon_Matrix(temp_params, self.model_type, self.scale_matrix)() )
        assert(len(self.matrix) > 0), "Matrices for a heterogeneous codon model were improperly constructed."





    def _assign_rates(self):
        '''
            Assign and sanity-check site heterogeneity rates for non-codon models. Note that heterogeneity is *ignored* for ECM models, because it's entirely unclear how/why this should be implemented.
        '''
        if "ECM" not in self.model_type:
            # Draw gamma rates if specified
            if self.alpha is not None:
                assert(self.pinv >= 0. and self.pinv <= 1.), "\n\nThe proportion of invariant sites must be a value between 0 and 1 (inclusive)."
                self._draw_gamma_rates()
            self._assign_rate_probs(self.rate_factors)
            self._sanity_rate_factors()
        else:
            self.rate_probs = np.ones(1)



 
    def _draw_gamma_rates(self):
        '''
            Function to draw and assign rates from a gamma distribution, if specified. By default, 4 categories are drawn.
            If a proportion of invariant sites has been specified, draw gamma rates for remaining probability.
        '''       
        if self.k_gamma is None and self.rate_probs is None:
            self.k_gamma = 4
        elif self.k_gamma is None and self.rate_probs is not None:
            self.k_gamma = len(self.rate_probs)                
        elif self.k_gamma is not None and self.rate_probs is not None:
            assert(self.k_gamma == len(self.rate_probs)), "\n\nWhen specifying custom probabilities for a gamma distribution, the length of your rate_probs list must equal the num_categories."
        self.rate_factors = list(np.random.gamma(self.alpha, scale = self.alpha, size = self.k_gamma))
        if self.pinv > 0.:
            self.rate_factors.append(0.)
        self.rate_factors = np.array(self.rate_factors)



    def _assign_gamma_pinv_rate_probs(self):
        '''
            Function to compute rate probabilities under the specific Gamma + Pinv rate heterogeneity scheme.
        '''
        # Default
        if self.rate_probs is None:
            remaining_prob = 1. - self.pinv
            self.rate_probs = list(np.repeat( remaining_prob/self.k_gamma, self.k_gamma ))        
        # Custom
        else:
            self.rate_probs = list(self.rate_probs)
        self.rate_probs.append(self.pinv)
        self.rate_probs = np.array( self.rate_probs ) 
            
       
     
         
    def _assign_rate_probs(self, category_variable):
        '''
            Compute probabilities for rate heterogeneity.
            By default, equal probabilities will be assigned to all rate categories (either scaling factors or dN/dS values and corresponding matrices for a heterogeneous codon model).
            
            If regular Model, the category_variable argument should be a list of rate factors
            If codon het Model, the category_variable argument should be a list of matrices.
        '''
        
        # Gamma + Pinv
        if self.alpha is not None and self.pinv > 0.:
            self._assign_gamma_pinv_rate_probs()
        
        # Nothing provided, default.
        if self.rate_probs is None:
            self.rate_probs = np.repeat( 1./len(category_variable), len(category_variable) )  

        ### Perform some checks ###
        
        # Ensure sums to 1
        assert(abs(1. - np.sum(self.rate_probs)) <= ZERO), "\n\nProvided rate probabilities (rate_probs list) must sum to 1.\nNote: if you are specifying Gamma+Pinv heteregeneity with custom probabilities, ensure that the sum of pinv and your rate_probs list is equal to 1."
            
        # Ensure numpy array
        try:
            self.rate_probs = np.array( self.rate_probs )
        except:
            raise AssertionError("\n\nRate probabilities improperly specified.")                
            
        # Size sanity check.
        assert( len(self.rate_probs) == len(category_variable) ), "\n\nDifferent numbers of probabilities and matrices."
            



    def _sanity_rate_factors(self):
        '''
            Perform sanity checks on rate heterogeneity set-up:
                1. Ensure that rate_factors is type np.array
                2. Ensure that rates are properly normalized with probabilities
        '''
        # Ensure numpy array
        try:
            self.rate_factors = np.array( self.rate_factors )
        except:
            raise AssertionError("\n Rate factors improperly specified.")                

        if abs( 1. - np.sum(self.rate_probs * self.rate_factors)) > ZERO:
            self.rate_factors /= np.sum(self.rate_factors * self.rate_probs)




    def num_classes(self):
        ''' 
            Return the number of rate classes associated with a given model.
        '''
        return len(self.rate_probs)   

            
    

    def assign_name(self, name):
        '''
            Assign name to a Model instance. 
            In cases of branch/temporal homogeneity, names are unneeded.
            However, in cases of **branch heterogeneity, each model must be named**. Names are used to map to model flags given in the phylogeny.
            NOTE that name can also be assigned as a keyword argument when initializing a Model object.
        '''
        self.name = name  
        
        
        
    def is_codon_model(self):
        '''
            Return True if the model is a heterogeneous codon model and return False otherwise.
        '''
        return self.codon_model
          
        
