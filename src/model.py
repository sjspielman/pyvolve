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
from .matrix_builder import *
from .genetics import *
from .parameters_sanity import *
from scipy.stats import gamma
from scipy.special import gammainc
ZERO      = 1e-8
MOLECULES = Genetics()

class Model():
    ''' 
        This class defines evolutionary model objects.
        All evolutionary models contain information about the substitution process (rate matrix) and information about rate heterogeneity.
        Note that, in cases of rate heterogeneity, non-dN/dS models use a single rate matrix and model heterogeneity using discrete scaling factors and associated probabilities.
        Alternatively, rate heterogeneity in dN/dS models is implemented using a set of matrices with distinct dN/dS values, and each matrix has an associated probability. 
    '''
    
    def __init__(self, model_type, parameters = None, **kwargs):
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
            If you supply equilibrium frequencies for your custom model ("state_freqs" key in parameters dictionary), pyvolve will require a corresponding symmetric matrix, and the final substitution matrix will be calculated from your symmetric matrix and the provided frequencies. If frequencies are not provided, then they will be calculated directly from your provided matrix. In this circumstance, your matrix will be taken at face value.
            
            A second positional argument, **parameters** may additionally be specified. This argument should be a dictionary of parameters pertaining to substitution process. Each individual evolutionary model will have its own parameters. Note that if this argument is not provided, default parameters for your selected model will be assigned. Note that this argument is **required** for mechanistic codon (dN/dS) models, as this rate ratio must be assigned!
            
                            
            Optional keyword arguments include, 
                1. **name**, the name for a Model object. Names are not needed in cases of branch homogeneity, but when there is **branch heterogeneity**, names are required to map the model to the model flags provided in the phylogeny.
                2. **rate_factors**, for specifying rate heterogeneity in nucleotide or amino acid models. This argument should be a list/numpy array of scalar factors for rate heterogeneity. Default: rate homogeneity.
                3. **rate_probs**, for specifying rate heterogeneity probabilities in nucleotide, amino acid, or codon models. This argument should be a list/numpy array of probabilities (which sum to 1!) for each rate category. Default: equal.
                4. **alpha**, for specifying rate heterogeneity in nucleotide or amino acid models if gamma-distributed heterogeneity is desired. This value indicates the alpha shape parameter which should be used to draw rates from a discrete gamma distribution.
                5. **num_categories**, for specifying the number of gamma categories to draw for rate heterogeneity in nucleotide or amino acid models. Should be used in conjunction with the "alpha" parameter. Default: 4.               
                6. **pinv**, for specifying a proportion of invariant sites when gamma heterogeneity is used. When specifying custom rate heterogeneity, a proportion of invariant sites can be specified simply with a rate factor of 0.
                7. **save_custom_frequencies**, for specifying a file name in which to save the state frequencies from a custom matrix. When necessary, pyvolve automatically computes the proper frequencies and will save them to a file named "custom_matrix_frequencies.txt", and you can use this argument to change the file name. Note that this argument is really only relevant for custom models.
                8. **neutral_scaling**, for specifying that **codon models** (GY, MG) be scaled such that the mean rate of neutral substitution per unit time is 1. By default, codon models are scaled according to set the mean substitution rate to be 1. Setting this parameter to True will scale codon models neutrally. Note that this argument is only relevant for codon models and is ignored in other cases. Default: False.
       '''
    
        
        
        self.model_type   = model_type.lower()
        if parameters is None:
            self.params = {}
        else:
            self.params = parameters
        self.name                      = kwargs.get('name', None)               # Can be overwritten through .assign_name()
        self.rate_probs                = kwargs.get('rate_probs', None)         # Default is no rate hetereogeneity
        self.rate_factors              = kwargs.get('rate_factors', np.ones(1)) # Default is no rate hetereogeneity
        self.alpha                     = kwargs.get('alpha', None)
        self.k_gamma                   = kwargs.get('num_categories', 4)
        self.pinv                      = kwargs.get('pinv', 0.)                 # If > 0, this will be the last entry in self.rate_probs, and 0 will be the last entry in self.rate_factors
        self._save_custom_matrix_freqs = kwargs.get('save_custom_frequencies', "custom_matrix_frequencies.txt")
        self.neutral_scaling           = kwargs.get('neutral_scaling', False)
        self.code                      = None
        
        # There are lots of these
        self.aa_models    = ['jtt', 'wag', 'lg', 'ab', 'mtmam', 'mtrev24', 'dayhoff']
        
        self._check_acceptable_model()
        self._check_hetcodon_model()
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
                raise ValueError("\n\nUnknown code to evolve.")



    def _check_acceptable_model(self):
        '''
            Check that the model type has been properly specified. If the model is custom and there is a code, they must be the same length.
        '''
        
        self.model_type = self.model_type.replace("94", "") # Allows users to give GY94, MG94 
        assert(type(self.neutral_scaling) is bool), "\n\nThe argument 'neutral_scaling' must be True or False."
        accepted_models = ['nucleotide', 'codon', 'gy', 'mg', 'mutsel', 'ecm', 'ecmrest', 'ecmunrest', 'custom'] + self.aa_models
        assert(self.model_type in accepted_models), "\n\nInappropriate model type specified."
        assert(type(self.params) is dict), "\n\nThe parameters argument must be a dictionary."
        
        # Assign default codon, ecm models
        if self.model_type == 'codon':
            self.model_type = 'gy'
            print("Using default codon model, GY-style.")
        if self.model_type == 'ecm':
            self.model_type = 'ecmrest'
            print("Using restricted ECM model.")
        if self.model_type == 'custom':
            assert("matrix" in self.params), "\n\nTo use a custom model, you must provide a matrix in your params dictionary under the key 'matrix'. Also note that pyvolve orders nucleotides, amino acids, and codons alphabetically by their abbreviations (e.g. amino acids are ordered A, C, D, ... Y)."          
            if "state_freqs" in self.params:
                # Matrix must be symmetric, dimensions must be compatible, and frequencies should sum to 1
                assert( abs(1. - np.sum(self.params["state_freqs"])) <= ZERO), "\n\nProvided state frequencies for custom model do not sum to 1."
                assert(np.all( self.params["matrix"] == self.params["matrix"].T)),"\n\nYou have provided a nonsymmetric matrix and state frequencies for your custom model. If you wish to provide state frequencies, your matrix must be symmetric."
                assert(len(self.params["state_freqs"]) == len(self.params["matrix"][0])),"\n\nThe dimensions of your provided state frequencies do not match the dimensions of your custom matrix."
            else:
                print("\nSince you have provided a custom matrix without state frequencies, pyvolve will calculate them for you directly from your provided custom matrix. These frequencies will be saved, for your convenience, to a file",self._save_custom_matrix_freqs,".")
        


    def _check_hetcodon_model(self):
        '''
            Determine if this is a heterogenous codon model and assign self.hetcodon_model accordingly.
        '''
        self.hetcodon_model = False
        # This must be done here in order to check the type of model. This code is also used in the sanity check, though.
        if "omega" in self.params:
            self.params["beta"] = self.params["omega"]
            self.params.pop("omega")
            
        # If iterable, codon model.
        try:
            (x for x in self.params["beta"])
            self.hetcodon_model = True
        except:
            self.hetcodon_model = False



    def _construct_model(self):
        '''
            Construct the evolutionary model.
        '''
        
        # Assign rate heterogeneity
        if self.hetcodon_model:
            self._assign_rate_probs()
        else:
            self._assign_rates()
        
        # Assign matrix (or matrices, for a heterogeneous codon model)
        self._assign_matrix() 

        # Finally, once all matrices are built, assign code (custom, nuc, aa, or codon)
        self._assign_code()



    def _assign_matrix(self):
        '''
            Construct the model rate matrix, Q, based on model_type by calling the matrix_builder module. Alternatively, call the method self._assign_codon_model_matrices() if we have a heterogenous codon model.
            Note that before matrix construction, we sanity check and update, as needed, all provided parameters.
        '''
        
        
        if self.model_type == 'nucleotide':
            self.params = Nucleotide_Sanity(self.model_type, self.params, size = 4)()
            self.matrix = Nucleotide_Matrix(self.model_type, self.params)()
                
                    
        elif self.model_type in self.aa_models:
            self.params = AminoAcid_Sanity(self.model_type, self.params, size = 20)()
            self.matrix = AminoAcid_Matrix(self.model_type, self.params)()
             
             
        elif self.model_type == 'gy' or self.model_type == 'mg':
            self.params = MechCodon_Sanity(self.model_type, self.params, size = 61, hetcodon_model = self.hetcodon_model )()
            self.params["neutral_scaling"] = self.neutral_scaling
            if self.hetcodon_model:
                self._assign_hetcodon_model_matrices()
            else:
                self.matrix = MechCodon_Matrix(self.model_type, self.params )()
        
        
        elif 'ecm' in self.model_type:
            self.params = ECM_Sanity(self.model_type, self.params, size = 61)()
            self.matrix = ECM_Matrix(self.model_type, self.params)()
 
 
        elif self.model_type == 'mutsel':
            self.params = MutSel_Sanity(self.model_type, self.params)()
            self.matrix = MutSel_Matrix(self.model_type, self.params)()
            
            # Need to construct and add frequencies to the model dictionary if the matrix was built with fitness values
            if not self.params["calc_by_freqs"]:
                self._calculate_state_freqs_from_matrix()
        
        
        elif self.model_type == 'custom':
            self._assign_custom_matrix()
            np.savetxt(self._save_custom_matrix_freqs, self.params["state_freqs"]) 

        else:
            raise ValueError("\n\nYou have reached this in error! Please file a bug report, with this error, at https://github.com/sjspielman/pyvolve/issues .")

        # Double check that state frequencies made it in. 
        assert("state_freqs" in self.params), "\n\nYour model has no state frequencies."
            


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
        assert ( np.allclose( np.zeros(dim), np.sum(custom_matrix, 1) , rtol=1e-5) ), "\n\nRows in custom transition matrix do not sum to 0."
        
        # "Re-normalize" matrix with higher tolerance and confirm, re-save
        for s in range(dim):
            temp_sum = np.sum(custom_matrix[s]) - np.sum(custom_matrix[s][s])
            custom_matrix[s][s] = -1. * temp_sum
            assert ( abs(np.sum(custom_matrix[s])) <= ZERO ), "\n\nRe-normalized row in custom transition matrix does not sum to 0."
        self.matrix = custom_matrix
        
        # Multiply by frequencies or extract frequencies if none were provided
        if "state_freqs" in self.params:
            diag_state_freqs = np.diag(self.params["state_freqs"])
            self.matrix = np.dot(diag_state_freqs, self.matrix) 
            # One more check..          
            assert ( abs(np.sum(self.matrix[s])) <= ZERO ), "\n\nAfter frequency multiplication, rows in custom transition matrix no longer sum to 0."
        else:
            self._calculate_state_freqs_from_matrix()
        
      
      
      
      
    def _assign_hetcodon_model_matrices(self):
        '''
            Construct each model rate matrix, Q, to create a list of codon-model matrices.
        '''
        # Determine mean dN/dS for scaling calculation. Note that these calculations will be effectively ignored if neutral scaling has been specified.
        dnds_values = np.array(self.params["beta"]) / np.array(self.params["alpha"])
        self.params["hetmodel_mean_dnds"] = np.average( dnds_values, weights = self.rate_probs )

        # Construct matrices
        self.matrix = []
        for i in range(len( self.params['beta'] )):
            temp_params = deepcopy(self.params)
            temp_params['beta'] = self.params['beta'][i]
            temp_params['alpha'] = self.params['alpha'][i]
            mb = MechCodon_Matrix(self.model_type, temp_params)
            self.matrix.append( mb() )
        assert(len(self.matrix) > 0), "Matrices for a heterogeneous codon model were improperly constructed."




    def _calculate_state_freqs_from_matrix(self):
        '''
            Determine the vector of state frequencies numerically from a matrix.
            This method is used when a MutSel model was built up using fitness values, or when a custom matrix was specified from which state frequencies must be calculated.
        ''' 
        size = self.matrix.shape[0]
        
        # Find maximum eigenvalue, index
        (w, v) = linalg.eig(self.matrix, left=True, right=False)
        max_i = np.argmax(w)
        max_w = w[max_i]
        assert( abs(max_w) <= ZERO ), "\n\nCould not extract dominant eigenvalue from matrix to determine state frequencies."
        max_v = v[:,max_i]
        max_v /= np.sum(max_v)
        eq_freqs = max_v.real # these are the stationary frequencies
        
        # Equaling zero gets numerically horrible, so clean and re-normalize.
        eq_freqs[eq_freqs == 0.] = ZERO
        #eq_freqs /= np.sum(eq_freqs) 
        assert(abs(1. - np.sum(eq_freqs)) <= ZERO), "\n\nState frequencies calculated calculated from matrix do not sum to 1."
        
        # Some sanity checks, many of which are overkill. 
        assert np.allclose(np.zeros(size), np.dot(eq_freqs, self.matrix)), "State frequencies not properly calculated." # should be true since eigenvalue of zero
        pi_inv = np.diag(1.0 / eq_freqs)
        s = np.dot(self.matrix, pi_inv)
        assert np.allclose(self.matrix, np.dot(s, np.diag(eq_freqs)), atol=ZERO, rtol=1e-5), "\n\nMatrix cannot be recovered from exchangeability and equilibrium when computing state frequencies from matrix."
        assert(not np.allclose(eq_freqs, np.zeros(size))), "\n\nState frequencies were not calculated from matrix at all."
        
        # Finally, assign the frequencies
        self.params["state_freqs"] = eq_freqs







    def _assign_rates(self):
        '''
            Assign and sanity-check site heterogeneity rates for non-codon models. Note that heterogeneity is *ignored* for ECM models, because it's entirely unclear how/why this should be implemented.
        '''
        if "ECM" in self.model_type:
            # Disallow heterogeneity for ECM
            self.rate_probs = np.ones(1)
        else:
            # Draw gamma rates if specified
            if self.alpha is not None:
                assert(self.pinv >= 0. and self.pinv <= 1.), "\n\nThe proportion of invariant sites must be a value between 0 and 1 (inclusive)."
                self._draw_gamma_rates() # Sanity done inside
            else:
                self._assign_rate_probs()
                self._sanity_rate_factors()            



 
    def _draw_gamma_rates(self):
        '''
            Function to draw and assign rates from a discretized gamma distribution, if specified. By default, 4 categories are drawn.
        '''       
        if self.rate_probs is not None:
            print("\nThe provided value for the `rate_probs` argument will be ignored since gamma-distributed heterogeneity has been specified with the alpha parameter.")        
        if type(self.k_gamma) is not int:
            raise TypeError("\nProvided argument `num_categories` must be an integer.")

        #### Note that this code is adapted from gamma.c in PAML ####
        rv = gamma(self.alpha, scale = 1./self.alpha)
        freqK = np.zeros(self.k_gamma)
        rK = np.zeros(self.k_gamma)

        for i in range(self.k_gamma-1):
            raw=rv.ppf( (i+1.)/self.k_gamma )
            freqK[i] = gammainc(self.alpha + 1, raw*self.alpha)

        rK[0] = freqK[0] * self.k_gamma
        rK[self.k_gamma-1] = (1-freqK[self.k_gamma-2]) * self.k_gamma
        for i in range(1,self.k_gamma-1):
            rK[i] = self.k_gamma * (freqK[i] -freqK[i-1])    
        #############################################################
        
        # Kindly note, only the gamma-distributed rates satify \sum(p*r) = 1, because +I really makes no sense as a model. Yet, here we are, by MS-reviewer demand... 
        if self.pinv <= ZERO:
            self.rate_probs = np.repeat(1./self.k_gamma, self.k_gamma)
            self.rate_factors = deepcopy(rK)
        else:
            freqK *= (1. - self.pinv)
            freqK = list(freqK)
            freqK.append(self.pinv)
            self.rate_probs = np.array(freqK)
            
            rK = list(rK)
            rK.append(0.)
            self.rate_factors = np.array(rK)            
        



# 
# 
#     def _assign_gamma_pinv_rate_probs(self):
#         '''
#             Function to compute rate probabilities under the specific Gamma + Pinv rate heterogeneity scheme.
#         '''
#         # Default
#         if self.rate_probs is None:
#             remaining_prob = 1. - self.pinv
#             self.rate_probs = list(np.repeat( remaining_prob/self.k_gamma, self.k_gamma ))        
#         # Custom
#         else:
#             self.rate_probs = list(self.rate_probs)
#         self.rate_probs.append(self.pinv)
#         self.rate_probs = np.array( self.rate_probs ) 
#             
       
     
         
    def _assign_rate_probs(self):
        '''
            Compute probabilities for rate heterogeneity.
            By default, equal probabilities will be assigned to all rate categories (either scaling factors or dN/dS values and corresponding matrices for a heterogeneous codon model).
        '''
        
        # Gamma + Pinv
        #if self.alpha is not None and self.pinv > 0.:
        #    self._assign_gamma_pinv_rate_probs()
        
        if self.hetcodon_model:
            num_probs = len(self.params["beta"])
        else:
            num_probs = len(self.rate_factors)
        
        # Nothing provided, default.
        if self.rate_probs is None:
            self.rate_probs = np.repeat( 1./num_probs, num_probs )  

        ### Perform some checks ###
        
        # Ensure sums to 1
        assert(abs(1. - np.sum(self.rate_probs)) <= ZERO), "\n\nProvided rate probabilities (rate_probs list) must sum to 1.\nNote: if you are specifying Gamma+Pinv heteregeneity with custom probabilities, ensure that the sum of pinv and your rate_probs list is equal to 1."
            
        # Ensure numpy array
        try:
            self.rate_probs = np.array( self.rate_probs )
        except:
            raise TypeError("\n\nRate probabilities improperly specified. You must supply either a list of numpy array of probabilities.")                
            
        # Size sanity check.
        assert( len(self.rate_probs) == num_probs ), "\n\nDifferent numbers of rates and associated probabilities."
            



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
            raise TypeError("\n Rate factors improperly specified. You must supply either a list of numpy array of probabilities.")      

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
        
        
        
    def is_hetcodon_model(self):
        '''
            Return True if the model is a heterogeneous codon model and return False otherwise.
        '''
        return self.hetcodon_model
        
 
    # Convenience functions for users to call up parameters easily. #
        
    def extract_mutation_rates(self):
        '''
            Convenience function for returning the mutation rate dictionary to users.
        '''
        try:
            return self.params["mu"]
        except:
            print("\nYour model does not include mutation rates, so you can't extract them.")
    
    
    def extract_rate_matrix(self):
        '''
            Convenience function for returning the rate matrix/matrices users.
        '''
        return self.matrix


    def extract_state_freqs(self):
        '''
            Convenience function for returning the stationary frequencies.
        '''
        return self.params["state_freqs"]
    
    
    def extract_parameters(self):
        '''
            Convenience function for returning the params dictionary, which contains all model parameters used to construct the rate matrix (except for nucleotide/amino-acid rate heterogeneity).
        '''
        return self.params

