#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################




'''
    This module performs sanity checks on model parameters, specifically those which will be used to build the rate matrix.
'''



import numpy as np
from scipy import linalg
from copy import deepcopy
from .genetics import *
from .state_freqs import *
import warnings
ZERO      = 1e-8
MOLECULES = Genetics()


class ParametersSanity(object):
    '''
        Module used to sanity check Model parameters and update them accordingly.
        
        Child class include the following:
            1. *AminoAcid_Sanity* 
                - Empirical amino acid models
            2. *Nucleotide_Sanity* 
                - Nucleotide models (GTR and nested)
            3. *MG_Sanity*
                - MG-style models
            4. *GY_Sanity*
                - GY-style models
            5. *MutSel_Sanity* 
                - Mutation-selection model (Halpern and Bruno 1998), extended for either codon or nucleotides
            6. *ECM_Sanity*
                - ECM (Kosiol2007) empirical codon model
        
    '''
    
    def __init__(self, model_type, parameters, **kwargs):
        '''
            Requires two positional argument:
                1. **model_type**, the type of model which will be built
                2. **parameters**, a dictionary containing parameters about the substitution process which will be checked.

            
            And one optional arguments:
                1. **size**, optional argument indicating the length of the code used (4,20,61). If not provided, it will be figured out..
        '''
        self.model_type = model_type        
        self.params     = parameters
        self.size       = kwargs.get("size", None)
        
        

    def __call__(self):
        '''
            Call method. Performs sanity check and returns the updated parameters.
        '''
        self.sanity_params()
        return self.params
        
        
    
    def _sanity_state_freqs_common(self):
        '''
            State frequency sanity checks common to all child classes.
        '''
        assert( len(self.params['state_freqs']) == self.size ), "\n\nThe value associated with the 'state_freqs' key in the provided parameters dictionary does not contain the correct number of values for your specified model."
        assert( abs(1. - np.sum(self.params['state_freqs']) <= ZERO) ), "\n\nProvided state frequencies do not sum to 1."



    def _sanity_state_freqs(self, empirical = False):
        '''
            Sanity-check specifically state_freqs key/value in the params dictionary.
            If state_freqs not provided, then they are either set to equal, or if the model is empirical, then they are set to the default empirical model's frequencies.
        '''
        # Fill in if missing
        if 'state_freqs' not in self.params:
            if not empirical:
                self.params['state_freqs'] = np.repeat(1./self.size, self.size)
            else:
                f = EmpiricalModelFrequencies(self.model_type)
                self.params['state_freqs'] = f.compute_frequencies() 
                 
        # If present, check the size and the sum
        else:
            self._sanity_state_freqs_common()

           

    def _sanity_mutation_rates(self):
        '''
            Sanity-check specifically mu key/value in params dictionary.
            NOTE: mutation-selection models may take asymmetric mutation rates. However, this function assumes that rates are symmetric.
            Thus, if A->C is present but C->A is not, the latter will take on the A->C value.
            Any missing mutation rates are given a value of 1.
        '''
        
        if "mu" in self.params:
            try:        
                # Dictionary of mu's provided. Make sure dictionary is full.
                for key in ['AC', 'AG', 'AT', 'CG', 'CT', 'GT']:
                    rev_key = str(key[1] + key[0])
                
                    # Neither key pair. Add both
                    if key not in self.params['mu'] and rev_key not in self.params['mu']:
                        self.params['mu'][key] = 1.
                        self.params['mu'][rev_key] = 1.
                
                    # If one key but not the other, fill in missing one symmetrically.
                    elif key not in self.params['mu'] and rev_key in self.params['mu']:
                        self.params['mu'][key] = self.params['mu'][rev_key]
                
                    elif key in self.params['mu'] and rev_key not in self.params['mu']:
                        self.params['mu'][rev_key] = self.params['mu'][key]
            except:
                raise KeyError("\n\nYou must provide a dictionary of mutation rates for nucleotide pairs to the key 'mu' in the model parameters dictionary. Alternatively, you may provide the key 'kappa' to specify only a TI/TV ratio.")
        
        # Nothing specified, so simply use equal mutation rates to construct matrix
        else:
            self.params['mu'] = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}

        # Apply kappa as needed.        
        if 'kappa' in self.params:
            temp_mu = deepcopy( self.params['mu'] )
            self.params['mu'] = {'AC': temp_mu['AC'], 'AG': temp_mu['AG'] * float(self.params['kappa']), 'AT': temp_mu['AT'], 'CG': temp_mu['CG'], 'CT': temp_mu['CT']*float(self.params['kappa']), 'GT': temp_mu['GT'], 'CA': temp_mu['CA'], 'GA': temp_mu['GA'] * float(self.params['kappa']), 'TA': temp_mu['TA'], 'GC': temp_mu['GC'], 'TC': temp_mu['TC']*float(self.params['kappa']), 'TG': temp_mu['TG']}
       



    def _sanity_dnds(self):
        if 'omega' in self.params:
            self.params['beta'] = self.params['omega']
            self.params.pop("omega")
        
        if 'beta' not in self.params:
            raise KeyError("\n\nYou must provide a dN value (using either the key 'beta' or 'omega') in your dictionary of model parameters to use this model.")
    
        
        # heterogeneous codon models
        if self.hetmodel:
            try:
                lenalpha = len(self.params["alpha"])
            except KeyError:
                self.params['alpha'] = np.ones(len(self.params["beta"])) 
            except TypeError:
                raise TypeError("To specify both dN and dS heterogeneity, provide lists (or numpy arrays), of the *same lengths*, for keys 'alpha' and 'beta'.\n To only specify heterogeneity in dN, use only the key 'beta' (or 'omega').")
            
        else:
            if "alpha" not in self.params:
                self.params['alpha'] = 1.








class AminoAcid_Sanity(ParametersSanity):
    ''' 
        Child class of ParametersSanity for amino acid models.
        Required checks:
            1. state frequencies
    '''        
    
    def __init__(self, *args, **kwargs):
        super(AminoAcid_Sanity, self).__init__(*args, **kwargs)
        
    
    def sanity_params(self):
        '''
            Perform required sanity checks.
        '''
        self._sanity_state_freqs(empirical = True)
        





class Nucleotide_Sanity(ParametersSanity):
    ''' 
        Child class of ParametersSanity for nucleotide models.
        Required checks:
            1. state frequencies
            2. mutation rates
    '''        
    
    def __init__(self, *args, **kwargs):
        super(Nucleotide_Sanity, self).__init__(*args, **kwargs)

 
    def sanity_params(self):
        '''
            Perform required sanity checks.
        '''
        self._sanity_state_freqs()
        self._sanity_mutation_rates()
      
      
      





class MechCodon_Sanity(ParametersSanity):
    ''' 
        Child class of ParametersSanity for MG- and GY-style codon models.
        Required checks:
            1. state frequencies, and also nucleotide frequencies for MG
            3. mutation rates
            4. omega (or beta, alpha)
    '''        
    
    def __init__(self, *args, **kwargs):
        self.hetmodel = kwargs.get("hetcodon_model", False)
        super(MechCodon_Sanity, self).__init__(*args, **kwargs)

     
    def sanity_params(self):
        '''
            Perform required sanity checks.
        '''
        if self.model_type == "mg":
            self._sanity_state_freqs_MG()
        else:   
            self._sanity_state_freqs()
        self._sanity_mutation_rates()
        self._sanity_dnds()



    def _sanity_state_freqs_MG(self):
        '''
            This checks state frequencies specifically for MG-style models, as both nucleotide and codon frequencies are required.
            Additionally considers the key 'nuc_freqs' which may be provided for MG-style models. If provided, then state frequencies are computed as F1x4.
            Otherwise, nuc_freqs are derived from provided state_freqs.
        '''
        if 'nuc_freqs' in self.params:
            assert(len(self.params['nuc_freqs']) == 4), "\n\nTo provide nucleotide frequencies for MG-style models, be sure to provide *4* values."
            
            if 'state_freqs' in self.params:
                warn("You have provided both 'nuc_freqs' and 'state_freqs' for your MG-style model. The codon frequencies ('state_freqs') will be overwritten with the correct stationary frequencies (F1x4).")
            self.params["state_freqs"] = self._calc_f1x4()
            
        else:
            self._sanity_state_freqs() # basic checks, and below extract the nucleotide frequencies from it
            f = CustomFrequencies(by = 'codon', freq_dict = dict(list(zip(MOLECULES.codons, self.params['state_freqs']))))
            self.params["nuc_freqs"] = f.compute_frequencies(type = 'nucleotide')


    def _calc_f1x4(self):
        '''
            Calculate state codon frequencies from nucleotide frequencies for use in an MG-style model.
            NOTE: F1x4 are calculated because these *are* the state frequencies for an MG-stye model.
        '''
        f1x4 = np.ones(61)
        nf = self.params["nuc_freqs"]
        assert( abs(np.sum(nf) - 1.) <= ZERO), "\n\nProvided nucleotide frequencies for an MG-style model do not sum to 1."
        pi_stop = (nf[3]*nf[0]*nf[2]) + (nf[3]*nf[2]*nf[0]) + (nf[3]*nf[0]*nf[0])
        for i in range(61):
            codon = MOLECULES.codons[i]
            for j in range(3):
                f1x4[i] *= nf[ MOLECULES.nucleotides.index(codon[j]) ]
        f1x4 /= (1. - pi_stop)
        assert( abs(np.sum(f1x4) - 1.) <= ZERO ), "\n\nCould not properly calculate F1x4 frequencies for an MG-style model."
        return f1x4 
















        



class ECM_Sanity(ParametersSanity):
    ''' 
        Child class of ParametersSanity for the empirical codon model (ECM).
        Required checks:
            1. state frequencies 
            3. "special" ECM parameters. See Kosiol et al. 2007 for details.
    '''        
    
    def __init__(self, *args, **kwargs):
        super(ECM_Sanity, self).__init__(*args, **kwargs)

 
    def sanity_params(self):
        '''
            Perform required sanity checks.
        '''
        self._sanity_state_freqs(empirical = True)
        
        if 'omega' in self.params:
            self.params['beta'] = self.params['omega']
            self.params.pop('omega')
        if 'beta' not in self.params:
            self.params['beta'] = 1.
        if 'alpha' not in self.params:
            self.params['alpha'] = 1.
        if 'k_ti' not in self.params:
            self.params['k_ti'] = 1.
        if 'k_tv' not in self.params:
            self.params['k_tv'] = 1.






class MutSel_Sanity(ParametersSanity):
    ''' 
        Child class of ParametersSanity for the empirical codon model (ECM).
        Required checks:
            1. state frequencies and/or fitness values
            2. mutation rates
    '''        
    
    def __init__(self, *args, **kwargs):
        super(MutSel_Sanity, self).__init__(*args, **kwargs)


    def sanity_params(self):
        '''
            Perform required sanity checks.
        '''
        self._sanity_state_freqs_fitness()
        self._sanity_mutation_rates()
      


    def _sanity_state_freqs_fitness(self):
        '''
            Check that either state_freqs or fitness have been properly provided for a MutSel model.
            Additionally, add these keys to the parameters dictionary:
                1. "codon_model". Boolean indicating if this is a "codon" (1) or "nucleotide" (0) MutSel model
                2. "calc_by_freqs". Boolean indicating if calculations will be done using "state_freqs" (1) or "fitness" (0) values
        '''
        self.params["codon_model"] = None
        self.params["calc_by_freqs"] = None
        
        if 'state_freqs' in self.params:
            self.params["calc_by_freqs"] = True
            
            if len(self.params['state_freqs']) == 61:
                self.size = 61
                self.params["codon_model"] = True
                
            elif len(self.params['state_freqs']) == 4:
                self.size = 4
                self.params["codon_model"] = False
                
            else:    
                raise ValueError("\n\n The list of stationary frequencies for a mutation-selection model must be of length 4 (nucleotides) or 61 (codons).")
            self._sanity_state_freqs()
        
        elif 'fitness' in self.params:
            self.params["calc_by_freqs"] = False
            
            if len(self.params['fitness']) == 61 or len(self.params['fitness']) == 20:
                self.params["codon_model"] = True
                
                # Replace length-20 fitness with length-61 fitness, assuming equal fitness for synonymous codons.
                if len(self.params['fitness']) == 20: 
                    self._amino_to_codon_fitness()
            
            elif len(self.params['fitness']) == 4:
                self.params["codon_model"] = False
            
            else:
                raise ValueError("\n\n Your provided fitness values should be in a vector of length 4, 20, or 61.")
        
        else:
            raise KeyError("\n\nYou must provide either state frequencies ('state_freqs') or fitness ('fitness') as parameters for a mutation-selection model.")
        
        assert(self.params["codon_model"] is not None), "\n\nUnable to determine if your mutation-selection model uses nucleotides or codons."
        assert(self.params["calc_by_freqs"] is not None), "\n\nUnable to determine if you specified frequencies or fitnesses for your mutation-selection model."
              
        


    def _amino_to_codon_fitness(self):
        '''
            Convert a vector of amino acid fitness values to codon fitness values, assuming equal fitness among synonymous codons.
        '''
        d = {}
        for i in range(20):
            syn_codons = MOLECULES.genetic_code[i]
            for syn in syn_codons:
                d[ syn ] = self.params['fitness'][i]   
        
        codon_fitness = np.zeros(61)
        count = 0
        for i in range(61):
            codon_fitness[i] = d[MOLECULES.codons[i]]
            
        self.params['fitness'] = codon_fitness
        
