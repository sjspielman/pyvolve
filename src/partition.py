#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    This module defines the Partition() class, which indicates a particular evolutionary unit. 
'''

from .model import * 

class Partition():

    def __init__(self, **kwargs):
        '''
            Required keyword arguments:
                
                1. **size**, integer giving the root length of this partition (only required when no root sequence is given)
                2. **models** (or **model**), either a single Model object (for cases of branch homogeneity), or a list of Model objects (for cases of branch heterogeneity).
        
            Optional keyword arguments:
            
                1. **root_sequence**, a string giving the ancestral sequence for this partition. Note that, when provided, the **size** argument is not needed.
                2. **root_model_name**, the *name attribute* of the model to be used at the root of the phylogeny. Applicable only to cases of *branch heterogeneity*.

            Examples:
                .. code-block:: python

                   >>> # Define a temporally homogeneous partition
                   >>> my_partition = Partition(models = my_model, size = 500) 
                   
                   >>> # Define a temporally heterogeneous partition, in which three models (model1, model2, and rootmodel) are used during sequence evolution, and rootmodel is the model at the root of the tree
                   >>> my_other_partition = Partition(models = [model1, model2, rootmodel], size = 134, root_model_name = rootmodel.name)       
                
        '''
                
                
        self.size              = kwargs.get('size', None)   # Will be converted to list of integers representing partition length. If there is no rate heterogeneity, then the list is length 1. Else, list is length k, where k is the number of rate categories.
        self.MRCA              = kwargs.get('root_sequence', None) # String of root sequence. If provided, all specified rate heterogeneity and the size argument *will be ignored*.
        self.models            = kwargs.get('models', None)  # List of models associated with this partition. When length 1 (or not provided as a list) temporally homogeneous.
        if self.models is None:
            self.model         = kwargs.get('model', None)
        self.root_model_name   = kwargs.get('root_model_name', None)  # NAME of Model beginning evolution at root of tree. Used under *branch heterogeneity*, and should be None or False if process is temporally homogeneous. If there is branch heterogeneity, this string *MUST* correspond to one of the Model() object's names.
        self._shuffle          = False # Shuffle sites after evolving?
        self._root_model       = None  # The actual root model object.

        self._partition_sanity()
        self._size_MRCA_sanity()
        if self.MRCA is None:
            self._divvy_partition_size()



    def _size_MRCA_sanity(self):
        '''
            Sanity checks and setup for the size and MRCA, if provided.
        '''
        assert(self.size is not None or self.MRCA is not None), "\n\nWhen defining a Partition object, you must specify either a root sequence or a partition size."
        
        if self.MRCA is not None:
            assert(type(self.MRCA) is str), "\n\nThe provided root sequence in your Partition object must be a string."
            if self.size is not None:
                print("\n\nWARNING: You provided both a size and a root sequence for your Partition. The size argument will be ignored.")
            code_step = len(self._root_model.code[0])
            self.size = [len(self.MRCA) / code_step]
        
            # Remove site-rate heterogeneity if MRCA was provided
            for model in self.models:
                model.rate_probs = np.array([1.])
                model.rate_factors = np.array([1.])
            



    def _partition_sanity(self):
        ''' 
            Sanity checks that Partition has been properly setup.
        '''
        assert(self.models is not None), "\n\nNo model(s) was/were provided to this Partition. Please check that you have specified a proper model object, or list of model objects, using the keyword 'model' or 'models' (both are accepted)."
        # Ensure that self.models is a list
        if type(self.models) is not list:
            self.models = [self.models]
        
        # Ensure that all models use the same code
        code1 = self.models[0].code
        if len(self.models) > 1:
            for m in self.models[1:]:
                assert(m.code == code1), "\n\nYour partitions are evolving according to different codes/alphabets. This is not allowed."

        # Assign _root_model
        if self.branch_het():
            for m in self.models:
                if m.name == self.root_model_name:
                    self._root_model = m
                    break
        else:
            self._root_model = self.models[0] 
        assert(self._root_model != None), "\n Root model not properly assigned in your partition. Make sure that you specified a root model name if you have branch heterogeneity! Do so with the argument root_model_name."

        # Ensure branch-site is ok - number of rate categories has to be the same across branches.
        if self.site_het():
            self._shuffle = True
            for model in self.models:
                assert( len(model.rate_probs) == len(self._root_model.rate_probs) ), "For branch-site models, the number of rate categories must remain constant over the tree in a given partition."




    def _divvy_partition_size(self):
        '''
            Turn size attribute into a list of different rate-heterogeneity size chunks (based on rate_probs to model object).
            Note that this is *probability-based!* So if you specified 25% for a category, it is not strictly true that 25% of sites will be in that category, but rather that category will occur with a probability of 25%. 
            
            If no rate heterogeneity, will simply be a list of length 1 containing full size.
        '''
        nc = self._root_model.num_classes()
        rate_occurrences = np.random.choice(nc, int(self.size), p = self._root_model.rate_probs)          
        new_size = np.bincount(rate_occurrences, minlength = nc)
        assert( sum(new_size) ==  self.size ), "\n\nImproperly divvied up rate heterogeneity."
        assert( len(new_size) == nc), "\n\nPartition size does not correspond to the number of rate categories. Please report this error."
        self.size = list(new_size)

    
    def branch_het(self):
        ''' 
            Return True if the partition uses branch heterogeneity, and False if homogeneous.
        '''

        if isinstance(self.models, Model) or len(self.models) == 1:
            return False
        elif len(self.models) > 1:
            return True
        else:
            raise AssertionError("\n\nPartition has no associated models.")



    def site_het(self):
        ''' 
            Return True if the partition uses site/rate heterogeneity, and False if homogeneous.
            Also returns False if an MRCA has been provided, as we do not support rate heterogeneity in this case.
        '''
        if self.MRCA is not None:
            return False
        if self.models[0].num_classes() > 1:
            return True
        else:
            return False
    
    
    def is_codon_model(self):
        '''
            Return True if the partition is evolving with dN/dS heterogeneity, and False otherwise.
        '''
        return self.models[0].is_codon_model()

