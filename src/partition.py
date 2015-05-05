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
                
                1. **size**, integer giving the root length of this partition
                2. **models**, either a single Model/CodonModel instance (for cases of branch homogeneity), or a list of Model/CodonModel instances (for cases of branch heterogeneity)
        
            Examples:
                .. code-block:: python

                   >>> # Define a temporally homogeneous partition
                   >>> my_partition = Partition(size = 500, my_model) 
                   
                   >>> # Define a temporally heterogeneous partition, in which three models (model1, model2, and model3) are used during sequence evolution.
                   >>> my_other_partition = Partition(size = 134, [model1, model2, model3])       
                
        '''
                
                
        self.size              = kwargs.get('size', [])   # List of integers representing partition length. If there is no rate heterogeneity, then the list is length 1. Else, list is length k, where k is the number of rate categories.
        self.models            = kwargs.get('models', None)  # List of models associated with this partition. When length 1, temporally homogeneous.
        self.root_model_name   = None  # NAME of of Model beginning evolution at root of tree. Used under *branch heterogeneity*, and should be None or False if process is temporally homogeneous. If there is branch heterogeneity, this string *MUST* correspond to one of the Model() object's names and also a corresponding phylogeny flag.
        self.shuffle           = False # Shuffle sites after evolving? The evolver module will set this up.
        self._root_model       = None  # The actual root model. Used internally in evolver module.

    
    def branch_het(self):
        ''' 
            Return True if the partition uses branch heterogeneity, and False if homogeneous.
        '''
        if isinstance(self.models, Model) or isinstance(self.models, CodonModel) or len(self.models) == 1:
            return False
        elif len(self.models) > 1:
            return True
        else:
            raise AssertionError("\n\nPartition has no associated models, so I have no clue what sort of heterogeneity there is...because there's no model...")

    def site_het(self):
        ''' 
            Return True if the partition uses site/rate heterogeneity, and False if homogeneous.
        '''
        m = self.models[0]
        if m.num_classes() > 1:
            return True
        else:
            return False
    
    
    def codon_model(self):
        '''
            Return True if Partition is evolving according to CodonModel objects (dN/dS model heterogeneity), and False if Model objects.
        '''
        if isinstance(self.models[0], CodonModel):
            return True
        elif isinstance(self.models[0], Model):
            return False
        else:
            raise AssertionError("\n\nPartition has no models so can't tell if codonmodel or not...")
