#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    Define Partition() class.
'''

from model import * 



class Partition():
    '''
        Defines a Partition() object.
    '''
    def __init__(self, **kwargs):
        self.size              = kwargs.get('size', [])   # List of integers representing partition length. If there is no rate heterogeneity, then the list is length 1. Else, list is length k, where k is the number of rate categories.
        self.models            = kwargs.get('models', None)  # List of models associated with this partition. When length 1, temporally homogeneous.
        self.root_model_name   = None  # NAME of of Model beginning evolution at root of tree. Used under *branch heterogeneity*, and should be None or False if process is temporally homogeneous. If there is branch heterogeneity, this string *MUST* correspond to one of the Model() object's names and also a corresponding phylogeny flag.
        self.root_seq          = None  # User may choose to provide a root sequence for each partition, and it'll be stored here. Totally optional - will otherwise be generated from steady-state frequencies.
        self.shuffle           = False # Shuffle sites after evolving? 
        self._root_model       = None  # The actual root model. Used internally in evolver module.

    
    def branch_het(self):
        ''' 
            Return False if there is no branch heterogeneity, True if there is.
        '''
        if isinstance(self.models, Model) or isinstance(self.models, CodonModel) or len(self.models) == 1:
            return False
        elif len(self.models) > 1:
            return True
        else:
            raise AssertionError("\n\nPartition has no associated models, so I have no clue what sort of heterogeneity there is...because there's no model...")

    def site_het(self):
        ''' 
            Return False if there is no site heterogeneity, True if there is.
        '''
        m = self.models[0]
        if m.num_classes() > 1:
            return True
        else:
            return False
    
    
    def codon_model(self):
        '''
            Return True if Partition is evolving with CodonModel objects (dN/dS model heterogeneity) and False otherwise.
        '''
        if isinstance(self.models[0], CodonModel):
            return True
        elif isinstance(self.models[0], Model):
            return False
        else:
            raise AssertionError("\n\nPartition has no models so can't tell if codonmodel or not...")
