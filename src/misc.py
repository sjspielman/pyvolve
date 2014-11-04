#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Class definitions and global variables used throughout pyvolve.
'''


ZERO = 1e-8

class Genetics():
    '''
        Molecular alphabet objects.
    '''
    
    def __init__(self):
        self.pyrims       = ["C", "T"]
        self.purines      = ["A", "G"]
        self.nucleotides  = ["A", "C", "G", "T"]
        self.amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        self.genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
        self.codon_dict   = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
        self.codons       = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.stop_codons  = ["TAA", "TAG", "TGA"]

class Tree():
    '''
        Defines a Tree() object. The final tree contains a series of nested Tree() objects.
    '''
    def __init__(self):
        self.name           = None # Internal node unique id or leaf name
        self.children       = []   # List of children, each of which is a Tree() object itself. If len(children) == 0, this tree is a tip.
        self.branch_length  = None # Branch length leading up to node
        self.model_flag     = None # Flag indicate that this branch evolves according to a distinct model from parent
        self.seq            = None # Contains sequence (represented by integers) for a given node. EVENTUALLY THIS WILL BE REPLACED BY A LIST OF Site() OBJECTS.



class EvoModels(object):
    ''' 
        Parent class for Model(), CodonModel()
    '''
    def __init__(self):
        self.name    = None   # Name of model. Must be used in cases of branch heterogeneity, otherwise may remain None. When used, the name *MUST* correspond to its respective flag in the phylogeny.
        self.probs   = [1.]   # Rate heterogeneity site class probabilities. Default is single class with probability of 1.

        ##### MAYBE MOVE SHUFFLE HERE??? ######
        # self.shuffle_map = [ [], [], [] ...]
        #######################################
  
    def num_classes(self):
        ''' 
            How many rate classes?
        '''
        return len(self.probs)   
    
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
        Defines a Model() object. Used for nucleotide, amino-acid, and mutation-selection balance models.   
    '''
    
    def __init__(self):
        super(Model, self).__init__()
        self.params     = {}     # Parameters pertaining to substitution process. For all models, this includes a vector of stationary frequencies. Each individual evolutionary model will have its own additional parameters.
        self.matrix     = None   # Instantaneous rate matrix
        self.rates      = [1.]   # Rate heterogeneity model. List of rate factors. Default 1.0 (homogeneous). Must be same length as self.probs!


class CodonModel(EvoModels):
    '''
        Defines a CodonModel() object. 
        These objects are reserved for cases of *codon model heterogeneity* where dN/dS (or other) varies, and hence matrices vary.
    '''
    def __init__(self):
        super(CodonModel, self).__init__()
        self.matrices = []   # List of rate matrices associated with this codon model.
        
        
class Partition():
    '''
        Defines a Partition() object.
    '''
    def __init__(self):
        self.size           = []    # List of integers representing partition length. If there is no rate heterogeneity, then the list is length 1. Else, list is length k, where k is the number of rate categories.
        self.models         = None  # List of models associated with this partition. When length 1, temporally homogeneous.
        self.root_model     = None  # Model to begin at root of tree. Used under *branch heterogeneity*, and should be None or False if process is temporally homogeneous. If there is branch heterogeneity, this string *MUST* correspond to one of the Model() object's names and also a corresponding phylogeny flag.
        self.root_seq       = None  # User may choose to provide a root sequence for each partition, and it'll be stored here. Totally optional - will otherwise be generated from steady-state frequencies.
        self.shuffle        = False # Shuffle sites after evolving? 
    
    
    def branch_het(self):
        ''' 
            Return False is len(self.models) == 1 and True if len(self.models) > 1
        '''
        if isinstance(self.models, Model) or isinstance(self.models, CodonModel) or len(self.models) == 1:
            return False
        elif len(self.models) > 1:
            return True
        else:
            raise AssertionError("\n\nPartition has no associated models, so I have no clue what sort of heterogeneity there is...because there's no model...")
    
    
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
        
        
class Site():
    '''
        Defines a Site() object.
    '''
    def __init__(self):
        self.int_seq      = None # integer sequence at a site
        self.postition    = None # location of site in full alignment size. <- shuffle.
        self.rate         = None # This is a tuple containing (dN, dS) for codon models and for nuc/amino models it's the relative rate. Used for storing info.








             
                
        
        