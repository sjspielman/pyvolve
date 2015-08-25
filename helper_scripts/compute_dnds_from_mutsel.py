'''
    SJS.
    This script implements a module to compute dN/dS (or dN, dS individually) from mutation-selection model parameters, as described in Spielman and Wilke 2015 (doi: 10.1093/molbev/msv003).
    **Please cite the paper linked above if you use this script or calculations derived from it!**
    
    Usage:
        from compute_dnds_from_mutsel import *
    
        c = dNdS_from_MutSel(<codon_frequencies>, <mutation_dictionary>)
        # codon_frequencies: either a list, numpy array, or dictionary of equilibrium codon frequencies. If list or numpy array, frequencies should be in alphabetical order, excluding stop codons. The frequencies provided *must* sum to 1.
        # mutation_dictionary: **optional** argument indicating the mutation rates between nucleotides. This argument is analogous to the "mu" dictionary provided to pyvolve when simulating with custom mutation rates. If this argument is not provided, this module assumes equal mutation rates (e.g. a JC69 situation).
        
        c.compute_dnds() # to return dN/dS
        c.compute_ds()   # to return dS
        c.compute_dn()   # to return dN    
    
    
    
    Please post all questions, bugs, etc. to https://github.com/sjspielman/pyvolve/Issues
    
'''

import re
from copy import deepcopy
import numpy as np
from scipy import linalg
import sys
from random import uniform, shuffle



class dNdS_from_MutSel():
    '''
        Class to compute dN, dS, and/or dN/dS from mutation-selection model parameters, as described in Spielman and Wilke 2015.

        Positional arguments:
            codon_frequencies = either a list, numpy array, or dictionary of equilibrium codon frequencies. If list or numpy array, frequencies should be in alphabetical order, excluding stop codons. Relatively little sanity checking is done here, so provide something reasonable..
            mutation_dictionary = (optional) dictionary of mutation rates. This argument is analogous to the "mu" dictionary provided to pyvolve when simulating with custom mutation rates. If this argument is not provided, assumes equal mutation rates.
    '''
    def __init__(self, codon_frequencies, mutation_dictionary = None):
        
        self.ZERO = 1e-10
        self.codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}

        # Frequency setup
        if type(codon_frequencies) is dict:
            self.codon_freqs_dict = codon_frequencies
        else:
            self.codon_freqs_dict = {}
            for c in range(len(self.codons)):
                self.codon_freqs_dict[self.codons[c]] = codon_frequencies[c]
        assert(abs(1. - np.sum(self.codon_freqs_dict.values())) < self.ZERO), "\n\nProvided codon frequencies must sum to 1."
        
        
        # Mutation setup        
        if mutation_dictionary is None:
            self.mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
        else:
            self.mu_dict = deepcopy(mutation_dictionary)
            for key in mutation_dictionary: 
                if key[1] + key[0] not in self.mu_dict:
                    self.mu_dict[key[1] + key[0]] = mutation_dictionary[key]



    def compute_dn(self):
        ''' 
            Compute dN from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        return self.compute_quantity("nonsyn")




    def compute_ds(self):
        ''' 
            Compute dS from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        return self.compute_quantity("syn")




    def compute_dnds(self):
        ''' 
            Compute dN/dS from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        dn = self.compute_dn()
        ds = self.compute_ds()
        
        if dn <= self.ZERO:
            return 0.
        elif ds <= self.ZERO:
            return np.inf
        else:
            return dn/ds


    def compute_quantity(self, type):
        '''
            Compute either dN or dS (type is 'nonsyn' or 'syn', respectively).
        '''
        numer = 0.
        denom = 0.

        for codon in self.codon_freqs_dict:
            rate = 0.
            sites = 0.
            rate, sites = self.calc_paths(codon, type)
            numer += rate
            denom += sites
        
        assert( denom > self.ZERO ), "\n\nProvided frequencies indicate no evolution is 'possible'."
        return numer/denom




     
    def calc_paths(self, source, type):
        ''' 
            Compute a term in the quantity numerator for a given source codon.
        '''
        rate = 0.
        sites = 0.
        source_freq = self.codon_freqs_dict[source]
        for target in self.codons:
            diff = self.get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
            if (type == 'nonsyn' and self.codon_dict[source] != self.codon_dict[target]) or (type == 'syn' and self.codon_dict[source] == self.codon_dict[target]):
                if len(diff) == 2:
                    rate  += self.calc_subst_prob( source_freq, self.codon_freqs_dict[target], self.mu_dict[diff], self.mu_dict[diff[1]+diff[0]] )
                    sites += self.mu_dict[diff]
        rate  *= source_freq
        sites *= source_freq
        return rate, sites



    def get_nuc_diff(self, source, target):
        '''
            Determine nucleotide difference between two codons.
        '''
        diff = ''
        for i in range(3):
            if source[i] != target[i]: 
                diff += source[i]+target[i]
        return diff

 
    
    def calc_subst_prob(self, pi, pj, mu_ij, mu_ji):
        '''
            Compute the substitution probability between two codons, as in Halpern and Bruno 1998.
        '''
        if abs(pi) <= self.ZERO or abs(pj) <= self.ZERO:
            fixation_rate = 0.
        else:
            p_mu = (mu_ji*pj)/(mu_ij*pi)
        
            # If p_mu == 1, L'Hopitals gives fixation rate of 1 (substitution probability is the forward mutation rate) 
            if abs(1. - p_mu) <= self.ZERO:
                fixation_rate = 1. 
            else:
                fixation_rate =  (np.log(p_mu))/(1. - (1./p_mu))
        return fixation_rate * mu_ij            










