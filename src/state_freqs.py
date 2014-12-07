#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Compute a vector of stationary frequencies.
'''


import os
import re
import time
import numpy as np
import random as rn
from Bio import SeqIO, AlignIO
from genetics import *
ZERO      = 1e-8
MOLECULES = Genetics()



class StateFrequencies(object):
    '''
    Parent class for stationary frequency calculations.
    
    Overall strategy for frequency calculations:
        User determines the alphabet in which computations should occur (by = <amino/codon/nuc>). In this way, stationary frequencies can be calculated with, say, codons, and then amino acid frequencies can ultimately be returned.
        
        
    Child classes include the following:
    1. *EqualFrequencies*          : Sets frequencies as equal (i.e. 1/4 for all nucleotides if by='nuc', and so on.) *DEFAULT FREQUENCY CALCULATIONS*
    2. *RandomFrequencies*         : Computes (semi-)random frequency values for a given alphabet.
    3. *CustomFrequencies*         : Computes frequencies based on a user-provided dictionary of frequencies.
    4. *ReadFrequencies*           : Computes frequencies from a sequence file. Contains an option to select specific columns from sequence file only, but this requires that the file is an alignemnt.
    
    REQUIRED ARGUMENTS:
    1. by = nuc/amino/codon. Specifies the alphabet used to compute frequencies.
    
    
    
    OPTIONAL ARGUMENTS:
    1. *type* represents the final frequencies to return. If not specified, then the 'by' frequencies will be returned. If specified, the 'by' frequencies will be converted to the 'type' alphabet.
    
    2. *restrict* may can be used in conjuction with either the EqualFrequencies or RandomFrequencies subclasses (else, ignored). 
        This argument should be a list (each element a string) which should have non-zero stationary frequencies. 
        For instance, if (by = 'amino', restrict = ['C', D', 'G']) is specified, then frequencies will be computed in amino-acid space, and all amino acids EXCEPT cysteine, aspartate, and glycine will automatically be given a frequency of 0.
        If this argument is not provided, all states are allowed to have non-zero frequencies.
        
    3. *codon_bias* may be used in conjuction with any child class, and is meant specifically (only!) for applying codon bias when by=amino is specified by final codon frequencies are desired [i.e. (by = 'amino', type = 'codon')]. This conversion will take place by randomly selected a preferred codon for each amino acid, and assigning it a higher frequency.
        This argument should be a float (decimal) value between (0,1]. The overall amino acid frequency is preserved, but partitioned among constituent codons such that a single codon is preferred and rest non-preferred.
        For instance, if codon_bias = 0.5 and a given 4-fold degenerate amino acid's overall frequency is 0.1, then a (*randomly selected*) preferred codon for this amino acid will have a frequency of 0.05, and the remaining three non-preferred codons will each have frequencues of 0.05/3.

    '''
    
    
    def __init__(self, **kwargs):
        
        # Frequency vectors "initialized". It is possible that not all of these will be used, but we set them up in case. 
        self.nuc_freqs    = np.zeros(4)     
        self.amino_freqs  = np.zeros(20)
        self.codon_freqs  = np.zeros(61)
        
        # Input parameters and general setup. 
        self._by = kwargs.get('by')
        assert(self._by =='amino' or self._by == 'codon' or self._by == 'nuc'), "\n\nYou have either no 'by' or a wrong 'by'. Remember, codon, amino, or nuc only!"
        self._set_code_size()

        self._byFreqs     = np.zeros(self._size)        
        self._restrict    = kwargs.get('restrict', self._code)
        self._codon_bias  = kwargs.get('codon_bias', None)    # To implement codon bias, can provide a decimal giving the percent usage of the preferred state. NOTE: CURRENTLY THE PREFERRED STATE IS RANDOMLY CHOSEN.
        self._savefile    = kwargs.get('savefile', None)      # for saving the equilibrium frequencies to a file

        if self._codon_bias:
            assert(ZERO < self._codon_bias <= 1.0), "*codon_bias* argument must be a float between (0,1]."
        if self._restrict is not self._code:
            assert(type(self._restrict) is list), "*restrict* must be a list of state strings corresponding to the 'by' argument. For instance, you may use (by = 'amino', restrict = ['A', 'C', 'G', 'P'])."
        
        
        
    def __call__(self, **kwargs):
        ''' 
        
            Calculate and return state frequencies. At this stage, the StateFrequency object must already have been initialized with a by = <amino/codon/nuc>.
            Will return a vector of stationary frequencies, based on the argument *type*.
            
            OPTIONAL ARGUMENTS:
                1. *type* is the type of frequencies to return. For instance, if by = amino was previously specified, the user can call calculate_freqs(type = 'codon') to obtain codon frequencies from those amino acid frequencies already computed.
                    If *type* is not provided, then calculate_freqs will return the same frequencies as "by".
                
                2. *savefile* is a file to which final frequencies will be saved. Output frequencies will be ordered alphabetically, i.e. for amino acids, the order will be A, C, D, etc. and for codons AAA, AAC, AAG, etc.                                 
        '''
        
        # Input arguments and general setup
        type = kwargs.get('type', self._by)
        assert(type =='amino' or type == 'codon' or type == 'nuc'), "Can only calculate codon, amino, or nuc frequencies."
        if type == 'amino' or type == 'codon':
            assert(self._by == 'amino' or self._by == 'codon'), "\n\nIncompatible *type* argument! If you would like to obtain amino acid or codon frequencies, your 'by' argument MUST be either codon or ami calculations must use either amino acids or codons, NOT nucleotides."
        savefile = kwargs.get('savefile', None)
        if self._codon_bias is not None:
            assert(self._by == 'amino' and type == 'codon')

        # Create the self._byFreqs, if does not already exist. Once created, assign as either amino, codon, nuc frequencies.
        if np.array_equal(self._byFreqs, np.zeros(self._size)):
            self._generate_byFreqs()  
            assert( abs(np.sum(self._byFreqs) - 1.) < ZERO), "State frequencies improperly generated. Do not sum to 1." 
            self._assign_byFreqs()
        
        # Convert frequencies if needed
        if type != self._by:
            conv_expr = "self._"+self._by+"_to_"+type+"()"
            eval(conv_expr)
        
        # Save if needed
        if savefile is not None:
            np.savetxt(savefile, eval("self."+type+"_freqs"), fmt='%.5e')
        return eval("self."+type+"_freqs")        


    def _set_code_size(self):
        ''' Set up the code (alphabet) and dimensionality for computing self._byFreqs '''
        if self._by == 'amino':
            self._code = MOLECULES.amino_acids
        elif self._by == 'codon':
            self._code = MOLECULES.codons
        elif self._by == 'nuc':
            self._code = MOLECULES.nucleotides
        self._size = len(self._code)
 
 

   
    
    
    ############################################# FREQUENCY CONVERSIONS ###############################################
    def _amino_to_codon(self):
        ''' 
            Calculate codon frequencies from amino acid frequencies. (by = 'amino', type = 'codon')
            Unless codon bias is specified, assumes equal frequencies for synonymous codons.
        '''
        
        for aa_count in range(20):
            syn = MOLECULES.genetic_code[aa_count]
            if self._codon_bias:
                self._apply_codon_bias(aa_count, syn)
            else:
                for synCodon in syn:
                    cind = MOLECULES.codons.index(synCodon)
                    self.codon_freqs[cind] = self.amino_freqs[aa_count]/float(len(syn))
        assert( abs(np.sum(self.codon_freqs) - 1.) < ZERO), "Codon state frequencies improperly calculated from amino acid frequencies. Do not sum to 1."                 
      
    
    
    def _apply_codon_bias(self, aa_count, syn):
        ''' 
            Implements codon bias for a given amino acid. (by = 'amino', type = 'codon', codon_bias = (0,1]).
            Note that if we are dealing with either M or W (non-degenerate), we'll just break out.
            Arguments:
                1. *aa_count* is the amino acid index we are working with
                2. *syn* the list of synonymous codons for this amino acid index        
        '''

        sum = 0.
        
        # M and W are non-degenerate, hence bias is not applicable.
        if len(syn) == 1:
            cind = MOLECULES.codons.index(syn[0])    
            self.codon_freqs[cind] = self.amino_freqs[aa_count]
            sum += self.amino_freqs[aa_count]
        
        # Degenerate amino acids
        else:
            pref_index    = rn.randint(0, len(syn)-1)
            pref_freq     = self.amino_freqs[aa_count] * self._codon_bias
            nonpref_freq  = (self.amino_freqs[aa_count] - pref_freq)/(len(syn) - 1.)  
            
            for s in range(len(syn)):
                cind = MOLECULES.codons.index(syn[s])
                if s == pref_index:
                    self.codon_freqs[cind] = pref_freq
                    sum += pref_req
                else:
                    self.codon_freqs[cind] = nonpref_freq
                    sum += nonpref_freq
                    
        assert( abs(sum - self.amino_freqs[aa_count]) < ZERO ), "Codon bias improperly implemented."
            
    def _codon_to_amino(self):
        ''' 
            Calculate amino acid frequencies from codon frequencies (by = 'codon', type = 'amino').
        '''
        
        for a in range(len(MOLECULES.amino_acids)):
            codons1 = MOLECULES.genetic_code[a]
            for c in codons1:
                ind = MOLECULES.codons.index(c)
                self.amino_freqs[a] += self.codon_freqs[ind]
        assert( abs(np.sum(self.amino_freqs) - 1.) < ZERO), "Amino acid state frequencies improperly generate_byFreqsd from codon frequencies. Do not sum to 1." 


    def _codon_to_nuc(self):
        ''' 
            Calculate nucleotide frequencies from codon frequencies (by = 'codon', type = 'nuc').
        '''        
        
        for i in range(61):
            codon_freq = self.codon_freqs[i]
            codon = MOLECULES.codons[i]
            for n in range(4):
                nuc =  MOLECULES.nucleotides[n]
                nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
                self.nuc_freqs[n] += codon_freq * nuc_freq
        assert( abs(np.sum(self.nuc_freqs) - 1.) < ZERO), "Nucleotide state frequencies improperly generate_byFreqsd. Do not sum to 1." 

        
    def _amino_to_nuc(self):
        ''' 
            Calculate nucleotide frequencies from amino acid frequencies (by = 'amino', type = 'nuc').
        '''

        self._amino_to_codon()
        self._codon_to_nuc()  
     
    
    #####################################################################################   

    def _assign_byFreqs(self):
        ''' 
            Called from within function calcuate_freqs, this function will assign a frequency vector to the appropriate attribute variable.
        '''
        if self._by == 'codon':
            self.codon_freqs = self._byFreqs
        elif self._by == 'amino':
            self.amino_freqs = self._byFreqs
        elif self._by == 'nuc':
            self.nuc_freqs = self._byFreqs
        else:
            raise AssertionError("WHAT ARE WE DOING HERE.")









class EqualFrequencies(StateFrequencies):
    ''' 
        DEFAULT child class of StateFrequencies. Computes equal state frequencies (amino = 1/20, codon = 1/61, nucleotide = 1/4).
        May be used in conjuction with self._restrict attribute.

    '''
    
    def __init__(self, **kwargs):
        super(EqualFrequencies, self).__init__(**kwargs)
    
    def _generate_byFreqs(self):
        '''
            Compute self._byFreqs
        '''
        fill = 1./float(len(self._restrict))
        for entry in self._restrict:
            self._byFreqs[self._code.index(entry)] = fill     
                    
                    
                    
                    
                    
                    
                    
class RandomFrequencies(StateFrequencies):
    ''' 
        Child class of StateFrequencies. Computes "random" state frequencies. 
        The distributions are not, in fact, completely random, but will return virtually flat distributions with some small amount of noise.
        May be used in conjuction with self._restrict attribute.
        
    '''
    def __init__(self, **kwargs):
        super(RandomFrequencies, self).__init__(**kwargs)
        self._partial_restrict = self._restrict[:-1] # all but last
        
      
    def _generate_byFreqs(self):
        '''
            Compute self._byFreqs. Since random sampling, we can run into timing issues. Make sure we don't get stuck!!
        '''
        max = 2./len(self._restrict)
        min = 1e-5
        abort_after_time = 0.001 + time.time()
        
        restart_search = True
        while restart_search:
            restart_search = False
            sum = 0.
            self._byFreqs = np.zeros(self._size)
            for entry in self._partial_restrict:
                freq = rn.uniform(min,max)
                while (sum + freq > 1):
                    freq = rn.uniform(min,max)
                    if time.time() > abort_after_time:
                        restart_search = True 
                        break
                if restart_search:
                    break
                sum += freq
                self._byFreqs[self._code.index(entry)] = freq
        self._byFreqs[self._code.index(self._restrict[-1])] = (1.-sum)    
        




class CustomFrequencies(StateFrequencies):
    ''' 
        Child class of StateFrequencies. Computes frequencies using directly with a user-input python dictionary of frequencies.
        Note that 'by' should correspond to the sort of frequencies that they've entered. 'type' should correspond to what they want at the end.
        For instance, it is possible to provide amino acid frequencies and ultimately obtain codon frequencies (with synonymous treated equally, in this circumstance).
        
        REQUIRED ARGUMENTS:
            1. *freq_dict* is a python dictionary of frequencies, in which keys are states (e.g. a codon key would be 'ATC', an amino acid key would be 'W', and a nucleotide key would be 'T'), and values are float frequencies.
                Any states not included in this dictionary are assumed to have an equal frequency. Hence, the dictionary values *MUST* sum to 1.
                This dictionary's keys (alphabet) must correspond to the alphabet specified in 'by'.
                Note that at this time we support only single-letter amino acid (NOT 3-letter!) codes.  
    
    '''
    def __init__(self, **kwargs):
        super(CustomFrequencies, self).__init__(**kwargs)    
        self.given_freqs = kwargs.get('freq_dict', {}) # Dictionary of desired frequencies.    
        self._sanity_freq_dict()                        # Quick sanity check on frequencies


    def _sanity_freq_dict(self):
        ''' 
            Sanity check to ensure the following:
                1. self._by is the same alphabet as the freq_dict keys
                2. freq_dict keys are consistent.
                3. frequencies sum to 1
        '''       
        prob_sum = 0.
        for entry in self.given_freqs:
            assert( len(entry) == len(self._code[0]) and entry in self._code ), "\n\n Your *freq_dict* keys are not properly format. Please ensure that your keys correspond to the *by* calculations, and that you only specify canonical amino acids/nucleotide, or  sense codons."  
            prob_sum += float(self.given_freqs[entry])
        assert( abs( 1. - prob_sum) < ZERO), "\n\nFrequencies provided in *freq_dict* do not sum to 1!"
      
 
    def _generate_byFreqs(self):
        '''
            Compute self._byFreqs.
        '''
        for i in range(self._size):
            element = self._code[i]
            if element in self.given_freqs:
                self._byFreqs[i] = self.given_freqs[element]







class ReadFrequencies(StateFrequencies):
    ''' 
        Child class of StateFrequencies. Computes frequencies from a given sequence file. 
        Frequencies may be computed globally (using entire file), or based on specific columns (i.e. site-specific frequencies) in the file, provided the file contains a sequence alignment.

        
        REQUIRED ARGUMENTS:
            1. *file* is the file containing sequences from which we will obtain state frequencies. 
                This file is assumed to be in FASTA format, although you can specify a different format with the argument *format*
            
        OPTIONAL ARGUMENTS:
            1. *format* is the sequence file format. We parse sequence files using Biopython, so we accept the following formats only: fasta, phylip, phylip-relaxed, nexus.
        
            2. *columns* is a list of integers giving the column(s) which should be considered in frequency calculations.
                This list should be indexed *from 1*.
                If this argument is not provided, all positions in sequence file will be considered. 
    
     ''' 
     
    def __init__(self, **kwargs):
        super(ReadFrequencies, self).__init__(**kwargs)
        
        # Input variables, options
        self.seqfile          = kwargs.get('file', None)
        self.format           = kwargs.get('format', 'fasta').lower()   # Biopython requires that this flag is lowercase.
        self.which_columns    = kwargs.get('columns', None)
        self._seqs            = []                                      # Sequence records obtained from sequence file
        self._make_seq_list()

    
    def _sanity_which_columns(self):
        ''' 
            Sanity check that *columns* argument has been properly specified.
                1. Should be a list
                2. Should not include columns outside of the length of the alignment [1,alnlen]
                3. Should be converted to a numpy array
        '''
        try:    
            AlignIO.read(self.seqfile, self.format)
        except:
            raise AssertionError("\n\nYour sequence file does not appear to be an *alignment.* If you would like to get frequencies from specific columns only, it must be an alignment!") 
        assert( type(self.which_columns) is list), "\n\nArgument *which_columns* should be a list of integers giving the column(s) (indexed from 1!) which should be considered for frequency calculations."
        self.which_columns = np.array(self.which_columns) - 1
        assert( np.all(self.which_columns) >= 0 and np.all(self.which_columns) <= self._alnlen - 1), "\n\nYour column indices specified in *which_columns* do not play well with alignment! Remember that column indexing starts at *1*, and you cannot specify columns that don't exist."
        
        
        
        
    def _make_seq_list(self):
        ''' 
            Read in sequence file and set up variables used in frequency calculations.
            Additionally performs some sanity checks on sequence file, sequences themselves, and the which_columns argument (if specified).
         '''
         
        assert(self.seqfile is not None), "\n\n You must provide a sequence/alignment file with the argument file=<my_file_name> to use the ReadFrequencies class."
        assert(os.path.exists(self.seqfile)), "\n\n Your input file does not exist! Check the path?"
        try:
            raw = list(SeqIO.parse(self.seqfile, self.format))
        except:
            raise AssertionError("\n\nYour sequence file could not be parsed. Note that if your sequence file is not in FASTA format, you must specify its format with the argument *format*.")  
        self._numseq = len(raw)
        self._alnlen = len(raw[0]) # This will only come into play if we're collecting columns.
        if self._by == 'codon':
             assert( self._alnlen%3 == 0), "\n\nThe length of your sequence alignment is not a multiple of three, so you don't seem to actually have codons."
        if self.which_columns is not None:
            self._sanity_which_columns()
        for entry in raw:
            self._seqs.append(str(entry.seq))  
    
    
    def _generate_byFreqs(self):
        ''' 
            Compute self._byFreqs
        ''' 
        
        total_characters = 0.
        for row in self._seqs: 
            if self.which_columns is not None:   
                for col in self.which_columns:
                    if self._by == "codon": 
                        char = row[col*3 : col*3 + 3]
                    else:
                        char = row[col]
                    if char in self._code:
                        total_characters += 1.
                        self._byFreqs[ self._code.index(char) ] += 1          
            else:
                for i in range(len(row)):
                    if self._by == "codon":
                        char = row[i*3 : i*3+3]
                    else:
                        char = row[i]
                    if char in self._code:
                        total_characters += 1.
                        self._byFreqs[ self._code.index(char) ] += 1     
        self._byFreqs = np.divide(self._byFreqs, total_characters)





class EmpiricalModelFrequencies():
    ''' 
        This class is used only to return default frequencies (i.e. those from original papers) of empirical models (all amino acid models and the ECM models). 
        Note that this is *not* a child class StateFrequencies, because it uses none of its functionality whatsoever, although it is used in the same way that any StateFrequencies class is used.
        Empirical model matrices and corresponding frequency vectors are stored in src/empirical_matrices.py . 
        Currently supported models:
            1. Amino acid: JTT, WAG, LG
            2. Codon:      ECM(un)rest
            [NOTE: We additionally support the empirical codon model SCG05, but this model does not have its own stationary frequencies.]
        
        REQUIRED ARGUMENTS:
            1. *model* is empirical model of choice. This argument should be specified as any of the following: JTT, WAG, LG, ECMrest, ECMunrest. Argument is case insensitive, so have at it.
     
     ''' 
    
    def __init__(self, **kwargs):
        try:
            self.empirical_model = kwargs.get('model', None).lower()
        except KeyError:
            print "\n\n You must specify an empirical model with *model* to obtain its frequencies."
        

    def __call__(self):    
        ''' 
            Simply return the default empirical frequencies.
        '''
        import empirical_matrices as em
        try:
            return np.array( eval("em."+self.empirical_model+"_freqs") )
        except:
            print "Couldn't figure out your empirical model specification! We only support the following empirical models (for frequency specification):"
            print "Amino acid: JTT, WAG, LG."
            print "Codon:      ECM restricted or unrestricted, which can be specified respectively as ECMrest and ECMunrest (case insensitive)."
            print "I'm quiting :/"
            sys.exit()




