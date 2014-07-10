#### TO DO IDEAS:
###### Set the GC content, for type codon or nucleotide frequencies ONLY.


import os
import re
import numpy as np
import random as rn
from Bio import SeqIO
from misc import ZERO, Genetics
MOLECULES = Genetics()




class StateFreqs(object):
    '''Will return frequencies. huzzah.'''
    def __init__(self, **kwargs):
        
        # Generics
        self.nuc_freqs    = np.zeros(4)     
        self.amino_freqs  = np.zeros(20)
        self.codon_freqs  = np.zeros(61)
        
        # Input options and some sanity checking
        self.by = kwargs.get('by') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. If nucleotide, well, yeah.
        assert(self.by =='amino' or self.by == 'codon' or self.by == 'nuc'), "\n\nYou have either no 'by' or a wrong 'by'. Remember, codon, amino, or nuc only!"
        self._set_code_size()
        self.byFreqs = np.zeros(self._size) # keep this naming convention since it's a weird variable.
        
        self.restrict   = kwargs.get('restrict', self._code) # For the equal, rand subclasses only.
        self.constraint = kwargs.get('constraint', 1.0) # For the user, read subclasses only. Constrain provided amino acids to be a certain percentage of total equilbrium frequencies. This allows for non-zero propensities throughout, but non-preferred will be exceptionally rare.
        self.bias       = kwargs.get('codonBias', None) # To implement codon bias, can provide a decimal giving the percent usage of the preferred state. NOTE: CURRENTLY THE PREFERRED STATE IS RANDOMLY CHOSEN.
        self.savefile   = kwargs.get('savefile', None) # for saving the equilibrium frequencies to a file

        if self.bias:
            assert(ZERO < self.bias <= 1.0), "Codon bias must be >0, <=1."
        if self.constraint:
            assert(ZERO <  self.constraint <= 1.0), "Constraint must be >0, <=1."
        if self.restrict is not self._code:
            assert(type(self.restrict) is list), "Restriction must be a list."
        




    def _set_code_size(self):
        ''' Set the codes and lengths ''' 
        if self.by == 'amino':
            self._code = MOLECULES.amino_acids
        elif self.by == 'codon':
            self._code = MOLECULES.codons
        elif self.by == 'nuc':
            self._code = MOLECULES.nucleotides
        self._size = len(self._code)

           
       
    def _unconstrain_frequencies(self):
        ''' This function will allow for some frequency constraints to be lessened for the self.byFreqs
            FUNCTION MAY BE USED BY USERFREQS AND READFREQS ONLY.
            If the constraint value is 0.95, then the preferred (non-zero frequency) entries should only sum to 0.95.
            The remaining 0.05 will be partitioned equally among the non-preferred (freq = 0) entries.
            Therefore, this function allows for some evolutionary "wiggle room" while still enforcing a strong preference.
        '''
        self.byFreqs = np.multiply(self.byFreqs, self.constraint)
        assert (self._size > np.count_nonzero(self.byFreqs)), "All state frequencies are 0! This is problematic for a wide variety of reasons."
        addToZero = float( (1.0 - self.constraint) / (self._size - np.count_nonzero(self.byFreqs)) )
        for i in range(self._size):
            if ( abs(self.byFreqs[i] - 0.0) < ZERO):
                self.byFreqs[i] = addToZero
        assert( abs( np.sum(self.byFreqs) - 1.0) < ZERO), "unconstraining frequencies did not work properly - freqs don't sum to 1."
        
        
    def _apply_codon_bias(self, aa_count, syn):
        ''' Implements codon bias. There is a self.bias param which gives a decimal indicating the frequency of the preferred codon.
            Requires by=amino, type=codon
            TO DO, FUTURE DEVELOPMENT: Allow users to provide the preferred codons. As of now, the preferred codon is RANDOM.
            Args: aa_count = the amino acid index we are working with
                  syn      = the list of synonymous codons for this amino acid index        
        '''
        # If we are dealing with a single-codon amino acid (M or W), simply fill the value since codon bias is not possible.
        sum = 0. # for some debugging
        if len(syn) == 1:
            cind = MOLECULES.codons.index(syn[0])    
            self.codon_freqs[cind] = self.amino_freqs[aa_count]
            sum += self.amino_freqs[aa_count]
        else:
            prefIndex = rn.randint(0, len(syn)-1)
            prefFreq = self.amino_freqs[aa_count] * self.bias
            nonprefFreq = (self.amino_freqs[aa_count] - prefFreq)/(len(syn) - 1.)  
            
            for s in range(len(syn)):
                cind = MOLECULES.codons.index(syn[s])
                if s == prefIndex:
                    self.codon_freqs[cind] = prefFreq
                    sum += prefFreq
                else:
                    self.codon_freqs[cind] = nonprefFreq
                    sum += nonprefFreq
        assert(abs (sum - self.amino_freqs[aa_count]) < ZERO), "Codon bias improperly implemented."

   
    
    
    ######################################### FREQUENCY CONVERSIONS #########################################
    def _amino_to_codon(self):
        ''' Calculate codon frequencies from amino acid frequencies. 
            Unless codon bias is specified, will assume equal synonymous frequencies.
        '''
        for aa_count in range(20):
            syn = MOLECULES.genetic_code[aa_count]
            if self.bias:
                self._apply_codon_bias(aa_count, syn)
            else:
                for synCodon in syn:
                    cind = MOLECULES.codons.index(synCodon)
                    self.codon_freqs[cind] = self.amino_freqs[aa_count]/float(len(syn))
        assert( abs(np.sum(self.codon_freqs) - 1.) < ZERO), "Codon state frequencies improperly generate_byFreqsd from amino acid frequencies. Do not sum to 1."                 
                
                
    def _codon_to_amino(self):
        ''' Calculate amino acid frequencies from codon frequencies. ''' 
        for a in range(len(MOLECULES.amino_acids)):
            codons1 = MOLECULES.genetic_code[a]
            for c in codons1:
                ind = MOLECULES.codons.index(c)
                self.amino_freqs[a] += self.codon_freqs[ind]
        assert( abs(np.sum(self.amino_freqs) - 1.) < ZERO), "Amino acid state frequencies improperly generate_byFreqsd from codon frequencies. Do not sum to 1." 
    
    def _codon_to_nuc(self):
        ''' Calculate nucleotide frequencies from codon frequencies. '''
        for i in range(61):
            codon_freq = self.codon_freqs[i]
            codon = MOLECULES.codons[i]
            for n in range(4):
                nuc =  MOLECULES.nucleotides[n]
                nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
                if nuc_freq > 0 :
                    self.nuc_freqs[n] += codon_freq * nuc_freq
        assert( abs(np.sum(self.nuc_freqs) - 1.) < ZERO), "Nucleotide state frequencies improperly generate_byFreqsd. Do not sum to 1." 

        
    def _amino_to_nuc(self):
        ''' Calculate nucleotide frequencies from amino acid frequencies. Lazy function, hurray!'''
        self._amino_to_codon()
        self._codon_to_nuc()  
    
    #####################################################################################   

    def _assign_byFreqs(self):
        ''' Assign self.byFreqs to either amino, codon, or nuc. '''
        if self.by == 'codon':
            self.codon_freqs = self.byFreqs
        elif self.by == 'amino':
            self.amino_freqs = self.byFreqs
        elif self.by == 'nuc':
            self.nuc_freqs = self.byFreqs
        else:
            raise AssertionError("WHAT ARE WE DOING HERE.")


    def calculate_freqs(self, **kwargs):
        ''' Calculate and return state frequencies.            
            State frequencies are calculated for whatever "by" specifies.
            The "type", as provided here, will be what we return to users.
            
        '''
        type = kwargs.get('type', self.by)
        assert(type =='amino' or type == 'codon' or type == 'nuc'), "Can only calculate codon, amino, or nuc frequencies."
        if type == 'amino' or type == 'codon':
            assert(self.by == 'amino' or self.by == 'codon'), "\n\nIncompatible by! For amino acid or codon frequencies, calculations must use either amino acids or codons, NOT nucleotides."
        save = kwargs.get('savefile', None)
        if self.bias is not None:
            assert(self.by == 'amino' and type == 'codon')

        # Create the self.byFreqs, if does not already exist. Once created, assign as either amino, codon, nuc frequencies.
        if np.array_equal(self.byFreqs, np.zeros(self._size)):
            self._generate_byFreqs() # generate_byFreqss self.byFreqs       
            assert( abs(np.sum(self.byFreqs) - 1.) < ZERO), "State frequencies improperly generate_byFreqsd. Do not sum to 1." 
            self._assign_byFreqs()
        
        # Convert frequencies if needed
        if type != self.by:
            conv_expr = "self._"+self.by+"_to_"+type+"()"
            eval(conv_expr)
        
        # Save if needed
        if save is not None:
            np.savetxt(save, eval("self."+type+"_freqs"), fmt='%.5e')
        return eval("self."+type+"_freqs")




class BoltzmannFreqs(StateFreqs):
    ''' Return state frequencies generate_byFreqsd by Boltzmann distribution.
        !!! This sub-class works only for by='amino !!!
        Default factor is 1.0. Increase for more constraint, decrease for less constraint.
        This option requires a ranking (list of aminos, must be complete). 
        If no ranking is provided, aminos randomly ranked.
    '''
    def __init__(self, **kwargs):
        super(BoltzmannFreqs, self).__init__(**kwargs)
        
        assert(self.by == 'amino'), "Boltzmann-type frequencies may only be done with amino acid calculations. Try again."    
        self.factor  = float( kwargs.get('factor', 1.0) )
        self.ranking = kwargs.get('rank', None)
        if self.ranking is None:
            self.ranking = self._code
        else:
            assert( type(self.ranking) is list ), "Ranking must be a full list of amino acids in order."
            assert( sorted(self.ranking) == self._code ), "Your ranking list does not appear to contain all amino acids, or has incorrect letters in it."
    
    def _set_Boltzmann(self):
        ddg_values = np.random.normal(size = self._size) 
        numer_list = np.zeros(self._size)
        denom = 0.
        for d in range(self._size):
            val = np.exp(-1. * self.factor * ddg_values[d])
            denom += val
            numer_list[d] = val
        return numer_list/(np.sum(numer_list)) 
    
    def _generate_byFreqs(self):
        tempFreqs = self._set_Boltzmann()
        if self.ranking is not self._code:
            tempFreqs = np.sort(tempFreqs)[::-1]
        count = 0
        for aa in self.ranking:
            self.byFreqs[self._code.index(aa)] = tempFreqs[count]
            count += 1
                





class EqualFreqs(StateFreqs):
    ''' Return equal state frequencies. 
        NOTE: THIS IS THE DEFAULT BEHAVIOR.
    '''
    
    def __init__(self, **kwargs):
        super(EqualFreqs, self).__init__(**kwargs)
    
    def _generate_byFreqs(self):
        fillValue = 1./float(len(self.restrict))
        for entry in self.restrict:
            self.byFreqs[self._code.index(entry)] = fillValue                
                    
                    
                    
                    
                    
class RandFreqs(StateFreqs):
    ''' Return random state frequencies.
        Will return essentially flat distributions, but with noise.
    '''
    def __init__(self, **kwargs):
        super(RandFreqs, self).__init__(**kwargs)
      
    def _generate_byFreqs(self):
        partial_restrict = self.restrict[:-1] # all but last
        max = 2./len(self.restrict)
        min = 1e-5
        sum = 0.
        for entry in partial_restrict:
            freq = rn.uniform(min,max)
            while (sum + freq > 1):
                freq = rn.uniform(min,max)
            sum += freq
            self.byFreqs[self._code.index(entry)] = freq
        self.byFreqs[self._code.index(self.restrict[-1])] = (1.-sum)    
    
    
    




class UserFreqs(StateFreqs):
    ''' Assign frequencies based on user input. Assume that if not specified, the frequency is zero. 
        Note that 'by' should correspond to the sort of frequencies that they've entered. 'type' should correspond to what they want at the end.
        For instance, it is possible to provide amino acid frequencies and ultimately obtain codon frequencies (with synonymous treated equally, in this circumstance).
        
        NOTE: UNCONSTRAINING IS POSSIBLE HERE.
    
    '''
    def __init__(self, **kwargs):
        super(UserFreqs, self).__init__(**kwargs)    
        self.givenFreqs = kwargs.get('freqs', {}) # Dictionary of desired frequencies.    
        self._check_by_keys() ######## this will likely be removed when formal sanity checking is implemented eventually ######


    def _check_by_keys(self):
        ''' To make sure that self.by is the same alphabet as provided in the dictionary and that keys are ok.'''
        keysize = len( str(self.givenFreqs.keys()[0]) ) # Size of first key. All other keys should be the same size as this one. NOTE THAT IF THIS IS REALLY NOT A STRING, IT WILL BE CAUGHT LATER!! Perhaps/definitely this is inelegant, but I'll deal w/ it later.
        for key in self.givenFreqs.keys():
            assert( len(key) == keysize), "\n\n BOOO keys not all same size."
        if keysize == 3:
            assert(self.by == 'codon'), "bad keys"
        elif keysize == 1:
            assert (self.by == 'nuc' or self.by == 'amino'), "bad keys"
        else:
            raise AssertionError("\n\nBad dictionary keys for userfreqs.")
    
    def _generate_byFreqs(self):
        for i in range(self._size):
            element = self._code[i]
            if element in self.givenFreqs:
                self.byFreqs[i] = self.givenFreqs[element]
        if self.constraint < 1.0:
            self._unconstrain_frequencies()


class ReadFreqs(StateFreqs):
    ''' Retrieve frequencies from a file. Can either do global or specify a particular column/group of columns.
        NOTE: UNCONSTRAINING IS POSSIBLE HERE.
        
        TO DO: SANITY CHECKING WILL NEED TO VERIFY THAT THE PROVIDED SEQUENCE FILE IS IN THE SAME ALPHABET AS THE BY IS.
     ''' 
    def __init__(self, **kwargs):
        super(ReadFreqs, self).__init__(**kwargs)
        
        self.seqfile  = kwargs.get('file', None)   # Can also read frequencies from a sequence file
        self.format   = kwargs.get('format', 'fasta') # Default for that file is fasta
        self.which_columns = kwargs.get('columns', None)     # Which columns we are collecting frequencies from. Default is all columns combined. IF YOU GIVE IT A NUMBER, INDEX AT 0!!!!
        self.seqs     = [] # Sequence records obtained from sequence file
        self._full_sequence  = '' # Single sequence string from which to obtain frequencies
        self.DNA_characters  = re.compile(r"[^ACGT]") # DNA regexp for what to keep
        self.AMINO_characters = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]") # protein regexp for what to keep
       
        
    def _make_seq_list(self):
        ''' Set up sequences and relevent variables for frequency collection. '''
        raw = list(SeqIO.parse(self.seqfile, self.format))
        self.seqs = []
        self.numseq = len(raw)
        self.alnlen = len(raw[0]) # This will only come into play if we're collecting columns.
        for entry in raw:
            self.seqs.append(str(entry.seq))  

    
    def _process_seq_list(self):
        ''' 
            Turns full sequence we want to grab frequencies from into a single string.
            NOTE: If we want columns, we must get a string of the specific columns we're collecting from.
        ''' 
        if self.which_columns:
            # can likely get rid of these assertions once sanity checking is formally implemented.
            assert(self.alnlen%3 == 0), "Are you sure this is an alignment? Number of columns is not multiple of three."
            for i in range(1, len(self.seqs)):
                assert ( len(self.seqs[i]) == self.alnlen),  "Are you sure this is an alignment? Number of columns differs among sequences."           
            
            # Loop in increments of 3 for codons
            if self.by == "codon":
                for col in self.which_columns:
                    start = col*3
                    for row in self.seqs:
                        self._full_sequence += row[start:start+3]
            else:
                for col in self.which_columns:
                    for row in self.seqs:
                        self._full_sequence += row[col]
        else:
            for entry in self.seqs:
                self._full_sequence += entry
        
        # Uppercase and processing.
        self._full_sequence = self._full_sequence.upper()
        if self.by == 'amino':
            self._full_sequence = re.sub(self.AMINO_characters, '', self._full_sequence)              
        else:
            self._full_sequence = re.sub(self.DNA_characters, '', self._full_sequence)
        
        # Quick check to ensure that there are actually sequences to use
        assert( len(self._full_sequence) >= len(self._code[0])), "No sequences from which to obtain equilibrium frequencies!"



    def _generate_byFreqs_nuc_amino(self):
        ''' Function for cases when self.by == nuc or self.by == amino '''
        for i in range(0, len(self._full_sequence)):
            try:
                ind = self._code.index(self._full_sequence[i])
            except:
                raise AssertionError("\n\nYour sequences contain non-canonical genetics. Sorry, I'm quitting!")
            self.byFreqs[ind]+=1
        self.byFreqs = np.divide(self.byFreqs, len(self._full_sequence))


    def _generate_byFreqs_codon(self):
        ''' Function for case when self.by == codon ''' 
        for i in range(0, len(self._full_sequence),3):
            codon = self._full_sequence[i:i+3]
            try:
                ind = self._code.index(codon)
            except:
                if codon in MOLECULES.stop_codons:
                    print "\nThere are stop codons in your dataset. I will ignore these, but you should double check your sequences if this was unexpected!"
                    continue
                else:
                    raise AssertionError("\n\nThere is a non-canonical codon triplet in your sequences. Sorry, I'm quitting!")
            self.byFreqs[ind]+=1
        self.byFreqs = np.divide(self.byFreqs, len(self._full_sequence)/3)


    def _generate_byFreqs(self):
        ''' Crux function extraordinaire! '''
        self._make_seq_list()    
        self._process_seq_list()
        if self.by == 'codon':
            self._generate_byFreqs_codon()
        else:
            self._generate_byFreqs_nuc_amino() 
        if self.constraint < 1.0:
            self._unconstrain_frequencies()        
        
















class EmpiricalFreqs(StateFreqs):
    ''' Return state frequencies for empirical models (ones originally used to develop those models).
        The state frequencies are stored in empiricalMatrices.py
        SUPPORTED:
            1. Amino acid: JTT, WAG, LG
            2. Codon:      ECM(un)rest
            NB: scg05 codon model is supported BUT IT DOES NOT HAVE OWN FREQUENCIES.
    '''
    
    def __init__(self, **kwargs):
        super(EmpiricalFreqs, self).__init__(**kwargs)
        try:
            self.empirical_model = kwargs.get('model', None).lower()
        except KeyError:
            print "Need to specify empirical model to get its freqs."
        

    def calculate_freqs(self):    
        ''' Overwrite of parent class function. Such an overwrite will happen only for the EmpiricalFreqs child class, as calculations are not needed.
            We are merely reading from a file to assign state frequencies.
            Currently, we do not support converting these frequencies to a different alphabet.
        '''
        import empiricalMatrices as em
        try:
            return eval("em."+self.empirical_model+"_freqs")
        except:
            print "Couldn't figure out your empirical matrix specification."
            print "Note that we currently support only the following empirical models:"
            print "Amino acid: JTT, WAG, LG."
            print "Codon:      ECM (restricted or unrestricted)."
            print "I'm quiting :/"
            sys.exit()
