#### TO DO IDEAS:
###### Boltzmann as a class
###### Set the GC content, for type codon or nucleotide frequencies ONLY.
###### 

import os
import re
import numpy as np
import random as rn
from Bio import SeqIO
from misc import Genetics


class StateFreqs(object):
    '''Will return frequencies. huzzah.'''
    def __init__(self, **kwargs):
        
        # Generic
        self.zero        = 1e-10
        self.molecules   = Genetics()   
        self.nucFreqs    = np.zeros(4)     
        self.aminoFreqs  = np.zeros(20)
        self.codonFreqs  = np.zeros(61)
        
        # Input options and some sanity checking
        self.by = kwargs.get('by') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. If nucleotide, well, yeah.
        assert(self.by =='amino' or self.by == 'codon' or self.by == 'nuc'), "\n\nYou have either no 'by' or a wrong 'by'. Remember, codon, amino, or nuc only!"
        self.setCodeLength()
        self.byFreqs = np.zeros(self.size)
        
        self.restrict   = kwargs.get('restrict', self.code) # For the equal, rand subclasses only.
        self.constraint = kwargs.get('constraint', 1.0) # For the user, read subclasses only. Constrain provided amino acids to be a certain percentage of total equilbrium frequencies. This allows for non-zero propensities throughout, but non-preferred will be exceptionally rare.
        self.bias       = kwargs.get('codonBias', None) # To implement codon bias, can provide a decimal giving the percent usage of the preferred state. NOTE: CURRENTLY THE PREFERRED STATE IS RANDOMLY CHOSEN.
        self.savefile   = kwargs.get('savefile', None) # for saving the equilibrium frequencies to a file

        if self.bias:
            assert(self.zero < self.bias <= 1.0), "Codon bias must be >0, <=1."
        if self.constraint:
            assert(self.zero <  self.constraint <= 1.0), "Constraint must be >0, <=1."
        if self.restrict is not self.code:
            assert(type(self.restrict) is list), "Restriction must be a list."
        




    def setCodeLength(self):
        ''' Set the codes and lengths ''' 
        if self.by == 'amino':
            self.code = self.molecules.amino_acids
        elif self.by == 'codon':
            self.code = self.molecules.codons
        elif self.by == 'nuc':
            self.code = self.molecules.nucleotides
        self.size = len(self.code)

           
       
    def unconstrainFreqs(self):
        ''' This function will allow for some frequency constraints to be lessened for the self.byFreqs
            FUNCTION MAY BE USED BY USERFREQS AND READFREQS ONLY.
            If the constraint value is 0.95, then the preferred (non-zero frequency) entries should only sum to 0.95.
            The remaining 0.05 will be partitioned equally among the non-preferred (freq = 0) entries.
            Therefore, this function allows for some evolutionary "wiggle room" while still enforcing a strong preference.
        '''
        self.byFreqs = np.multiply(self.byFreqs, self.constraint)
        assert (self.size > np.count_nonzero(self.byFreqs)), "All state frequencies are 0! This is problematic for a wide variety of reasons."
        addToZero = float( (1.0 - self.constraint) / (self.size - np.count_nonzero(self.byFreqs)) )
        for i in range(self.size):
            if ( abs(self.byFreqs[i] - 0.0) < self.zero):
                self.byFreqs[i] = addToZero
        assert( abs( np.sum(self.byFreqs) - 1.0) < self.zero), "unconstraining frequencies did not work properly - freqs don't sum to 1."
        
        
    def codonBias(self, aa_count, syn):
        ''' Implements codon bias. There is a self.bias param which gives a decimal indicating the frequency of the preferred codon.
            Requires by=amino, type=codon
            TO DO, FUTURE DEVELOPMENT: Allow users to provide the preferred codons. As of now, the preferred codon is RANDOM.
            Args: aa_count = the amino acid index we are working with
                  syn      = the list of synonymous codons for this amino acid index        
        '''
        # If we are dealing with a single-codon amino acid (M or W), simply fill the value since codon bias is not possible.
        sum = 0. # for some debugging
        if len(syn) == 1:
            cind = self.molecules.codons.index(syn[0])    
            self.codonFreqs[cind] = self.aminoFreqs[aa_count]
            sum += self.aminoFreqs[aa_count]
        else:
            prefIndex = rn.randint(0, len(syn)-1)
            prefFreq = self.aminoFreqs[aa_count] * self.bias
            nonprefFreq = (self.aminoFreqs[aa_count] - prefFreq)/(len(syn) - 1.)  
            
            for s in range(len(syn)):
                cind = self.molecules.codons.index(syn[s])
                if s == prefIndex:
                    self.codonFreqs[cind] = prefFreq
                    sum += prefFreq
                else:
                    self.codonFreqs[cind] = nonprefFreq
                    sum += nonprefFreq
        assert(abs (sum - self.aminoFreqs[aa_count]) < self.zero), "Codon bias improperly implemented."

   
    
    
    ######################################### FREQUENCY CONVERSIONS #########################################
    def amino2codon(self):
        ''' Calculate codon frequencies from amino acid frequencies. 
            Unless codon bias is specified, will assume equal synonymous frequencies.
        '''
        for aa_count in range(20):
            syn = self.molecules.genetic_code[aa_count]
            if self.bias:
                self.codonBias(aa_count, syn)
            else:
                for synCodon in syn:
                    cind = self.molecules.codons.index(synCodon)
                    self.codonFreqs[cind] = self.aminoFreqs[aa_count]/float(len(syn))
        assert( abs(np.sum(self.codonFreqs) - 1.) < self.zero), "Codon state frequencies improperly generated from amino acid frequencies. Do not sum to 1."                 
                
                
    def codon2amino(self):
        ''' Calculate amino acid frequencies from codon frequencies. ''' 
        for a in range(len(self.molecules.amino_acids)):
            codons1 = self.molecules.genetic_code[a]
            for c in codons1:
                ind = self.molecules.codons.index(c)
                self.aminoFreqs[a] += self.codonFreqs[ind]
        assert( abs(np.sum(self.aminoFreqs) - 1.) < self.zero), "Amino acid state frequencies improperly generated from codon frequencies. Do not sum to 1." 
    
    def codon2nuc(self):
        ''' Calculate nucleotide frequencies from codon frequencies. '''
        for i in range(61):
            codon_freq = self.codonFreqs[i]
            codon = self.molecules.codons[i]
            for n in range(4):
                nuc =  self.molecules.nucleotides[n]
                nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
                if nuc_freq > 0 :
                    self.nucFreqs[n] += codon_freq * nuc_freq
        assert( abs(np.sum(self.nucFreqs) - 1.) < self.zero), "Nucleotide state frequencies improperly generated. Do not sum to 1." 

        
    def amino2nuc(self):
        ''' Calculate nucleotide frequencies from amino acid frequencies. Lazy function, hurray!'''
        self.amino2codon()
        self.codon2nuc()  
    
    #####################################################################################   

    def assign_byFreqs(self):
        ''' Assign self.byFreqs to either amino, codon, or nuc. '''
        if self.by == 'codon':
            self.codonFreqs = self.byFreqs
        elif self.by == 'amino':
            self.aminoFreqs = self.byFreqs
        elif self.by == 'nuc':
            self.nucFreqs = self.byFreqs
        else:
            raise AssertionError("WHAT ARE WE DOING HERE.")


    def calcFreqs(self, **kwargs):
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
        if np.array_equal(self.byFreqs, np.zeros(self.size)):
            self.generate() # generates self.byFreqs       
            assert( abs(np.sum(self.byFreqs) - 1.) < self.zero), "State frequencies improperly generated. Do not sum to 1." 
            self.assign_byFreqs()
        
        # Convert frequencies if needed
        if type != self.by:
            conv_expr = "self."+self.by+"2"+type+"()"
            eval(conv_expr)
        
        # Save if needed
        if save is not None:
            np.savetxt(save, eval("self."+type+"Freqs"), fmt='%.5e')
        return eval("self."+type+"Freqs")




class BoltzmannFreqs(StateFreqs):
    ''' Return state frequencies generated by Boltzmann distribution.
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
            self.ranking = self.code
        else:
            assert( type(self.ranking) is list ), "Ranking must be a full list of amino acids in order."
            assert( sorted(self.ranking) == self.code ), "Your ranking list does not appear to contain all amino acids, or has incorrect letters in it."
    
    def setBoltzmann(self):
        ddg_values = np.random.normal(size = self.size) 
        numer_list = np.zeros(self.size)
        denom = 0.
        for d in range(self.size):
            val = np.exp(-1. * self.factor * ddg_values[d])
            denom += val
            numer_list[d] = val
        return numer_list/(np.sum(numer_list)) 
    
    def generate(self):
        tempFreqs = self.setBoltzmann()
        if self.ranking is not self.code:
            tempFreqs = np.sort(tempFreqs)[::-1]
        count = 0
        for aa in self.ranking:
            self.byFreqs[self.code.index(aa)] = tempFreqs[count]
            count += 1
                





class EqualFreqs(StateFreqs):
    ''' Return equal state frequencies. 
        NOTE: THIS IS THE DEFAULT BEHAVIOR.
    '''
    
    def __init__(self, **kwargs):
        super(EqualFreqs, self).__init__(**kwargs)
    
    def generate(self):
        fillValue = 1./float(len(self.restrict))
        for entry in self.restrict:
            self.byFreqs[self.code.index(entry)] = fillValue                
                    
                    
                    
                    
                    
class RandFreqs(StateFreqs):
    ''' Return random state frequencies.
        Will return essentially flat distributions, but with noise.
    '''
    def __init__(self, **kwargs):
        super(RandFreqs, self).__init__(**kwargs)
      
    def generate(self):
        partial_restrict = self.restrict[:-1] # all but last
        max = 2./len(self.restrict)
        min = 1e-5
        sum = 0.
        for entry in partial_restrict:
            freq = rn.uniform(min,max)
            while (sum + freq > 1):
                freq = rn.uniform(min,max)
            sum += freq
            self.byFreqs[self.code.index(entry)] = freq
        self.byFreqs[self.code.index(self.restrict[-1])] = (1.-sum)    
    
    
    




class UserFreqs(StateFreqs):
    ''' Assign frequencies based on user input. Assume that if not specified, the frequency is zero. 
        Note that 'by' should correspond to the sort of frequencies that they've entered. 'type' should correspond to what they want at the end.
        For instance, it is possible to provide amino acid frequencies and ultimately obtain codon frequencies (with synonymous treated equally, in this circumstance).
        
        NOTE: UNCONSTRAINING IS POSSIBLE HERE.
    
    '''
    def __init__(self, **kwargs):
        super(UserFreqs, self).__init__(**kwargs)    
        self.givenFreqs = kwargs.get('freqs', {}) # Dictionary of desired frequencies.    
        self.checkByKeys() ######## this will likely be removed when formal sanity checking is implemented eventually ######


    def checkByKeys(self):
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
    
    def generate(self):
        for i in range(self.size):
            element = self.code[i]
            if element in self.givenFreqs:
                self.byFreqs[i] = self.givenFreqs[element]
        if self.constraint < 1.0:
            self.unconstrainFreqs()










class ReadFreqs(StateFreqs):
    ''' Retrieve frequencies from a file. Can either do global or specify a particular column/group of columns.
        NOTE: UNCONSTRAINING IS POSSIBLE HERE.
        
        TO DO: SANITY CHECKING WILL NEED TO VERIFY THAT THE PROVIDED SEQUENCE FILE IS IN THE SAME ALPHABET AS THE BY IS.
     ''' 
    def __init__(self, **kwargs):
        super(ReadFreqs, self).__init__(**kwargs)
        
        self.seqfile  = kwargs.get('file', None)   # Can also read frequencies from a sequence file
        self.format   = kwargs.get('format', 'fasta') # Default for that file is fasta
        self.whichCol = kwargs.get('columns', None)     # Which columns we are collecting frequencies from. Default is all columns combined. IF YOU GIVE IT A NUMBER, INDEX AT 0!!!!
        self.seqs     = [] # Sequence records obtained from sequence file
        self.fullSeq  = '' # Single sequence string from which to obtain frequencies
        self.keepDNA  = re.compile(r"[^ACGT]") # DNA regexp for what to keep
        self.keepPROT = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]") # protein regexp for what to keep
       
        
    def makeSeqList(self):
        ''' Set up sequences and relevent variables for frequency collection. '''
        raw = list(SeqIO.parse(self.seqfile, self.format))
        self.seqs = []
        self.numseq = len(raw)
        self.alnlen = len(raw[0]) # This will only come into play if we're collecting columns.
        for entry in raw:
            self.seqs.append(str(entry.seq))  

    
    def processSeqList(self):
        ''' 
            Turns full sequence we want to grab frequencies from into a single string.
            NOTE: If we want columns, we must get a string of the specific columns we're collecting from.
        ''' 
        if self.whichCol:
            # can likely get rid of these assertions once sanity checking is formally implemented.
            assert(self.alnlen%3 == 0), "Are you sure this is an alignment? Number of columns is not multiple of three."
            for i in range(1, len(self.seqs)):
                assert ( len(self.seqs[i]) == self.alnlen),  "Are you sure this is an alignment? Number of columns differs among sequences."           
            
            # Loop in increments of 3 for codons
            if self.by == "codon":
                for col in self.whichCol:
                    start = col*3
                    for row in self.seqs:
                        self.fullSeq += row[start:start+3]
            else:
                for col in self.whichCol:
                    for row in self.seqs:
                        self.fullSeq += row[col]
        else:
            for entry in self.seqs:
                self.fullSeq += entry
        
        # Uppercase and processing.
        self.fullSeq = self.fullSeq.upper()
        if self.by == 'amino':
            self.fullSeq = re.sub(self.keepPROT, '', self.fullSeq)              
        else:
            self.fullSeq = re.sub(self.keepDNA, '', self.fullSeq)
        
        # Quick check to ensure that there are actually sequences to use
        assert( len(self.fullSeq) >= len(self.code[0])), "No sequences from which to obtain equilibrium frequencies!"



    def generate_nuc_amino(self):
        ''' Function for cases when self.by == nuc or self.by == amino '''
        for i in range(0, len(self.fullSeq)):
            try:
                ind = self.code.index(self.fullSeq[i])
            except:
                raise AssertionError("\n\nYour sequences contain non-canonical genetics. Sorry, I'm quitting!")
            self.byFreqs[ind]+=1
        self.byFreqs = np.divide(self.byFreqs, len(self.fullSeq))


    def generate_codon(self):
        ''' Function for case when self.by == codon ''' 
        numstop = 0
        for i in range(0, len(self.fullSeq),3):
            codon = self.fullSeq[i:i+3]
            try:
                ind = self.code.index(codon)
            except:
                if codon in self.molecules.stop_codons:
                    numstop += 3
                    print "\nThere are stop codons in your dataset. I will ignore these, but you should double check your sequences if this was unexpected!"
                    print "stop at pos",i
                    continue
                else:
                    raise AssertionError("\n\nThere is a non-canonical codon triplet in your sequences. Sorry, I'm quitting!")
            self.byFreqs[ind]+=1
        self.byFreqs = np.divide(self.byFreqs, (len(self.fullSeq) - numstop)/3)


    def generate(self):
        ''' Crux function extraordinaire! '''
        self.makeSeqList()    
        self.processSeqList()
        if self.by == 'codon':
            self.generate_codon()
        else:
            self.generate_nuc_amino() 
        if self.constraint < 1.0:
            self.unconstrainFreqs()        
        
















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
            self.empiricalModel = kwargs.get('model', None).lower()
        except KeyError:
            print "Need to specify empirical model to get its freqs."
        

    def calcFreqs(self):    
        ''' Overwrite of parent class function. Such an overwrite will happen only for the EmpiricalFreqs child class, as calculations are not needed.
            We are merely reading from a file to assign state frequencies.
            Currently, we do not support converting these frequencies to a different alphabet.
        '''
        import empiricalMatrices as em
        try:
            return eval("em."+self.empiricalModel+"_freqs")
        except:
            print "Couldn't figure out your empirical matrix specification."
            print "Note that we currently support only the following empirical models:"
            print "Amino acid: JTT, WAG, LG."
            print "Codon:      ECM (restricted or unrestricted)."
            print "I'm quiting :/"
            sys.exit()
