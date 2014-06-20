import os
import re
import numpy as np
import random as rn
from Bio import SeqIO

from misc import Genetics


class StateFreqs(object):
    '''Will return frequencies. '''
    def __init__(self, **kwargs):
        self.type       = kwargs.get('type') # Type of frequencies to RETURN to user. Either amino, codon, nuc.
        self.by         = kwargs.get('by', self.type) # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. If nucleotide, well, yeah.
        self.debug      = kwargs.get('debug', False) # debug mode. some printing. Can likely be removed once parser and more formal sanity checks are implemented.
        self.savefile   = kwargs.get('savefile', None) # for saving the equilibrium frequencies to a file
        self.constraint = kwargs.get('constraint', 1.0) # Constrain provided amino acids to be a certain percentage of total equilbrium frequencies. This allows for non-zero propensities throughout, but non-preferred will be exceptionally rare. Really only used for ReadFreqs and UserFreqs
        self.bias       = kwargs.get('codonBias', None) # To implement codon bias, can provide a decimal giving the percent usage of the preferred state. NOTE: CURRENTLY THE PREFERRED STATE IS RANDOMLY CHOSEN.
        
        self.molecules   = Genetics()
        self.aminoFreqs  = np.zeros(20)
        self.codonFreqs  = np.zeros(61)
        self.nucFreqs    = np.zeros(4)
        self.zero        = 1e-10

        # Set up immediately
        self.sanityByType()
        self.setCodeLength()
        if self.bias is not None:
            assert(self.zero < self.bias <= 1.0), "Codon bias must be >0, <=1."
            assert(self.by == 'amino' and self.type == 'codon'), "If you want codon bias, must have by amino, type codon. Otherwise, I ignore."

    def sanityByType(self):
        ''' Confirm that by and type are compatible, and reassign as needed. 
            RULES:
                1. by=amino      any type
                2. by=codon      any type
                3. by=nuc        type = nuc
               =============================================              
                1. type=amino    by = codon, amino
                2. type=codon    by = codon, amino
                3. type=nuc      by = any by        
        '''
        # This case must raise an assertion error.
        assert(self.by =='amino' or self.by == 'codon' or self.by == 'nuc'), "Codon, amino, or nuc by only!"
        assert(self.type =='amino' or self.type == 'codon' or self.type == 'nuc'), "Codon, amino, or nuc type only!"
        if self.type == 'amino' or self.type == 'codon':
            assert(self.by == 'amino' or self.by == 'codon'), "Incompatible by! For amino acid or codon frequencies, calculations must use either amino acids or codons, NOT nucleotides."
        if self.by == 'amino' or self.by == 'codon':
            assert(self.type == 'amino' or self.type == 'codon'), "Incompatible type! When performing calculations using amino acid or codon frequencies, can return only one of those, NOT nucleotides." 
        if self.by == 'nuc':
            assert(self.type == 'nuc'), "No dice, nuc goes with nuc."
       
  
    def setCodeLength(self):
        ''' Set the codes and lengths once all, if any, "by" issues are resolved ''' 
        if self.by == 'amino':
            self.code = self.molecules.amino_acids
        elif self.by == 'codon':
            self.code = self.molecules.codons
        elif self.by == 'nuc':
            self.code = self.molecules.nucleotides
        self.size = len(self.code)
        
       
    def unconstrainFreqs(self, freqs):
        ''' This function will allow for some frequency constraints to be lessened.
            FUNCTION MAY BE USED BY USERFREQS AND READFREQS ONLY.
            If the constraint value is 0.95, then the preferred (non-zero frequency) entries should only sum to 0.95.
            The remaining 0.05 will be partitioned equally among the non-preferred (freq = 0) entries.
            Therefore, this function allows for some evolutionary "wiggle room" while still enforcing a strong preference.
            
            NB: MAY NOT BE USED IN CONJUNCTION WITH POSITIONAL NUCLEOTIDE FREQUENCIES.
        '''
        freqs = np.multiply(freqs, self.constraint)
        assert (self.size > np.count_nonzero(freqs)), "All state frequencies are 0! This is problematic for a wide variety of reasons."
        addToZero = float( (1.0 - self.constraint) / (self.size - np.count_nonzero(freqs)) )
        for i in range(len(freqs)):
            if ( abs(freqs[i] - 0.0) < self.zero):
                freqs[i] = addToZero
        assert( abs( np.sum(freqs) - 1.0) < self.zero), "unconstraining frequencies did not work properly - freqs don't sum to 1."
        return freqs
        
        
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

   
    
    
    ############################# FREQUENCY CONVERSIONS #########################################
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
        ''' Calculate nucleotide frequencies from amino acid frequencies.'''
        self.amino2codon()
        self.codon2nuc()  
    #####################################################################################   

    def assignFreqs(self, freqs):
        ''' For generate() functions when frequencies are created generally, assign to a specific type with this function. '''
        if self.by == 'codon':
            self.codonFreqs = freqs
        elif self.by == 'amino':
            self.aminoFreqs = freqs
        elif self.by == 'nuc':
            self.nucFreqs = freqs
        else:
            raise AssertionError("I don't know how to calculate state frequencies! I'm quitting.")



    def calcFreqs(self):
        ''' Calculate and return state frequencies.            
            State frequencies are calculated for whatever "by specifies. If "type" is different, convert before returning. 
        '''
        freqs = self.generate()
        assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly generated. Do not sum to 1." 
        self.assignFreqs(freqs)
        if self.type == 'codon':
            if self.by == 'amino':
                self.amino2codon()
            return2user = self.codonFreqs       
        elif self.type == 'amino':
            if self.by == 'codon':
                self.codon2amino()
            return2user = self.aminoFreqs       
        elif self.type == 'nuc':
            if self.by == 'codon':
                self.codon2nuc()
            if self.by == 'amino':
                self.amino2nuc()
            return2user = self.nucFreqs
        if self.savefile:
            self.save2file()    
        return return2user    
        
    
    def convert(self, newtype): 
        ''' Convert between types to return to user.
            TO DO: IN THE FUTURE, OVERHAUL CLASS SUCH THAT TYPE IS REMOVED, MAYBE.
        '''
        conv_expr = "self."+self.by+"2"+newtype+"()"
        if newtype == "nuc":
            if np.array_equal(self.nucFreqs, np.zeros(4)):
                eval(conv_expr)
            return self.nucFreqs
        elif newtype == "codon":
            if np.array_equal(self.codonFreqs, np.zeros(61)):
                try:
                    eval(conv_expr)
                except:
                    raise AssertionError("Can't convert to what you've specified. Quitting.")
            return self.codonFreqs
        elif newtype == "amino":
            if np.array_equal(self.aminoFreqs, np.zeros(20)):
                try:
                    eval(conv_expr)
                except:
                    raise AssertionError("Can't convert to what you've specified. Quitting.")
            return self.aminoFreqs
    
    def save2file(self):
        if self.type == 'codon':
            np.savetxt(self.savefile, self.codonFreqs)
        elif self.type == 'amino':
            np.savetxt(self.savefile, self.aminoFreqs)
        elif self.type == 'nuc':
            np.savetxt(self.savefile, self.nucFreqs)
        else:
            raise AssertionError("This error should seriously NEVER HAPPEN. If it does, someone done broke everything. Please email Stephanie.")



    def freq2dict(self):
        ''' Return a dictionary of frequencies, based on self.type .
        '''
        self.freqDict = {}  # based on TYPE
        if self.type == 'amino':
            freqs = self.aminoFreqs
            code = self.molecules.amino_acids
        elif self.type == 'codon':
            freqs = self.codonFreqs
            code = self.molecules.codons
        else:
            freqs = self.nucFreqs
            code = self.molecules.codons
        for i in range(len(code)):
            if freqs[i] == 0.:
                continue
            else:
                self.freqDict[code[i]] = round(freqs[i], 5)
        return self.freqDict
            










class EqualFreqs(StateFreqs):
    ''' Return equal state frequencies. 
        NOTE: THIS IS THE DEFAULT BEHAVIOR.
    '''
    
    def __init__(self,     **kwargs):
        super(EqualFreqs, self).__init__(**kwargs)
        self.restrict = kwargs.get('restrict', self.code) # Default is all allowed
        
        # TO DO: REPLACE THIS WITH A MORE THOROUGH SANITY CHECK.
        ### Includes: type is list. same alphabet as self.by
        if self.restrict is not self.code:
            assert(type(self.restrict) is list),"restriction needs to be a list."
        
    
    def generate(self):
        fillValue = 1./float(len(self.restrict))
        freqs = np.zeros(self.size)
        for entry in self.restrict:
            freqs[self.code.index(entry)] = fillValue
        return freqs
                 
                    
                    
class RandFreqs(StateFreqs):
    ''' Return random state frequencies.
        Will return essentially flat distributions, but with noise.
    '''
    def __init__(self, **kwargs):
        super(RandFreqs, self).__init__(**kwargs)
        self.restrict = kwargs.get('restrict', self.code) # Default is all allowed
        
        # TO DO: REPLACE THIS WITH A MORE THOROUGH SANITY CHECK.
        ### Includes: type is list. same alphabet as self.by
        if self.restrict:
            assert(type(self.restrict) is list),"restriction needs to be a list."


    def generate(self):
		freqs = np.zeros(self.size)
		
		partial_restrict = self.restrict[:-1] # all but last
		max = 2./len(self.restrict)
		min = 1e-5
		sum = 0.
		for entry in partial_restrict:
			freq = rn.uniform(min,max)
			while (sum + freq > 1):
				freq = rn.uniform(min,max)
			sum += freq
			freqs[self.code.index(entry)] = freq
		freqs[self.code.index(self.restrict[-1])] = (1.-sum)	
		return freqs
    
    
    




class UserFreqs(StateFreqs):
    ''' Assign frequencies based on user input. Assume that if not specified, the frequency is zero. 
        Note that 'by' should correspond to the sort of frequencies that they've entered. 'type' should correspond to what they want at the end.
        For instance, it is possible to provide amino acid frequencies and ultimately obtain codon frequencies (with synonymous treated equally, in this circumstance).
        
        NOTE: UNCONSTRAINING IS POSSIBLE HERE.
    
    '''
    def __init__(self, **kwargs):
        super(UserFreqs, self).__init__(**kwargs)    
        self.givenFreqs = kwargs.get('freqs', {}) # Dictionary of desired frequencies.    
        self.checkBy()


    
    def checkBy(self):
        ''' To make sure that self.by is the same alphabet as provided in the dictionary.
            This function will probably eventually be replaced in a parser/sanity check mechanism.
        '''
        keysize = len( str(self.givenFreqs.keys()[0]) ) # Size of first key. All other keys should be the same size as this one. NOTE THAT IF THIS IS REALLY NOT A STRING, IT WILL BE CAUGHT LATER!! Perhaps/definitely this is inelegant, but I'll deal w/ it later.
        assert(keysize == 1 or keysize == 3), "Bad dictionary keys for userfreqs."
        if keysize == 3:
            self.by == 'codon'
        elif keysize == 1:
            if self.type == 'nuc':
                self.by == 'nuc'
            else:
                self.by == 'amino' 
        
    
    def generate(self):
        freqs = np.zeros(self.size)
        for i in range(self.size):
            element = self.code[i]
            if element in self.givenFreqs:
                freqs[i] = self.givenFreqs[element]
        if self.constraint < 1.0:
            freqs = self.unconstrainFreqs(freqs)
        return freqs
   









class ReadFreqs(StateFreqs):
    ''' Retrieve frequencies from a file. Can either do global or specify a particular column/group of columns.
        NOTE: UNCONSTRAINING IS POSSIBLE HERE.
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
        ''' If we want columns, we must get a string of the specific columns we're collecting from.
            Otherwise, we can just turn the whole alignment into a single string.
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
        if self.by != 'amino':
            self.fullSeq = re.sub(self.keepDNA, '', self.fullSeq)              
        else:
            self.fullSeq = re.sub(self.keepPROT, '', self.fullSeq)
        
        # Quick check to ensure that there are actually sequences to use
        assert( len(self.fullSeq) >= len(self.code[0])), "No sequences from which to obtain equilibrium frequencies!"



    def generate_nuc_amino(self, freqs, sequence):
        ''' Function for cases when self.by == nuc or self.by == amino '''
        for i in range(0, len(sequence)):
            try:
                ind = self.code.index(sequence[i])
            except:
                raise AssertionError("Your sequences contain non-canonical genetics. Sorry, I'm quitting!")
            freqs[ind]+=1
        return np.divide(freqs, len(sequence))


    def generate_codon(self, freqs):
        ''' Function for case when self.by == codon ''' 
        for i in range(0, len(self.fullSeq),3):
            codon = self.fullSeq[i:i+3]
            try:
                ind = self.code.index(codon)
            except:
                if codon in self.molecules.stop_codons:
                    if self.debug:
                        print "There are stop codons in your dataset. I will ignore these, but you should double check your sequences if this was unexpected!"
                        continue
                    else:
                        raise AssertionError("There is a non-canonical codon triplet in your sequences. Sorry, I'm quitting!")
            freqs[ind]+=1
        return np.divide(freqs, len(self.fullSeq)/3)



    def generate(self):
        ''' Crux function extraordinaire! '''
        # Create fullSeq (a single string) for frequency calculations. 
        self.makeSeqList()    
        self.processSeqList()

        freqs = np.zeros(self.size)
        if self.by == 'codon':
            freqs = self.generate_codon(freqs)
        else:
            freqs = self.generate_nuc_amino(freqs, self.fullSeq) # Note that we provide an attribute as an argument because that function is also used for generate_posNuc, and for that function we will NOT be providing an attribute.
            
        if self.constraint < 1.0:
            freqs = self.unconstrainFreqs(freqs)        
        return freqs
        
















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
            freqs = eval("em."+self.empiricalModel+"_freqs")
        except:
            print "Couldn't figure out your empirical matrix specification."
            print "Note that we currently support only the following empirical models:"
            print "Amino acid: JTT, WAG, LG."
            print "Codon:      ECM (restricted or unrestricted)."
        return freqs
