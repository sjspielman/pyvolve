import os
import re
import numpy as np
import random as rn
from Bio import AlignIO

from misc import Genetics


class StateFreqs(object):
	'''Will return frequencies. '''
	def __init__(self, **kwargs):
		self.by       = kwargs.get('by', 'amino') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid.
		self.save     = kwargs.get('save', True) # save statefreqs or not
		self.savefile = kwargs.get('savefile', 'stateFreqs.txt') #default file for saving 
		self.molecules = Genetics()
		self.codonFreqs = np.zeros(61)
		self.aminoFreqs = np.zeros(20)
		self.zero = 1e-10

		# Set length based on type
		if self.by == 'amino':
			self.code = self.molecules.amino_acids
		elif self.by == 'codon':
			self.code = self.molecules.codons
		else:
			raise AssertionError("type must be either 'amino' or 'codon'")
		self.length = len(self.code)

		
	def amino2codon(self):
		count = 0
		for codon in self.molecules.codons:
			ind = self.molecules.amino_acids.index(self.molecules.codon_dict[codon])	
			if codon in self.molecules.genetic_code[ind]:
				numsyn = float(len(self.molecules.genetic_code[ind]))
				self.codonFreqs[count] = self.aminoFreqs[ind]/numsyn
			count += 1
						
				
	def codon2amino(self):
		for a in range(len(self.molecules.amino_acids)):
			codons1 = self.molecules.genetic_code[a]
			for c in codons1:
				ind = self.molecules.codons.index(c)
				self.aminoFreqs[a] += self.codonFreqs[ind]
	
	def getCodonFreqs(self):
		self.setFreqs()
		return self.codonFreqs


	def getAminoFreqs(self):
		self.setFreqs()
		return self.aminoFreqs				
			
			
	def setFreqs(self):
		if self.by == 'codon':
			self.codonFreqs = self.generate()
			self.codon2amino()
		elif self.by == 'amino':
			self.aminoFreqs = self.generate()
			self.amino2codon()
	
	def save2file(self):
		np.savetxt(self.savefile, self.codonFreqs)
	
	
	def generate(self):
		return 0
	
	
	def getNucFreqs(self):
		''' Gets the nucleotide frequencies from the codon frequencies. Really silly to do this way instead of user providing, but whatever for now. ''' 
		self.setFreqs() # This will get us the codon frequencies. Now convert those to nucleotide
		self.nucFreqs = np.zeros(4) ## ACGT
		for i in range(61):
			codon_freq = self.codonFreqs[i]
			codon = self.molecules.codons[i]
			for n in range(4):
				nuc =  self.molecules.nucleotides[n]
				nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
				if nuc_freq > 0 :
					self.nucFreqs[n] += codon_freq * nuc_freq
		return self.nucFreqs
		
	

class ReadFreqs(StateFreqs):
	''' Retrieve frequencies from a file. Can either do global specify a particular column/group of columns ''' 
	def __init__(self, **kwargs):
		super(ReadFreqs, self).__init__(**kwargs)
		self.alnfile     = kwargs.get('alnfile', None) # Can also read frequencies from an alignment file
		self.format      = kwargs.get('format', 'fasta') # Default for that file is fasta
		
		# Make sure the file exists and that it also has NUCLEOTIDE (not amino acid) data!!
		assert (os.path.exists(self.alnfile)), ("Alignment file,", self.alnfile, ", does not exist. Check path?")
		### TO DO: verify that we are giving it nucleotide data
		
		self.whichCol  = kwargs.get('which', None) # Which columns we are collecting frequencies from. Default is all columns combined. IF YOU GIVE IT A NUMBER, INDEX AT 0!!!!
		
		# make sure self.whichCol is a list. If it's an integer, convert it.
		if type(self.whichCol) is int:
			self.whichCol = [self.whichCol]
	
		## Set up input alignment for use
		tempaln = AlignIO.read(self.alnfile, self.format)
		self.aln = [] 
		self.numseq = len(tempaln)
		self.alnlen = len(tempaln[0]) 
		
		for entry in tempaln:
			self.aln.append(str(entry.seq))
		
	
	def getSeq(self):
		''' Creates a string of the specific columns we are collecting frequencies from '''
		seq = ''
		if self.by == "codon":
			assert(self.alnlen%3 == 0), "\n\nAre you sure this is a codon alignment? Number of columns is not multiple of three."

			if self.whichCol:
				for col in self.whichCol:
					start = col*3
					for row in self.aln:
						seq += row[start:start+3]
			else:
				for entry in self.aln:
					seq += entry
			
			# Remove ambiguities and gaps
			seq = seq.upper()
			seq = re.sub('[^ACGT]', '', seq)
		
		elif self.by == "amino":
			if self.whichCol:
				for col in self.whichCol:
					for row in self.aln:
						seq += row[col]
			else:
				for entry in self.aln:
					seq += entry
			#Remove ambig, nonstandard, gaps
			seq = seq.upper()
			seq = seq.translate(None, '-?.*BJOUXZ')
		return seq
		

	def generate(self):
		seq = self.getSeq()	
		freqs = np.zeros(self.length)
		if self.by == 'codon':
			for i in range(0, len(seq),3):
				codon = seq[i:i+3]
				ind = self.code.index(codon)
				freqs[ind]+=1
			self.freqs = np.divide(self.freqs, len(seq)/3)
		elif self.by == 'amino':
			for i in range(0, len(seq)):
				ind = self.code.index(seq[i])
				freqs[ind]+=1
			freqs = np.divide(freqs, len(seq))
		assert(np.sum(freqs) - 1 < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		return freqs





class EqualFreqs(StateFreqs):
	def __init__(self, 	**kwargs):
		super(EqualFreqs, self).__init__(**kwargs)

	def generate(self):
		freqs = np.array(np.repeat(1./float(self.length), self.length))
		assert(np.sum(freqs) - 1 < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		return freqs
					
class RandFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(RandFreqs, self).__init__(**kwargs)
		
	def generate(self):
		freqs = np.zeros(self.length)
		freq=float(1)
		max=0.05
		sum=float(0)
		for i in range(int(self.length)):
			freq = rn.uniform(0,max)
			while ((freq==0) or (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			freqs[i] = freq
		freqs[-1] = (1.-sum)	
		assert(np.sum(freqs) - 1 < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		return freqs
		
		
class UserFreqs(StateFreqs):
	''' Assign frequencies based on user input. Assume that if not specified, the frequency is zero. '''
	def __init__(self, **kwargs):
		super(UserFreqs, self).__init__(**kwargs)	
		self.givenFreqs = kwargs.get('freqs', {}) # Dictionary of desired frequencies. Example, if by='codon', could provide 
		
		###### Check that user provided frequencies correctly ####### 
		# They should be three letters for codons and one for amino, and they should have actual corresponding molecules.
		# Provided frequencies must also sum to 1.
		if self.by == 'codon':
			keylen = 3
		elif self.by == 'amino':
			keylen = 1
		
		sum = 0
		for entry in self.givenFreqs:
			self.givenFreqs[entry.upper()] = self.givenFreqs.pop(entry) # Ensure upper-case key
			assert ( len(entry) == keylen ), ("\n\nThis key,", entry, "is an unaccepted format. Please use three-letter codes for codons and one-letter codes for amino acids. No ambiguities or stop codons allowed.")
			assert ( entry.upper() in self.code ), ("\n\nThis key,", entry, "is not part of the genetic code. Remember, no ambiguous genetic code is allowed.")
			sum += self.givenFreqs[entry.upper()]
		assert ( abs(sum - 1.) < self.zero), ("\n\nIf you provide frequencies, they must sum to 1. The provided frequencies sum to",sum,".")
		
			
	def generate(self):
		freqs = np.zeros(self.length)
		for i in range(self.length):
			element = self.code[i]
			if element in self.givenFreqs:
			 	freqs[i] = self.givenFreqs[element]
		assert(np.sum(freqs) - 1 < self.zero), "State frequencies improperly provided. Do not sum to 1." 
		return freqs	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
