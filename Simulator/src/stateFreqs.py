import os
import re
import numpy as np
import random as rn
from Bio import SeqIO

from misc import Genetics


class StateFreqs(object):
	'''Will return frequencies. '''
	def __init__(self, **kwargs):
		self.type  = kwargs.get('type') # Type of frequencies to RETURN to user. Either amino, codon, nuc.
		self.by    = kwargs.get('by', self.type) # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. If nucleotide, well, yeah.
		self.molecules  = Genetics()
		self.aminoFreqs = np.zeros(20)
		self.codonFreqs = np.zeros(61)
		self.nucFreqs   = np.zeros(4)
		self.zero = 1e-10



	def sanityByType(self):
		''' Confirm that by and type are compatible, and reassign as needed. '''
		if self.by == 'nuc' and self.type != 'nuc' and self.type is not None:
			print "CAUTION: If calculations are performed with nucleotides, you can only retrieve nucleotide frequencies."
			print "I'm going to calculate nucleotide frequencies for you."
			self.type = 'nuc'
		if self.by == 'nuc' and self.type == 'amino':
			print "CAUTION: Amino acid frequencies cannot be calculated from nucleotide frequencies."
			print "I'm going to calculate your frequencies using amino acid frequencies."
			self.by = 'amino'
		if self.by == 'amino' and self.type == 'nuc':
			print "CAUTION: Nucleotide frequencies cannot be calculated from amino acid frequencies."
			print "I'm going to calculate nucleotide frequencies for you."
			self.by = 'nuc'
		assert(self.type is not None), "I don't know what type of frequencies to calculate! I'm quitting."

	def setCodeLength(self):
		''' Set the codes and lengths once all, if any, "by" issues are resolved ''' 
		if self.by == 'amino':
			self.code = self.molecules.amino_acids
		elif self.by == 'codon':
			self.code = self.molecules.codons
		elif self.by == 'nuc':
			self.code = self.molecules.nucleotides
		self.length = len(self.code)
		
	def generate(self):
		''' BASE CLASS. NOT IMPLEMENTED. '''  
	
	def amino2codon(self):
		''' Calculate codon frequencies from amino acid frequencies. CAUTION: assumes equal frequencies for synonymous codons!! '''
		count = 0
		for codon in self.molecules.codons:
			ind = self.molecules.amino_acids.index(self.molecules.codon_dict[codon])	
			if codon in self.molecules.genetic_code[ind]:
				numsyn = float(len(self.molecules.genetic_code[ind]))
				self.codonFreqs[count] = self.aminoFreqs[ind]/numsyn
			count += 1
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
		''' Calculate the nucleotide frequencies from the codon frequencies. ''' 
		self.generate() # This will get us the codon frequencies. Now convert those to nucleotide
		self.nucFreqs = np.zeros(4) ## ACGT
		for i in range(61):
			codon_freq = self.codonFreqs[i]
			codon = self.molecules.codons[i]
			for n in range(4):
				nuc =  self.molecules.nucleotides[n]
				nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
				if nuc_freq > 0 :
					self.nucFreqs[n] += codon_freq * nuc_freq
		assert( abs(np.sum(self.nucFreqs) - 1.) < self.zero), "Nucleotide state frequencies improperly generated. Do not sum to 1." 


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

	def calcFreqs(self, save=False, savefile=None):
		''' Calculate and return state frequencies. type = the type of frequency to return (nuc, codon, amino). 
			Users can save to file if they would like. If a name for this file is not provided by user still wants to save, a default name is applied in fxn save2file.
		'''
		
				
		# Generate state frequencies for whatever the 'by' is. If the 'type' is different, convert below before returning.
		self.sanityByType()
		self.setCodeLength()
		self.generate() 
		
		if self.type == 'codon':
			if self.by == 'amino':
				self.amino2codon()
			return self.codonFreqs	
		elif self.type == 'amino':
			if self.by == 'codon':
				self.codon2amino()
			return self.aminoFreqs	
		elif self.type == 'nuc':
			if self.by == 'codon':
				self.codon2nuc()
			return self.nucFreqs
		else:
			raise AssertionError("The final type of frequencies you want must be either amino, codon, or nucleotide. I don't know which to calculate, so I'm quitting.")
		
		# As needed.
		if save:
			self.save2file(savefile)	
			

	def save2file(self, savefile):
		if savefile is None:
			savefile = self.type+'_equilibrium_frequencies.txt'
		if self.type == 'codon':
			np.savetxt(savefile, self.codonFreqs)
		elif self.type == 'amino':
			np.savetxt(savefile, self.aminoFreqs)
		elif self.type == 'nuc':
			np.savetxt(savefile, self.nucFreqs)
		else:
			raise AssertionError("This error should seriously NEVER HAPPEN. If it does, someone done broke everything. Please email Stephanie.")



class EqualFreqs(StateFreqs):
	def __init__(self, 	**kwargs):
		super(EqualFreqs, self).__init__(**kwargs)

	def generate(self):
		freqs = np.array(np.repeat(1./float(self.length), self.length))
		assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		self.assignFreqs(freqs)
					
		
					
class RandFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(RandFreqs, self).__init__(**kwargs)

	def generate(self):
		freqs = np.zeros(self.length)
		max = 2. / (self.length) # times 2 for ease/speed/slightly less ridiculousness.
		sum = 0.
		for i in range(int(self.length) - 1):
			freq = rn.uniform(0,max)
			while ((freq==0) or (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			freqs[i] = freq
		freqs[-1] = (1.-sum)	
		assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		self.assignFreqs(freqs)
		
	

class UserFreqs(StateFreqs):
	''' Assign frequencies based on user input. Assume that if not specified, the frequency is zero. 
		Note that 'by' should correspond to the sort of frequencies that they've entered. 'type' should correspond to what they want at the end.
		For instance, it is possible to provide amino acid frequencies and ultimately obtain codon frequencies (with synonymous treated equally, in this circumstance).
	'''
	def __init__(self, **kwargs):
		super(UserFreqs, self).__init__(**kwargs)	
		self.givenFreqs = kwargs.get('freqs', {}) # Dictionary of desired frequencies.	
	
	
	def generate(self):
		freqs = np.zeros(self.length)
		for i in range(self.length):
			element = self.code[i]
			if element in self.givenFreqs:
			 	freqs[i] = self.givenFreqs[element]
		assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly converted from provided to internal object. Do not sum to 1." 
		self.assignFreqs(freqs)
		


########## NEEDS EXTRAORDINARY AMOUNTS OF DEBUGGING ############## 5/10/14.
class ReadFreqs(StateFreqs):
	''' Retrieve frequencies from a file. Can either do global or specify a particular column/group of columns ''' 
	def __init__(self, **kwargs):
		super(ReadFreqs, self).__init__(**kwargs)
		self.seqfile     = kwargs.get('file', None)   # Can also read frequencies from a sequence file
		self.format      = kwargs.get('format', 'fasta') # Default for that file is fasta
		self.whichCol    = kwargs.get('columns', None)     # Which columns we are collecting frequencies from. Default is all columns combined. IF YOU GIVE IT A NUMBER, INDEX AT 0!!!!
	
	def setUpSeqs(self):
		''' Set up sequences and relevent variables for frequency collection. '''
		self.seqs = []
		self.numseq = len(self.rawrecords)
		self.alnlen = len(self.rawrecords[0]) # This will only come into play if we're collecting columns.
		for entry in self.rawrecords:
			self.seqs.append(str(entry.seq))
		print self.seqs
						
	'''
	def getSeq(self):
		'''''' If we want columns, we must get a string of the specific columns we're collecting from.
			Otherwise, we can just turn the whole alignment into a single string.
		''''''
		seq = ''
		if self.by == "codon":
			print "i'm a codon"

		elif self.by == "amino":
			if self.whichCol:
				for col in self.whichCol:
					for row in self.seqs:
						seq += row[col]
			else:
				for entry in self.seqs:
					seq += entry
			#Remove ambig, nonstandard, gaps
			seq = seq.upper()
			seq = seq.translate(None, '-?.*BJOUXZ')
		return seq
		

	def generate(self):
	
		seq = self.getSeq()	
		
		if self.by == 'codon':
			for i in range(0, len(seq),3):
				codon = seq[i:i+3]
				## check for stop codons.
				ind = self.code.index(codon)
				self.codonFreqs[ind]+=1
			self.codonFreqs = np.divide(self.codonFreqs, len(seq)/3)
		
		elif self.by == 'amino':
			for i in range(0, len(seq)):
				ind = self.code.index(seq[i])
				self.aminoFreqs[ind]+=1
			self.aminoFreqs = np.divide(self.aminoFreqs, len(seq))		
		
	'''
		
		
		
		
		
		
		
		
		
		
		
		
		
		
