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
		self.debug = kwargs.get('debug', False) # debug mode. some printing.
		self.molecules  = Genetics()
		self.aminoFreqs = np.zeros(20)
		self.codonFreqs = np.zeros(61)
		self.nucFreqs   = np.zeros(4)
		self.zero = 1e-10
		



	def sanityByType(self):
		''' Confirm that by and type are compatible, and reassign as needed. '''
		if self.by == 'nuc' and self.type != 'nuc' and self.type is not None:
			if self.debug:
				print "CAUTION: If calculations are performed with nucleotides, you can only retrieve nucleotide frequencies."
				print "I'm going to calculate nucleotide frequencies for you."
			self.type = 'nuc'
		if self.by == 'nuc' and self.type == 'amino':
			if self.debug:
				print "CAUTION: Amino acid frequencies cannot be calculated from nucleotide frequencies."
				print "I'm going to calculate your frequencies using amino acid frequencies."
			self.by = 'amino'
		if self.by == 'amino' and self.type == 'nuc':
			if self.debug:
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
		''' Calculate and return state frequencies.			
			State frequencies are calculated for whatever "by specifies. If "type" is different, convert before returning. 
			Users can save to file if they would like. If a name for this file is not provided by user still wants to save, a default name is applied in fxn save2file.
		'''
		self.sanityByType()
		self.setCodeLength()
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
			return2user = self.codonFreqs
		elif self.type == 'nuc':
			if self.by == 'codon':
				self.codon2nuc()
			return2user = self.codonFreqs
		else:
			raise AssertionError("The final type of frequencies you want must be either amino, codon, or nucleotide. I don't know which to calculate, so I'm quitting.")
		print "here?????????"
		# As needed.
		if save:
			print "am i getting here"
			self.save2file(savefile)	
		return return2user	

	def save2file(self, savefile):
		print "what about ehre"
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
		return freqs
					
		
					
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
		return freqs
	

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
		return freqs
		


########## NEEDS EXTRAORDINARY AMOUNTS OF DEBUGGING ############## 5/10/14.
class ReadFreqs(StateFreqs):
	''' Retrieve frequencies from a file. Can either do global or specify a particular column/group of columns ''' 
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
			if self.by == "codon":	
				# Can probably get rid of this assertion later when implement parsing/sanity class.
				assert(self.alnlen%3 == 0), "Are you sure this is an alignment? Number of columns is not multiple of three."
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
		if self.by == 'codon' or self.by == 'nuc':
			self.fullSeq = re.sub(self.keepDNA, '', self.fullSeq)
		else:
			self.fullSeq = re.sub(self.keepPROT, '', self.fullSeq)
		
		# Quick check to ensure that there are actually sequences to use
		if self.by == 'codon':
			assert( len(self.fullSeq) >=3 ), "No sequences from which to obtain equilibrium frequencies!"
		else:
			assert( len(self.fullSeq) >=1 ), "No sequences from which to obtain equilibrium frequencies!"
		
	def generate(self):
	
		# Create fullSeq (a single string) for frequency calculations. 
		self.makeSeqList()	
		self.processSeqList()

		freqs = np.zeros(self.length)
		if self.by == 'codon': # loop in triplets for codon data
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
			freqs = np.divide(freqs, len(self.fullSeq)/3)
		else: #loop in increments of 1 for amino and nucleotide data
			for i in range(0, len(self.fullSeq)):
				try:
					ind = self.code.index(self.fullSeq[i])
				except:
					raise AssertionError("Your sequences contain non-canonical genetics. Sorry, I'm quitting!")
				freqs[ind]+=1
			freqs = np.divide(freqs, len(self.fullSeq))		
		return freqs
		
		
		
		
		
		
		
		
		
		
