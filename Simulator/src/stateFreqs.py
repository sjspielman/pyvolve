import os
import numpy as np
import random as rn
from Bio import AlignIO

from misc import Genetics


class StateFreqs(object):
	'''Will return codon frequencies'''
	def __init__(self, **kwargs):
		self.by       = kwargs.get('by', 'amino') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid.
		self.save     = kwargs.get('save', True) # save statefreqs or not
		self.savefile = kwargs.get('savefile', 'stateFreqs.txt') #default file for saving 
		self.molecules = Genetics()
		self.codonFreqs = np.zeros(61)
		self.aminoFreqs = np.zeros(20)

		# Set length based on type
		if self.by == 'amino':
			self.length = 20.
		elif self.by == 'codon':
			self.length = 61.
		else:
			raise AssertionError("type must be either 'amino' or 'codon'")
		
		
		
	def amino2codon(self):
		count = 0
		for codon in self.molecules.codons:
			ind = self.molecules.amino_acids.index(self.molecules.codon_dict[codon])	
			if codon in self.molecules.genetic_code[ind]:
				numsyn = float(len(self.molecules.genetic_code[ind]))
				self.codonFreqs[count] = self.aminoFreqs[ind]/numsyn
			count += 1
						
				
	def codon2amino(self):
		for a in len(self.molecules.amino_acids):
			codons1 = self.molecules.genetic_code[a]
			for c in codons1:
				ind = self.molecules.codons.index(c)
				self.aminoFreqs[a] += self.codonFreqs[ind]
	
	def save2file(self):
		np.savetxt(self.savefile, self.codonFreqs)
	
	
	def getCodonFreqs(self):
		self.setFreqs()
		return self.codonFreqs
	def getAminoFreqs(self):
		self.setFreqs()
		return self.aminoFreqs
	
	def getNucFreqs(self):
		''' Gets the nucleotide frequencies from the codon frequencies ''' 
		self.setFreqs() # This will get us the codon frequencies. Now convert those to nucleotide
		self.nucFreqs = np.zeros(4) ## ACGT
		for i in range(61):
			codon_freq = self.codonFreqs[i]
			codon = self.molecules.codons[i]
			print codon
			for n in range(4):
				nuc =  self.molecules.nucleotides[n]
				nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
				if nuc_freq > 0 :
					self.nucFreqs[n] += codon_freq * nuc_freq
				print self.nucFreqs[n]
		return self.nucFreqs
				
	
			
	def setFreqs(self):
		return 0
	
	def generate(self):
		return 0
		
	def writeFreqs(self):
		return
	

class ReadFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(ReadFreqs, self).__init__(**kwargs)
		self.alnfile     = kwargs.get('alnfile', None) # Can also read frequencies from an alignment file
		self.format      = kwargs.get('format', 'fasta') # Default for that file is fasta
		assert (os.path.exists(self.alnfile)), ("Alignment file,", self.alnfile, ", does not exist. Check path?")

		
		
	def setFreqs(self):
		aln = AlignIO.read(self.alnfile, self.format)
		bigSeq = ''
		for entry in aln:
			bigSeq += str(entry.seq).upper()
		# Remove ambig
		bigSeq = bigSeq.translate(None, '-?NXnx') #remove all gaps and ambiguous
		if self.by == 'codon':
			for i in range(0, len(bigSeq),3):
				codon = bigSeq[i:i+3]
				ind = self.molecules.codons.index(codon)
				self.stateFreqs[ind]+=1
			self.codonFreqs = np.divide(self.codonFreqs, len(bigSeq)/3)

		elif self.by == 'amino':
			for i in range(0, len(bigSeq)):
				ind = self.molecules.amino_acids.index(bigSeq[i])
				self.aminoFreqs[ind]+=1
			self.aminoFreqs = np.divide(self.aminoFreqs, len(bigSeq))
			self.amino2codon()

	
	
	
class EqualFreqs(StateFreqs):
	def __init__(self, 	**kwargs):
		super(EqualFreqs, self).__init__(**kwargs)
		
	def setFreqs(self):
		print self.by
		if self.by == 'codon':
			self.codonFreqs = self.generate()
		elif self.by == 'amino':
			self.aminoFreqs = self.generate()
			self.amino2codon()
	
	def generate(self):
		eqFreqs = np.zeros(self.length)
		for i in range(int(self.length)):
			eqFreqs[i] = 1./self.length
		return eqFreqs
			
	
			
			
			
					
class RandFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(RandFreqs, self).__init__(**kwargs)
		
	def setFreqs(self):
		if self.by == 'codon':
			self.codonFreqs = self.generate()
			print self.codonFreqs
		
		if self.by == 'amino':
			self.aminoFreqs = self.generate()
			self.amino2codon()	

		
	def generate(self):
		randFreqs = np.zeros(self.length)
		freq=float(1)
		max=0.05
		sum=float(0)
		for i in range(int(self.length)):
			freq = rn.uniform(0,max)
			while ((freq==0) or (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			randFreqs[i] = freq
		randFreqs[-1] = (1.-sum)	
		return randFreqs
		
