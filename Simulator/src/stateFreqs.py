import os
import numpy as np
import random as rn
from Bio import AlignIO

from misc import Genetics


class StateFreqs(object):
	'''Will return codon frequencies'''
	def __init__(self, **kwargs):
		self._type     = kwargs.get('type', 'amino') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid.
		self._save     = kwargs.get('save', True) # save statefreqs or not
		self._savefile = kwargs.get('savefile', 'stateFreqs.txt') #default file for saving 
		
		self._codonFreqs = np.zeros(61)
		self._aminoFreqs = np.zeros(20)
		
		# Check file if one was provided
		if self._file:
			assert (os.path(file)), ("Alignment file,", file, ", does not exist. Check path?")

		# Set length based on type
		if self._type == 'amino':
			self._length = 20.
		elif self._type == 'codon':
			self._length = 61.
		else:
			raise AssertionError("type must be either 'amino' or 'codon'")
		
		
		
	def amino2codon(self):
		count = 0
		for codon in molecules.codons:
			ind = molecules.amino_acids.index(molecules.codon_dict[codon])		
			if codon in molecules.genetic_code[ind]
				numsyn = len(molecules.genetic_code[ind])
				self._codonFreqs[count] = self._aminoFreqs[ind]/numsyn
			count += 1
						
				
	def codon2amino(self):
		for a in len(molecules.amino_acids):
			codons1 = molecules.genetic_code[a]
			for c in codons1:
				ind = molecules.codons.index(c)
				self._aminoFreqs[a] += self._codonFreqs[ind]
	
	def save2file(self):
		np.savetxt(self._savefile, self._codonFreqs)
	
			
	def setFreqs(self):
		return 0
	
	def generate(self):
		return 0
		
	def writeFreqs(self)
		return
	

class ReadFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(ReadFreqs, self).__init__(**kwargs)
		self._file     = kwargs.get('file', None) # Can also read frequencies from an alignment file
		self._format   = kwargs.get('format', 'fasta') # Default for that file is fasta
	
	def setFreqs(self):
		aln = AlignIO.read(self._file, self._format)
		bigSeq = ''
		for entry in aln:
			bigSeq += str(entry.seq)
		# Remove ambig
		bigSeq = bigSeq.translate(None, '-?NX') #remove all gaps and ambiguous
	

		if self._type == 'codon':
			for i in range(0, len(bigSeq),3):
				codon = bigSeq[i:i+3]
				ind = molecules.codons.index(codon)
				self._stateFreqs[ind]+=1
			self._codonFreqs = np.divide(self._codonFreqs, len(bigSeq)/3)

		elif self._type == 'amino':
			for i in range(0, len(bigSeq)):
				ind = molecules.amino_acids.index(bigSeq[i])
				self._aminoFreqs[ind]+=1
			self._aminoFreqs = np.divide(self._aminoFreqs, len(bigSeq))
			self.amino2codon()
			
		return self._codonFreqs
	
	
class EqualFreqs(StateFreqs):
	def __init__(self, 	**kwargs):
		super(EqualFreqs, self).__init__(**kwargs)
		
	def setFreqs(self):
		eqFreqs=np.zeros(self._length)
		for i in range(int(self._length)):
			eqFreqs[i] = 1./self._length
		if self._type == 'amino':
			eqFreqs = self.amino2codon(eqFreqs)
		return self._codonFreqs
	
	def generate(self):
		
		
		
class RandFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(RandFreqs, self).__init__(**kwargs)
		
	def setFreqs(self):
		if self._type == 'codon':
			self._codonFreqs = generate()
		
		if self._type == 'amino':
			self._aaFreqs = generate()
			self.amino2codon()	
		return self._codonFreqs
		
	def generate(self):
		randFreqs = np.zeros(self._length)
		freq=float(1)
		max=0.5 # meh
		sum=float(0)
		for i in range(int(self._length)):
			freq = rn.uniform(0,max)
			while ((freq!=0) & (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			randFreqs[i] = freq
		randFreqs[-1] = (1.-sum)	
		return randFreqs
		
