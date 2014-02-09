import os
import numpy as np
import random as rn
from Bio import AlignIO

from misc import Genetics


class StateFreqs(object):
	'''Will return codon frequencies'''
	def __init__(self, **kwargs):
		self._type   = kwargs.get('type', 'amino') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid.
		self._file   = kwargs.get('file', None) # Can also read frequencies from an alignment file
		self._format = kwargs.get('format', 'fasta') # Default for that file is fasta
		
		# Check file if one was provided
		if self._file:
			assert (os.path(file)), ("File", file, "does not exist. Check path?")
	

		# Set length based on type
		if self._type == 'amino':
			self._length = 20.
		elif self._type == 'codon':
			self._length = 61.
		else:
			raise AssertionError("type must be either 'amino' or 'codon'")
		
		
		
	def amino2codon(self, aaFreqs):
		for c in range(61):
			for a in range(len(molecules.genetic_code)):
				numsyn = float(len(molecules.genetic_code[a]))
				# We've found the correct amino acid
				if molecules.codons[c] in molecules.genetic_code[a]:
					codonFreqs[c] = (aaFreqs[a]/numsyn)
					continue
		return codonFreqs	
	
	def setFreqs(self):
		return 0
	

class ReadFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(ReadFreqs, self).__init__(**kwargs)
	
	def setFreqs:
		aln = AlignIO.read(self._file, self._format)
		bigSeq = ''
		for entry in aln:
			bigSeq += str(entry.seq)
		# Remove ambig
		bigSeq = bigSeq.translate(None, '-?NX') #remove all gaps and ambiguous
	
		stateFreqs = np.zeros(self._length)
		if self._type == 'codon':
			for i in range(0, len(bigSeq),3):
				codon = bigSeq[i:i+3]
				ind = molecules.codons.index(codon)
				stateFreqs[ind]+=1
			stateFreqs = np.divide(stateFreqs, len(bigSeq)/3)
		elif self._type == 'amino':
			for i in range(0, len(bigSeq)):
				ind = molecules.amino_acids.index(bigSeq[i])
				stateFreqs[ind]+=1
			stateFreqs = np.divide(stateFreqs, len(bigSeq))
			stateFreqs = self.amino2codon(stateFreqs)
		return stateFreqs
	
	
class EqualFreqs(StateFreqs):
	def __init__(self, 	**kwargs):
		super(EqualFreqs, self).__init__(**kwargs)
		
	def setFreqs(self):
		eqFreqs=np.zeros(self._length)
		for i in range(self._length):
			eqFreqs[i] = 1./self._length
		if self._type == 'amino':
			eqFreqs = self.amino2codon(eqFreqs)
		return eqFreqs
		
		
class RandFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(RandFreqs, self).__init__(**kwargs)
		
	def setFreqs(self):
		randFreqs = np.empty(self._length)
		freq=float(1)
		max=0.5 # meh
		sum=float(0)
		for i in range(self._length):
			freq = rn.uniform(0,max)
			while ((freq!=0) & (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			randFreqs[i] = freq
		randFreqs[-1] = (1.-sum)
		if self._type == 'amino':
			randFreqs = self.amino2codon(randFreqs)	
		return randFreqs	