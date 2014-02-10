######### 2/10/14. IN CASE I GET EXCITABLE. ########

import os
import numpy as np
import random as rn
from Bio import AlignIO

from misc import Genetics


## HIERARCHY <- OK ACTUALLY I HAVE NO IDEA
## Modeler(type, state_method, savestate, savestate_file, statefile, mu, kappa, omega)
##    StateFreqs(type, save, savefile)
##       EqualFreqs
##       RandFreqs
##       GivenFreqs(statefile)
##    Matrix
##       SellaHirsh(mu, kappa)
##       GY94(omega, kappa)


'''

		# Check file if one was provided
		if self._file:
			assert (os.path(file)), ("Alignment file,", file, ", does not exist. Check path?")

'''

class Modeler(object):
	def __init__(self, **kwargs):
		self._type      = kwargs.get('type', 'amino') # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. default amino.
		self._method    = kwargs.get('state_method', None) # To generate frequencies as equal, random, or from a file. default, equal
		self._matrix    = kwargs.get('matrix', 'GY94') # model matrix. default, GY94 
		self._savestate = kwargs.get('save_state', True) # whether or not to save state frequencies.
		if self._savestate:
			self._savestate_file = kwargs.get('savestate_file', 'codonFreqs.txt'):
		
	

class StateFreqs(Modeler):
	def __init__(self, **kwargs):
		self._codonFreqs = np.zeros(61)
		self._aminoFreqs = np.zeros(20)
		
		# Set length based on type
		if self._type == 'amino':
			self._length = 20.
		elif self._type == 'codon':
			self._length = 61.
		
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
			
	def setFreqs(self):
		return 0
	
	def generate(self):
		return 0
		
	def writeFreqs(self)
		return
	

class FileFreqs(StateFreqs):
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
































class matrixBuilder(object):
	def __init__(self, codonFreqs):
		
		# Need to be provided by user
		self._STATE = codonFreqs

		# Genetics variables		
		self._molecules = Genetics()

	def isTI(self, source, target):
		''' Returns True for transition, False for transversion.'''
		if (source in self._molecules.pyrims and target in self._molecules.purines):
			return False
		else:
			return True	
	
	
	def isSyn(self, source, target):
		''' Given a source codon and target codon, return True if the change is synonymous.'''
		if (self._molecules.codon_dict[source] == self._molecules.codon_dict[target]):
			return True
		else:
			return False
	
	
	def getCodonFreq(self, codon):
		''' Get the frequency for a given codon. '''
		Freq = self._STATE[self._molecules.codons.index(codon)]
		return Freq	
	
		
	def calcMutProb(self, source, target):
		''' Calculates a substitution probability between two codons. If single mutation, return the probability/rate. Else, return 0. ''' 
		mydiff=''
		for i in range(3):
			if source[i] == target[i]:	
				continue
			else:	
				mydiff+=source[i]+target[i]
		
		# Either no change >1 mutations. It's probability is zero. We will correct the diagonal later.	
		if len(mydiff)!=2:
			return 0
		
		# If a single mutational step away, return the probability
		else:		
			# Transitions
			if self.isTI(mydiff[0], mydiff[1]):
				if self.isSyn(source, target):
					return self.synTI(source, target)
				else:
					return self.nonSynTI(source, target)
			# Transversions
			else:
				if self.isSyn(source, target):
					return self.synTV(source, target)
				else:
					return self.nonSynTV(source, target)
				
				
	def buildQ(self):
		''' Builds the 61x61 matrix Q '''
		
		transMatrix = np.ones([61,61]) # Look at me, hardcoding that there are 61 codons!
		source_count=0
		for s in range(61):
			source = self._molecules.codons[s]
			for t in range(61):
				target = self._molecules.codons[t]
				rate = self.calcMutProb(source, target)
				transMatrix[s][t] = rate
			
			# Fill in the diagonal position so the row sums to 0. Confirm.
			transMatrix[s][s]= -1*(np.sum( transMatrix[s] ))
			assert (np.sum(transMatrix[s]==0)), "Row in matrix does not sum to 0."
		
		transMatrix = self.scaleMatrix(transMatrix)
		return transMatrix	
	
	
	
	
	def scaleMatrix(self, mat):
		''' Scale Q matrix so -Sum(pi_iQ_ii)=1 (Goldman and Yang 1994). '''
		scale_factor = 0
		for i in range(61):
			scale_factor += (mat[i][i] * self._STATE[i])
		scale_factor*=-1.
		mat = np.divide(mat, scale_factor)
		
		######## CHECK THAT THE SCALING WORKED OUT ##############
		sum=0.
		for i in range(61):
			sum += (mat[i][i] * self._STATE[i])
		assert(round(sum) == -1.0), "Matrix scaling was a bust."
		return mat		
		
		
		
		

	########################################################## 
	## Base functions for computing rates. Not implemented. ## 
	def synTI(self, source, target):
		return 0
	def synTV(self, source, target):
		return 0
	def nonSynTI(self, source, target):
		return 0
	def nonSynTV(self, source, target):
		return 0
	###########################################################
		
	
		
class SellaModel(Modeler):
	def __init__(self, stateFreqs, mu, kappa):
		''' Implement the Sella (2005) model '''
		super(SellaModel, self).__init__(stateFreqs)
		self._MU = mu
		self._KAPPA = kappa
	
	def fix(self, source_freq, target_freq):
		''' Given pi(i) and pi(j), where pi() is the equilibrium a given codon in that column, return probability_of_fixation_(i->j). '''
		if target_freq == 0 or source_freq == 0:
			return 0 # If either has 0 frequency, we should never reach it.
		elif source_freq == target_freq:
			return 1 # confirmed correct
		else:
			return ( (np.log(target_freq) - np.log(source_freq)) / (1 - source_freq/target_freq) )


	def synTI(self, source, target):
		''' Probability of synonymous transition '''
		return ( self._MU * self._KAPPA )
	
	
	def synTV(self, source, target):
		''' Probability of synonymous tranversion '''
		return ( self._MU )
	
	
	def nonSynTI(self, source, target):
		''' Probability of nonsynonymous transition '''
		sFreq = self.getCodonFreq(source)
		tFreq = self.getCodonFreq(target)
		return ( self._MU * self._KAPPA * self.fix(sFreq, tFreq) )				
	
		
	def nonSynTV(self, source, target):
		''' Probability of nonsynonymous tranversion '''
		sFreq = self.getCodonFreq(source)
		tFreq = self.getCodonFreq(target)
		return ( self._MU * self.fix(sFreq, tFreq) )	

				


class GY94Model(Modeler):
	def __init__(self, stateFreqs, kappa, omega):
		'''Implement the GY94 model '''
		super(GY94Model, self).__init__(stateFreqs)
		self._KAPPA = kappa
		self._OMEGA = omega
		
	
	def synTI(self, source, target):
		''' Probability of synonymous transition '''
		return ( self.getCodonFreq(target) * self._KAPPA )
	
	
	def synTV(self, source, target):
		''' Probability of synonymous tranversion '''
		return ( self.getCodonFreq(target) )
	
	
	def nonSynTI(self, source, target):
		''' Probability of nonsynonymous transition '''
		return ( self.getCodonFreq(target) * self._KAPPA * self._OMEGA )				
	
		
	def nonSynTV(self, source, target):
		''' Probability of nonsynonymous tranversion '''
		return ( self.getCodonFreq(target) * self._OMEGA )	

				
				
				
				
				
				
			
			