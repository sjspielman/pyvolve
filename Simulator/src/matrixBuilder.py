import numpy as np
from misc import Genetics, Model



class MatrixBuilder(object):
	def __init__(self, model):
		
		# Need to be provided by user
		self.EQFREQS = model.stateFreqs
		self.ZERO  = 1e-8

		# Genetics variables		
		self.molecules = Genetics()

	def isTI(self, source, target):
		''' Returns True for transition, False for transversion.'''
		ti_py = source in self.molecules.pyrims and target in self.molecules.pyrims
		ti_pu = source in self.molecules.purines and target in self.molecules.purines	
		if ti_py or ti_pu:
			return True
		else:
			return False
	
	
	def isSyn(self, source, target):
		''' Given a source codon and target codon, return True if the change is synonymous.'''
		if (self.molecules.codon_dict[source] == self.molecules.codon_dict[target]):
			return True
		else:
			return False
	
	
	def getCodonFreq(self, codon):
		''' Get the frequency for a given codon. '''
		Freq = self.EQFREQS[self.molecules.codons.index(codon)]
		return Freq	
	
	
	def calcMutProb(self, source, target):
		''' Calculate instantaneous probabilities based on this model ''' 
		mydiff=''
		for i in range(3):
			if source[i] == target[i]:	
				continue
			else:	
				mydiff+=source[i]+target[i]
		
		# Either no change, >1 mutations. We will correct the diagonal later.	
		if len(mydiff)!=2:
			return 0
		else:
			return ( self.getProb(mydiff, source, target) )
			
			
						
	def buildQ(self, paramFlag):
		''' Builds the 61x61 matrix Q. paramFlag specifies the type of parameters we have (kappa, nuc) '''
		
		instMatrix = np.ones([61,61]) # Look at me, hardcoding that there are 61 codons!
		source_count=0
		for s in range(61):
			source = self.molecules.codons[s]
			for t in range(61):
				target = self.molecules.codons[t]
				
				rate = self.calcMutProb(source, target)				
				instMatrix[s][t] = rate
			
			# Fill in the diagonal position so the row sums to 0. Confirm.
			instMatrix[s][s]= -1*(np.sum( instMatrix[s] ))
			assert ( np.sum(instMatrix[s]) < self.ZERO ), "Row in matrix does not sum to 0."
		
		instMatrix = self.scaleMatrix(instMatrix)
		return instMatrix	
	
	
	
	
	def scaleMatrix(self, mat):
		''' Scale Q matrix so -Sum(pi_iQ_ii)=1 (Goldman and Yang 1994). '''
		scale_factor = 0
		for i in range(61):
			scale_factor += (mat[i][i] * self.EQFREQS[i])
		scale_factor*=-1.
		mat = np.divide(mat, scale_factor)
		
		######## CHECK THAT THE SCALING WORKED OUT ##############
		sum=0.
		for i in range(61):
			sum += (mat[i][i] * self.EQFREQS[i])
		assert( abs(sum + 1.) <  self.ZERO ), "Matrix scaling was a bust."
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
	def syn(self, source, target):
		return 0
	def nonsyn(self, source, target):
		return 0
	###########################################################
		
	
		
class SellaMatrixKappa(MatrixBuilder):
	''' Implement the Sella (2005) model '''
	def __init__(self, model):
		super(SellaMatrix, self).__init__(model)
		self.params = model.params
		self.MU     = model.params["mu"]
		self.KAPPA  = model.params["kappa"]
	
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
		return ( self.MU * self.KAPPA )
	
	
	def synTV(self, source, target):
		''' Probability of synonymous tranversion '''
		return ( self.MU )
	
	
	def nonSynTI(self, source, target):
		''' Probability of nonsynonymous transition '''
		sFreq = self.getCodonFreq(source)
		tFreq = self.getCodonFreq(target)
		return ( self.MU * self.KAPPA * self.fix(sFreq, tFreq) )				
	
		
	def nonSynTV(self, source, target):
		''' Probability of nonsynonymous tranversion '''
		sFreq = self.getCodonFreq(source)
		tFreq = self.getCodonFreq(target)
		return ( self.MU * self.fix(sFreq, tFreq) )	


class GY94MatrixKappa(MatrixBuilder):
	'''Implement the GY94 model '''
	def __init__(self, model):
		super(GY94Matrix, self).__init__(model)
		self.params = model.params
		self.OMEGA  = model.params["omega"]
		self.KAPPA  = model.params["kappa"]
		
	
	def synTI(self, source, target):
		''' Probability of synonymous transition '''
		return ( self.getCodonFreq(target) * self.KAPPA )
	
	def synTV(self, source, target):
		''' Probability of synonymous tranversion '''
		return ( self.getCodonFreq(target) )
	
	def nonSynTI(self, source, target):
		''' Probability of nonsynonymous transition '''
		return ( self.getCodonFreq(target) * self.KAPPA * self.OMEGA )				
		
	def nonSynTV(self, source, target):
		''' Probability of nonsynonymous tranversion '''
		return ( self.getCodonFreq(target) * self.OMEGA )	

	def getProb(self, mydiff, source, target):
		''' Calculate instantaneous probabilities for GY94 Matrix ''' 
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


class Rodrigue(MatrixBuilder):
	''' Nicolas Rodrigue's 2010 model. Note that it uses nucleotide frequencies and NOT codon frequencies. '''
	def __init__(self, model):
		super(GY94Matrix, self).__init__(model)
		
		self.params = model.params
		
		### Mutational parameters.
		self.nucMut = {} # note that each pair is ordered alphabetically.
		self.nucMut["AC"] = model.params["AC"]
		self.nucMut["AG"] = model.params["AG"]
		self.nucMut["AT"] = model.params["AT"]
		self.nucMut["CG"] = model.params["CG"]
		self.nucMut["CT"] = model.params["CT"]
		self.nucMut["GT"] = model.params["GT"]
		self.nucFreqs = model.params["nucFreqs"] # state frequencies of nucleotides, in order ACGT
		self.aaVector = model.params["aaVector"] # amino acid propensity vector, in alphabetical order (as in molecules.amino_acids)
			
			
			
	def syn(self, target_nuc, diff):
		''' Probability of synonymous change '''
		return ( self.nucMut[diff] * self.nucFreqs[target_nuc] )
	
	def nonsyn(self, target_nuc, diff, source, target):
		''' Probability of nonsynonymous change '''
		selCoeff = self.getSelCoeff(source, target)
		return ( self.nucMut[diff] * self.nucFreqs[target_nuc] * self.fixationFactor(selfCoeff) )

	def getSelCoeff(self, source, target)
	
		# Amino acid propensities for source and target codons.
		aa_source_prop = self.aaVector[ self.amino_acids.index( self.molecules.codon_dict[source] ) ]
		aa_target_prop = self.aaVector[ self.amino_acids.index( self.molecules.codon_dict[target] ) ]
		
		# If either of them has a 0 propensity it should never occur
		if aa_source_prop == 0 or aa_target_prop == 0:
			return 0
		
		# If equal propensities then there should be no selective pressure
		elif aa_source_prop == aa_target_prop:
			return 1. #checked. Will simplify this way in the end.
		
		else:
			return np.log( self.aa_target_prop / self.aa_source_prop )
	
	def fixationFactor(self, selCoeff):
		return (selCoeff / (1. - np.exp(selCoeff) )
			
		
		
	def calcMutProb(self, mydiff, source, target):
		''' Calculate instantaneous probabilities for Rodrigue's model (similar to NielsenYang2008) ''' 
		
		sorted_diff = "".join(sorted(mydiff)) # Alphabetize the source/target nucleotides in order to easily retrieve the mutation probability
		
		if self.isSyn(mydiff[0], mydiff[1]):
			return self.syn( mydiff[1], sorted_diff )
		else:
			return self.nonsyn( mydiff[1], sorted_diff, source, target )
		
