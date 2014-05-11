import numpy as np
from misc import Genetics, Model


class MatrixBuilder(object):
	def __init__(self, model):
		
		# Need to be provided by user
		self.params = model.params
		self.molecules = Genetics()
		self.zero  = 1e-10
		
		### COMMENTING OUT FOR NOW, BUT MUST RETURN TO THIS!!!!!!!! #####
		##### THIS ACTUALLY IS OK, BUT I WILL NEED TO WORK THIS CONDITION IN MORE CAREFULLY. INVARIANTS ARE PERMITTED, BUT THEY DON'T REALLY HAVE A MATRIX. THEREFORE THEY NEED TO BE SOMEHOW CODED AS INVARIANT. ########
		# Double check that stateFreqs is ok (more than 1 character). 
		#for entry in self.stateFreqs:
		#	assert (1. - entry > self.zero), "You must permit evolution to occur!! Can't only allow one character at a site."	
	
		
	def isTI(self, source, target):
		''' Returns True for transition, False for transversion.'''
		ti_py = source in self.molecules.pyrims and target in self.molecules.pyrims
		ti_pu = source in self.molecules.purines and target in self.molecules.purines	
		if ti_py or ti_pu:
			return True
		else:
			return False
	
	def buildQ(self):
		''' Builds instantaneous matrix, Q. 
			For nucleotides, self.size = 4. Amino acids, self.size = 20. Codons, self.size = 61.
		'''	
		self.instMatrix = np.ones([self.size, self.size])
		source_count=0
		for s in range(self.size):
			source = self.code[s]
			for t in range(self.size):
				target = self.code[t]
				rate = self.calcInstProb(source, target)				
				self.instMatrix[s][t] = rate
				
			# Fill in the diagonal position so the row sums to 0.
			if np.sum(self.instMatrix[s]) > self.zero: # This check ensures that there are no -0 values in the matrix.
				self.instMatrix[s][s]= -1*(np.sum( self.instMatrix[s] ))
			assert ( np.sum(self.instMatrix[s]) < self.zero ), "Row in matrix does not sum to 0."
		self.scaleMatrix()

		
	def scaleMatrix(self):
		''' Scale the instantaneous matrix Q so -Sum(pi_iQ_ii)=1 (Goldman and Yang 1994). Ensures branch lengths meaningful for evolving. '''
		scale_factor = 0
		for i in range(self.size):
			scale_factor += (self.instMatrix[i][i] * self.stateFreqs[i]) ##### IS THIS OK FOR EMPIRICAL MODELS? CHECK THIS!!!
		scale_factor*=-1.
		self.instMatrix = np.divide(self.instMatrix, scale_factor)
		######## CHECK THAT THE SCALING WORKED OUT ##############
		sum=0.
		for i in range(self.size):
			sum += (self.instMatrix[i][i] * self.stateFreqs[i])
		assert( abs(sum + 1.) <  self.zero ), "Matrix scaling was a bust."
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

class codonModel_MatrixBuilder(MatrixBuilder):	
	''' This class implements functions relevant to constructing codon model instantaneous matrices (Q).
		Note that the GY and MG models are essentially nested versions of MGREV (which could also be MGHKY, really anything), so can include them all here. Nested versions will merely have extraneous variables fixed to 1 (done elsewhere).
		Model citations:
			GY94:      Yang Z. 1998. Likelihood ratio tests for detecting positive selection and application to primate lysozyme evolution. Mol Biol Evol. 15:568–573.
			MG94:      Muse SV, Gaut BS. 1994. A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome. Mol Biol Evol. 11:715–724.
			MG94(REV): Kosakovsky Pond SL, Muse SV. 2005. Site-to-site variation of synonymous substitution rates. Mol Biol Evol. 22:2375–2385.
	'''		
	def __init__(self, model):
		super(codonModel_MatrixBuilder, self).__init__(model)
		self.size = 61
		self.code = self.molecules.codons
		# PARAMETERS: alpha, beta, mu (this is a dictionary with keys as sourcetarget str, eg "AG")
		# Kappa is not needed. When assigning where I do that later, just make sure that the mu's for transitions are double.
	
	
	def calcInstProb(self, source, target):
		''' Calculate instantaneous probabilities for codon model matrices.	''' 
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
			if self.isSyn(source, target):
				return self.syn(source, target, mydiff[0], mydiff[1])
			else:
				return self.nonSyn(source, target, mydiff[0], mydiff[1])
	
	def getCodonFreq(self, codon):
		''' Get the frequency for a given codon. '''
		return self.stateFreqs[self.molecules.codons.index(codon)]
		
	def syn(self, source_codon, target_codon, source_nuc, target_nuc ):
		''' Calculate the probability of synonymous change.  '''
		return ( self.getCodonFreq(target) * self.params.alpha * self.params.mu[source_nuc + target_nuc] )
	
	
	def nonSyn(self, source, target, source_nuc, target_nuc ):
		''' Calculate the probability of synonymous change. TI/TV issues are taken care of via the mutation parameters. '''
		return ( self.getCodonFreq(target) * self.params.beta * self.params.mu[source_nuc + target_nuc] )

	
	
			