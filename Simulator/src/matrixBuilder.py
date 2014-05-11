import numpy as np
from misc import Genetics, Model



class MatrixBuilder(object):
	def __init__(self, model):
		
		# Need to be provided by user
		self.stateFreqs = model.params["stateFreqs"]
		self.zero  = 1e-10
		
		### COMMENTING OUT FOR NOW, BUT MUST RETURN TO THIS!!!!!!!! #####
		##### THIS ACTUALLY IS OK, BUT I WILL NEED TO WORK THIS CONDITION IN MORE CAREFULLY. INVARIANTS ARE PERMITTED, BUT THEY DON'T REALLY HAVE A MATRIX. THEREFORE THEY NEED TO BE SOMEHOW CODED AS INVARIANT. ########
		# Double check that stateFreqs is ok (more than 1 character). 
		#for entry in self.stateFreqs:
		#	assert (1. - entry > self.zero), "You must permit evolution to occur!! Can't only allow one character at a site."	

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
	
	# This function may or may not have to be generalized to accomodate additional types (aa, nuc, codon). OR, they could each have their own function.
	# TBD until overall organization considered.
	def getCodonFreq(self, codon):
		''' Get the frequency for a given codon. '''
		Freq = self.stateFreqs[self.molecules.codons.index(codon)]
		return Freq	
		
	
	########################################################## 
	## Base class functions for computing rates. Not implemented. ## 
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
	

	
	
	
class codonModel_MatrixBuilder(MatrixBuilder):
	''' Functions for building the instantaneous probability matrix, Q, for codon models of sequence evolution. 
		All implemented codon models assume that only a single mutation is allowed at a given codon along a branch.
		MODELS: GY, MG, ECM_w_k. The final model is empirical codon model plus kappa, omega parameters. 
	'''
	def __init__(self, model):
		super(codonModel_MatrixBuilder, self).__init__(model)
		# INSERT HERE: any general model parameters that all codon models share.
		self.kappa  = model.params["kappa"]
	
	def calcInstProb(self, source, target):
		''' Calculate instantaneous probabilities for substituting source -> target. 
			This function is codon-specific, and will consider only circumstances for a single nucleotide change.
		''' 
		mydiff=''
		for i in range(3):
			if source[i] == target[i]:	
				continue
			else:	
				mydiff+=source[i]+target[i]
		if len(mydiff)!=2:
			return 0
		else:
			return ( self.getProb(mydiff, source, target) )
	