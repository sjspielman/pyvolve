import numpy as np
from scipy import linalg
from Bio.Seq import Seq
from Bio.Alphabet import *
import re

class Model(object):
	def __init__(self, mu, kappa, aminoFreqs):
		
		# Need to be provided by user
		self._MU = mu
		self._KAPPA = kappa
		self._STATE = aminoFreqs
		
		# Unless we have a dramatic and unfortunate shift in the genetic code....
		self._PYR     = ['C', 'T']
		self._PUR     = ['A', 'G']
		self._AA_CODE = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		self._GENCODE = [['GCA', 'GCC', 'GCG', 'GCT'], ['TGC','TGT'], ['GAC', 'GAT'], ['GAA', 'GAG'], ['TTC', 'TTT'], ['GGA', 'GGC', 'GGG', 'GGT'], ['CAC', 'CAT'], ['ATA', 'ATC', 'ATT'], ['AAA', 'AAG'], ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], ['ATG'], ['AAC', 'AAT'], ['CCA', 'CCC', 'CCG', 'CCT'], ['CAA', 'CAG'], ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'] , ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], ['ACA', 'ACC', 'ACG', 'ACT'], ['GTA', 'GTC', 'GTG', 'GTT'], ['TGG'], ['TAC', 'TAT']]
		self._CODONS  = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
		
	
	def isTI(self, source, target):
		''' Returns True for transition, False for transversion.'''
		if (source in self._PYR and target in self._PUR):
			return False
		else:
			return True	
	
	
	def isSyn(self, source, target):
		''' Given a source codon and target codon, return True if the change is nonsynonymous.'''
		if (self._CODONS[source] == self._CODONS[target]):
			return True
		else:
			return False
	
	
	def getCodonFreq(self, codon):
		''' Get the frequency for a given codon. This function assumes that all synonymous codons have the same frequency. '''
		aa = self._CODONS[codon]
		Freq = self._STATE[self._AA_CODE.index(aa)]
		return Freq	
	
		
	def calcMutProb(self, source, target):
		''' Calculates a substitution probability between two codons. If single mutation, return the probability/rate. Else, return 0. ''' 
		mydiff=''
		count=0
		for i in range(3):
			if source[i] == target[i]:	
				continue
			else:	
				mydiff+=source[i]+target[i]
		
		# Either no change >1 mutations. It's probability is zero.	
		if len(mydiff)!=2:
			return 0
		
		# If a single mutational step away, return the probability
		else:		
			if self.isTI(mydiff[0], mydiff[1]):
				
				if self.isSyn(source, target):
					return self.synTI(source, target)
				else:
					return self.nonSynTI(source, target)
			else:
				
				if self.isSyn(source, target):
					return self.synTV(source, target)
				else:
					return self.nonSynTV(source, target)
				
				
	def buildCodonTransitionMatrix(self):
		''' Builds the matrix Q '''
		
		transMatrix = np.empty([61,61]) # Look at me, hardcoding that there are 61 codons!
		source_count=0
		
		for source in self._CODONS:
			target_count=0
			for target in self._CODONS:
				rate = self.calcMutProb(source, target)
				assert (rate>=0 and rate<1), "Fixation probability is exceedingly incorrect."
				transMatrix[source_count][target_count] = rate
				target_count+=1
				
			# Fill in the diagonal position so the row sums to 0
			transMatrix[source_count][source_count]= -1*(np.sum( transMatrix[source_count] ))
			assert (np.sum(transMatrix[source_count]==0)), "Row in matrix does not sum to 0."
			source_count+=1
		
		return transMatrix	

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
		
	
		
class SellaModel(Model):
	def __init__(self, *args, **kwargs):
		super(SellaModel, self).__init__(*args, **kwargs)
	
	def fix(self, source_freq, target_freq):
		''' Given pi(i) and pi(j), where pi() is the equilibrium a given codon in that column, return probability_of_fixation_(i->j). '''
		if source_freq == target_freq:
			return 1 # THIS IS WRONG. FIX!!!
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


	
			
			
	
	
	
	
	

codonFreqs = [0.006035686218971376, 0.006035686218971376, 0.006035686218971376, 0.006035686218971376, 0.04041378723808315, 0.04041378723808315, 0.032098342321156166, 0.032098342321156166, 0.043116404375596995, 0.043116404375596995, 0.009440547545140699, 0.009440547545140699, 0.011632846243747412, 0.011632846243747412, 0.011632846243747412, 0.011632846243747412, 0.04573777513398951, 0.04573777513398951, 0.024498759238302317, 0.024498759238302317, 0.024498759238302317, 0.009148227493706296, 0.009148227493706296, 0.013470358817488972, 0.013470358817488972, 0.013470358817488972, 0.013470358817488972, 0.013470358817488972, 0.013470358817488972, 0.004875895626235194, 0.01762993572322193, 0.01762993572322193, 0.004788541107242323, 0.004788541107242323, 0.004788541107242323, 0.004788541107242323, 0.04826387436585986, 0.04826387436585986, 0.0036601640920194563, 0.0036601640920194563, 0.0036601640920194563, 0.0036601640920194563, 0.0036601640920194563, 0.0036601640920194563, 0.0032299254578469906, 0.0032299254578469906, 0.0032299254578469906, 0.0032299254578469906, 0.0032299254578469906, 0.0032299254578469906, 0.022423097043176427, 0.022423097043176427, 0.022423097043176427, 0.022423097043176427, 0.01446788129594178, 0.01446788129594178, 0.01446788129594178, 0.01446788129594178, 0.016631129371868093, 0.02687200552651542, 0.02687200552651542]
myModel=SellaModel(1e-2, 4.5, codonFreqs)

Q = myModel.buildCodonTransitionMatrix()
bl = 0.353212809381
Qt = np.multiply(mat, bl)


			

				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		
	
	
	
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
	
	
