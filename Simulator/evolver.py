import numpy as np
from scipy import linalg
import random as rn

class Node(object):
	def __init__(self):
		self.isBase = False # Assigned to True if base of tree. Note that this should have name=None and BL=None. This will serve as an extra check.
		self.name = None # this can either be None (internal) or a leaf name.
		self.children = [] # list of children, each of which is a node
		self.BL = None # Branch length leading up to node
		self.seq = "" # Sequence can be stored here when simulating




class Evolver(Node):
	def __init__(self, codonFreqs, Q_matrix):
		super(Evolver, self).__init__(codonFreqs, Q_matrix)
		
		#Provided by user
		self._STATE = codonFreqs
		self._Q     = Q_matrix
		
		# Barring intelligent design suddenly taking hold
		self._AMINOS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		self._GENCODE = [['GCA', 'GCC', 'GCG', 'GCT'], ['TGC','TGT'], ['GAC', 'GAT'], ['GAA', 'GAG'], ['TTC', 'TTT'], ['GGA', 'GGC', 'GGG', 'GGT'], ['CAC', 'CAT'], ['ATA', 'ATC', 'ATT'], ['AAA', 'AAG'], ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], ['ATG'], ['AAC', 'AAT'], ['CCA', 'CCC', 'CCG', 'CCT'], ['CAA', 'CAG'], ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'] , ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], ['ACA', 'ACC', 'ACG', 'ACT'], ['GTA', 'GTC', 'GTG', 'GTT'], ['TGG'], ['TAC', 'TAT']]
		self._CODONS  = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
		
	
	def retrieveProb(self, codon):
		''' Given a codon, retrieve its row of probabilities from the probMatrix P(t). '''
		ind = self._CODONS[codon]
		row = self._Q[ind]
		return 
		
		
	def generateCodon(self):
		''' Sample a codon. '''
		r = rn.uniform(0,1)
		i=0
		sum=stateFreqs[i]
		while sum <= r:
			i+=1
			sum+=stateFreqs[i]	
		codNum = rn.randint(0,len(gencode[i])-1)
		
		return gencode[i][codNum]
	
	
	
	
	
	def sim_subtree(self, self.tree, baseSeq="<provide root sequence here>"):
		''' Args: (node we are evolving towards which contains the branch length,  parent sequence) '''
		
		# We are at the base of the tree and there is nothing to evolve (the root sequence was provided)
		if (self.tree.BL is None and self.tree.name is None):
			self.tree.seq = baseSeq # only maybe. Think about this some more.
			continue
		
		self.evolve_branch(self.tree, baseSeq)
				
		# We are at an internal node. Keep evolving
		if len(self.tree.children)>0:
			for self.node in self.tree.children:
				sim_subtree(self.node, self.tree.seq)
				
		# We are at a leaf. Print the resulting sequence
		else: 
			print ">"+self.tree.name+"\n"+self.tree.seq+"\n"
			break
		
		
	
	
	
	def evolve_branch(self, self.node, baseSeq):
		''' Node is the node towards which we evolve. baseSeq is the starting sequence for this branch.'''
		
		# Retrieve branch length
		bl = float( self.node.BL )
		assert (bl > 0), "Branch length is negative. Must be >= 0."
		
		# If there is no branch length then there is nothing to evolve. Attach baseSeq to node
		if bl == 0:
			self.node.seq = baseSeq
		
		else:
			# Generate probability matrix for evolution along this branch
			probMatrix = linalg.expm( np.multiply(self._Q, bl) )
			
			# Move along baseSeq and evolve
			newSeq = ""
			
			for i in range(0,len(baseSeq),3):
				codon = baseSeq[i:i+3]
				probRow = retrieveProb(probMatrix, codon)
				




































	























