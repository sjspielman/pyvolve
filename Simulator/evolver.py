import numpy as np
from scipy import linalg
import random as rn
import misc


		
class Evolver():
	def __init__(self, seqSize, codonFreqs, Q_matrix, tree):
		
		#Provided by user
		self._SEQLEN = seqSize
		self._STATE  = codonFreqs
		self._Q      = Q_matrix
		
		# Genetics variables
		self._molecules = misc.Genetics()



	def retrieveProbRow(self, probMatrix, codon):
		''' Given a codon, retrieve its row of probabilities from the probMatrix P(t). '''
		row = probMatrix[self._molecules.codons.index(codon)]
		return row
		
		
	def generateCodon(self, probArray):
		''' Sample a codon. probArray can be any list/numpy array of probabilities that sum to 1.'''
		r = rn.uniform(0,1)
		i=0
		sum=probArray[i]
		while sum <= r:
			i+=1
			sum+=probArray[i]
		return self._molecules.codons[i]
	
	
	def generateRootSeq(self):
		rootSeq = ''
		for i in range(self._SEQLEN):
			rootSeq += self.generateCodon(self._STATE)
		print "rootseq:", rootSeq
		return rootSeq	
	
	
	def sim_sub_tree(self, tree, baseSeq=''):
		''' Args: (node we are evolving towards which contains the branch length,  parent sequence) '''
		
		# We are at the base and must generate root sequence
		if (baseSeq==''):
			#print "BASE", "name", tree.name, "branch", tree.branch, "seq", tree.seq, "children", len(tree.children) 
			tree.seq = self.generateRootSeq()
					
		else:
			#print "INTERNAL", "name", tree.name, "branch", tree.branch, "seq", tree.seq, "children", len(tree.children) 
			self.evolve_branch(tree, baseSeq)
				
		# We are at an internal node. Keep evolving
		if len(tree.children)>0:
			for node in tree.children:
				self.sim_sub_tree(node, tree.seq)
				
		# We are at a leaf. Print the resulting sequence
		else: 
			print ">"+tree.name+"\n"+tree.seq+"\n"	
	
	def evolve_branch(self, node, baseSeq):
		''' Node is the node towards which we evolve. baseSeq is the starting sequence for this branch.'''
		
		# Retrieve branch length
		bl = float( node.branch )
		if bl==0:
			print node.name, node.branch
		assert (bl > 0), "Branch length is negative. Must be >= 0."
		
		# If there is no branch length then there is nothing to evolve. Attach baseSeq to node
		if bl == 0:
			node.seq = baseSeq
		
		else:
			# Generate probability matrix for evolution along this branch
			probMatrix = linalg.expm( np.multiply(self._Q, bl) )
			
			# Move along baseSeq and evolve
			newSeq = ""
			for i in range(0,len(baseSeq),3):
				codon = baseSeq[i:i+3]
				probRow = self.retrieveProbRow(probMatrix, codon)
				newSeq += self.generateCodon(probRow)
			
			# Attach final sequence to node
			node.seq = newSeq
				




































	























