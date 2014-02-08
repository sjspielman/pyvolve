import numpy as np
from scipy import linalg
import random as rn
import misc


		
class Evolver():
	def __init__(self, seqSize, codonFreqs, Q_matrix, tree, outfile):
		
		#Provided by user
		self._SEQLEN  = seqSize
		self._STATE   = codonFreqs
		self._Q       = Q_matrix
		self._OUTFILE = outfile
		
		#Internals
		self._ALNDICT = {}
		
		# Genetics variables
		self._molecules = misc.Genetics()



	def retrieveProbRow(self, probMatrix, codon):
		''' Given a codon, retrieve its row of probabilities from the probMatrix P(t). '''
		row = probMatrix[self._molecules.codons.index(codon)]
		return row
		
		
	def generateCodon(self, probArray):
		''' Sample a codon. probArray can be any list/numpy array of probabilities that sum to 1.'''
		#### CHECKED FXN ON 2/6/14. WORKS AS INTENDED #####
		# Assertion is overkill but who cares
		assert (round(np.sum(probArray)) == 1.), "Probabilities do not sum to 1. Cannot generate a codon."
		
		r = rn.uniform(0,1)
		i=0
		sum=probArray[i]
		while sum < r:
			i+=1
			sum+=probArray[i]
		return self._molecules.codons[i]
	
	
	def generateRootSeq(self):
		rootSeq = ''
		for i in range(self._SEQLEN):
			rootSeq += self.generateCodon(self._STATE)
		return rootSeq	
	
	
	def sim_sub_tree(self, tree, baseSeq='', final=''):
		''' Traverse the tree and simulate. '''
		
		# We are at the base and must generate root sequence
		if (baseSeq==''):
			tree.seq = self.generateRootSeq()		
		else:
			self.evolve_branch(tree, baseSeq)
				
		# We are at an internal node. Keep evolving
		if len(tree.children)>0:
			for node in tree.children:
				self.sim_sub_tree(node, tree.seq)
				
		# We are at a leaf. Save the final sequence
		else: 
			self._ALNDICT[tree.name]=tree.seq
			
					
	
	def evolve_branch(self, node, baseSeq):
		''' Node is the node towards which we evolve. baseSeq is the starting sequence for this branch (parent's sequence).'''
		
		## Check that there is a sequence to evolve from
		assert (baseSeq !=''), "There is no parent sequence."
		
		# Retrieve branch length
		bl = float( node.branch )
		assert (bl >= 0), "Branch length is negative. Must be >= 0."
		
		# If there is no branch length then there is nothing to evolve. Attach baseSeq to node
		if round(bl) == 0.:
			node.seq = baseSeq
		
		else:
			# Generate probability matrix for evolution along this branch and assert correct
			Qt = np.multiply(self._Q, bl) # Matrix has already been scaled properly.
			probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
			for i in range(61):
				assert( round(np.sum(probMatrix[i])) == 1.0 ), "Row in P(t) matrix does not sum to 1."
	
					
			# Move along baseSeq and evolve
			newSeq = ""
			for i in range(0,len(baseSeq),3):
				codon = baseSeq[i:i+3]
				probRow = self.retrieveProbRow(probMatrix, codon) #checked. it does retrieve the correct codon row
				newSeq += self.generateCodon(probRow) 
			# Attach final sequence to node
			node.seq = newSeq
	
	
	def writeAlignment(self):
		''' Write resulting alignment to a file'''
		out_handle=open(self._OUTFILE, 'w')
		for entry in self._ALNDICT:
			out_handle.write(">"+entry+"\n"+self._ALNDICT[entry]+"\n")
		out_handle.close()
				




































	























