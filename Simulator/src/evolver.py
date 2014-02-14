import numpy as np
from scipy import linalg
import random as rn
import misc


		
class Evolver():
	def __init__(self, seqSize, codonFreqs, Q_matrix, tree, outfile):
		
		#Provided by user
		self.SEQLEN  = seqSize
		self.STATE   = codonFreqs
		self.Q       = Q_matrix
		self.OUTFILE = outfile
		
		#Internals
		self.ALNDICT = {}
		self.ACCURACY = 1e-8
		
		# Genetics variables
		self.molecules = misc.Genetics()


	def codon2int(self, codon):
		''' Take a codon and return its integer index 0-60 '''
		codind = self.molecules.codons.index(codon)
		return codind
	
	def int2codon(self, codind):
		''' Take a codon index (0-60) and return its corresponding codon ''' 
		codon = self.molecules.codon[codind]
		return codon
	
	def intseq_to_string(self, intseq):
		''' Take a sequence coded as ints and turn to actual codon string '''
		stringseq = ''
		for i in intseq:
			codon = self.int2codon(i)
			stringseq += codon
		return stringseq
	
		
	def generateCodon(self, probArray):
		''' Sample a codon. probArray can be any list/numpy array of probabilities that sum to 1.'''
		#### CHECKED FXN ON 2/6/14. WORKS AS INTENDED #####
		# Assertion is overkill but who cares
		assert ( abs(np.sum(probArray) - 1.) < self.ACCURACY), "Probabilities do not sum to 1. Cannot generate a codon."
		
		r = rn.uniform(0,1)
		i=0
		sum=probArray[i]
		while sum < r:
			i+=1
			sum+=probArray[i]
		return i
	
	
	def generateRootSeq(self):
		rootSeq = []
		for i in range(self.SEQLEN):
			rootSeq.append( self.generateCodon(self.STATE) )
		return rootSeq	
	
	
	def sim_sub_tree(self, tree, baseSeq=''):
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
			self.ALNDICT[tree.name]=tree.seq
			
					
	
	def evolve_branch(self, node, baseSeq):
		''' Node is the node towards which we evolve. baseSeq is the starting sequence for this branch (parent's sequence).'''
		
		## Check that there is a sequence to evolve from
		assert (baseSeq != None), "There is no parent sequence."
		
		# Retrieve branch length
		bl = float( node.branch )
		assert (bl >= 0), "Branch length is negative. Must be >= 0."
		
		# If there is no branch length then there is nothing to evolve. Attach baseSeq to node
		if bl < self.ACCURACY:
			node.seq = baseSeq
		
		else:
			# Generate probability matrix for evolution along this branch and assert correct
			Qt = np.multiply(self.Q, bl) # Matrix has already been scaled properly.
			probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
			
			for i in range(61):
				assert( abs(np.sum(probMatrix[i]) - 1.) > self.ACCURACY ), "Row in P(t) matrix does not sum to 1."
	
			# Move along baseSeq and evolve
			newSeq = []
			for codind in baseSeq:
				probRow = probMatrix[codind]
				newSeq.append( self.generateCodon(probRow) )
			# Attach final sequence to node
			node.seq = newSeq


	def writeAlignment(self):
		''' Write resulting alignment to a file'''
		out_handle=open(self.OUTFILE, 'w')
		for entry in self.ALNDICT:
			seq = self.intseq_to_string(entry)
			out_handle.write(">"+seq+"\n"+self.ALNDICT[entry]+"\n")
		out_handle.close()
				




































	























