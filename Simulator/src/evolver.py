import numpy as np
from scipy import linalg
import random as rn
import misc
import stateFreqs
import matrixBuilder

		
class Evolver(object):
	def __init__(self, partitions, tree, outfile):
		
		#Provided by user
		self.PARTS    = partitions
		self.NUMPARTS = len(partitions)
		self.OUTFILE  = outfile
		
		self.SEQLEN   = 0
		for i in range(self.NUMPARTS):
			self.SEQLEN += self.PARTS[i][0]
				
		#Internals
		self.ALNDICT = {}
		self.ZERO = 1e-8
		self.molecules = misc.Genetics()

	def codon2int(self, codon):
		''' Take a codon and return its integer index 0-60 '''
		codind = self.molecules.codons.index(codon)
		return codind
	
	def int2codon(self, codind):
		''' Take a codon index (0-60) and return its corresponding codon ''' 
		codon = self.molecules.codons[codind]
		return codon
	
	def intseq_to_string(self, intseq):
		''' Take a sequence coded as ints and turn to actual codon string '''
		stringseq = ''
		for i in intseq:
			codon = self.int2codon(i)
			stringseq += codon
		return stringseq
	
		
	def generateRootSeq(self):
		rootSeq = np.empty(self.SEQLEN, dtype=int)
		index=0
		for i in range(self.NUMPARTS):
			seqlen = self.PARTS[i][0]
			freqs  = self.PARTS[i][1].stateFreqs
			for j in range(seqlen):
				rootSeq[index] = self.generateCodon(freqs)
				index += 1
		return rootSeq	
	

	def sim_sub_tree(self, tree, baseSeq = None):
		''' Traverse the tree and simulate. '''
		
		# We are at the base and must generate root sequence
		if (baseSeq is None):
			print "Evolving along tree"
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
	
	
	def generateCodon(self, probArray):
		''' Sample a codon. probArray can be any list/numpy array of probabilities that sum to 1.'''
		#### CHECKED FXN ON 2/6/14. WORKS AS INTENDED #####
		# Assertion is overkill but who cares
		assert ( abs(np.sum(probArray) - 1.) < self.ZERO), "Probabilities do not sum to 1. Cannot generate a codon."
		r = rn.uniform(0,1)
		i=0
		sum=probArray[i]
		while sum < r:
			i+=1
			sum+=probArray[i]
		return i		

	def checkBranch(self, node, baseSeq):
		''' Check that the baseSeq exists and the branch length is reasonable ''' 
		## Check that there is a sequence to evolve from
		assert (baseSeq != None), "There is no parent sequence."
		
		# Retrieve branch length
		bl = float( node.branch )
		assert (bl >= 0), "Branch length is negative. Must be >= 0."
		return bl			


	def writeAlignment(self):
		''' Write resulting alignment to a file'''
		print "Writing alignment to file"
		out_handle=open(self.OUTFILE, 'w')
		for entry in self.ALNDICT:
			seq = self.intseq_to_string(self.ALNDICT[entry])
			out_handle.write(">"+entry+"\n"+seq+"\n")
		out_handle.close()	
		
		
		
	def evolve_branch(self, node, baseSeq):
		'''Base class function. Not implemented.'''
		return 0
	


class StaticEvolver(Evolver):
	''' Evolve according to a static landscape (no temporal variation) ''' 
	def __init__(self, *args):
		super(StaticEvolver, self).__init__(*args)
	
	
	def evolve_branch(self, node, baseSeq):
		
		bl = self.checkBranch(node, baseSeq)
		
		# If there is no branch length then there is nothing to evolve. Attach baseSeq to node
		if bl < self.ZERO:
			print bl, "branch length of 0 detected"
			node.seq = baseSeq
		
		else:
			## Evolve for each partition and then join together
			newSeq = np.empty(self.SEQLEN, dtype=int)
			index = 0
			for i in range(self.NUMPARTS):
			
				# set the length and the instantaneous rate matrix for this partition
				seqlen  = self.PARTS[i][0]
				instMat = self.PARTS[i][1].Q
				
				# Generate probability matrix for evolution along this branch and assert correct
				Qt = np.multiply(instMat, bl) # Matrix has already been scaled properly.
				probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
				for i in range(61):
					assert( abs(np.sum(probMatrix[i]) - 1.) < self.ZERO ), "Row in P(t) matrix does not sum to 1."
	
				# Move along baseSeq and evolve. 
				for j in range(seqlen):
					newSeq[index] = self.generateCodon( probMatrix[baseSeq[index]] )
					index+=1
					
			# Attach final sequence to node
			node.seq = newSeq 
			
			
#### CREATED BUT NOT DIFFERENT AT THIS TIME FROM STATICEVOLVER CLASS ####
class ShiftingEvolver(Evolver):
	''' Evolve according to a changing landscape (temporal variation) ''' 
	def __init__(self, *args):
		super(ShiftingEvolver, self).__init__(*args)

	def evolve_branch(self, node, baseSeq):
		
		bl = self.checkBranch(node, baseSeq)
		
		# If there is no branch length then there is nothing to evolve. Attach baseSeq to node
		if bl < self.ZERO:
			print bl, "branch length of 0 detected"
			node.seq = baseSeq
		
		else:
			## Evolve for each partition and then join together
			newSeq = np.empty(self.SEQLEN, dtype=int)
			index = 0
			for i in range(self.NUMPARTS):
			
				# set the length and the instantaneous rate matrix for this partition
				seqlen  = self.PARTS[i][0]
				instMat = self.PARTS[i][1].Q
				
				# Generate probability matrix for evolution along this branch and assert correct
				Qt = np.multiply(instMat, bl) # Matrix has already been scaled properly.
				probMatrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
				for i in range(61):
					assert( abs(np.sum(probMatrix[i]) - 1.) < self.ZERO ), "Row in P(t) matrix does not sum to 1."
	
				# Move along baseSeq and evolve. 
				for j in range(seqlen):
					newSeq[index] = self.generateCodon( probMatrix[baseSeq[index]] )
					index+=1
					
			# Attach final sequence to node
			node.seq = newSeq 







