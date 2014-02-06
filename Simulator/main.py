import misc
from treeReader import Newick
import modelBuilder
from evolver import Evolver

import numpy as np
import random as rn

################################################################################
def generateStateFreqs(aaFreqs=None):
	''' Numpy array of randomly-generated reasonable codon frequencies, in alphabetical order.'''
	
	# no amino acid frequencies provided. generate random ones first
	if aaFreqs==None:
		aaFreqs = np.empty(20)
		freq=float(1)
		max=0.1
		sum=float(0)
		for i in range(19):
			freq = rn.uniform(0,max)
			while ((freq!=0) & (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			aaFreqs[i] = freq
		aaFreqs[19] = (1.-sum)
	
	stateFreqs = np.empty(61)
	for c in range(61):
		for a in range(len(molecules.genetic_code)):
			numsyn = float(len(molecules.genetic_code[a]))
			
			# We've found the correct amino acid
			if molecules.codons[c] in molecules.genetic_code[a]:
				stateFreqs[c] = (aaFreqs[a]/numsyn)
				continue
	return stateFreqs
################################################################################
def printTree(tree, level=0):
	indent=''
	for i in range(level):
		indent+='\t'
	print indent, tree.name, tree.branch, tree.seq
	if len(tree.children)>0:
		for node in tree.children:
			print tree.seq
			printTree(node, level+1)	
###################################################################################
	
# get dna/amino/codon lists
molecules = misc.Genetics()

# read in tree
my_tree  = Newick("giant.tre")
my_tree  = my_tree.readTree() # set True to print out the tree

# Build model
mu=1
kappa=1
stateFreqs = generateStateFreqs() # generates random codon state frequencies such that all synonymous codons have same frequency (neutral)
myModel = modelBuilder.SellaModel(mu, kappa, stateFreqs)
Q = myModel.buildQ()

# Evolve along the tree
sequenceLength=50 # number of codons
myEvolver = Evolver(sequenceLength, stateFreqs, Q, my_tree)
myEvolver.sim_sub_tree(my_tree)









				
				
