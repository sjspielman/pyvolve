import numpy as np
import random as rn
from Bio import AlignIO
import sys
sys.path.append("src/")

import misc
from treeReader import Newick
import modelBuilder
from evolver import Evolver



################################################################################
def generateStateFreqs(aaFreqs=None):
	''' Numpy array of randomly-generated reasonable codon frequencies, in alphabetical order.'''
	
	equal=True
	if equal:
		stateFreqs=np.zeros(61)
		for i in range(len(stateFreqs)):
			stateFreqs[i] = 1./61.
		return stateFreqs
	else:
		fix_to_zero=[]
		# no amino acid frequencies provided. generate random ones first
		if aaFreqs==None:
			aaFreqs = np.empty(20)
			freq=float(1)
			max=0.1
			sum=float(0)
			for i in range(19):
				if i in fix_to_zero:
					aaFreqs[i]=0
					continue	
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
					if aaFreqs[a] == 0:
						print c
					stateFreqs[c] = (aaFreqs[a]/numsyn)
					continue
		return stateFreqs

def readStateFreqs(alnfile, alnformat="fasta"):
	''' Read in codon state frequencies from an alignment file '''
	
	bigSeq = ''
	aln = AlignIO.read(alnfile, alnformat)
	for entry in aln:
		bigSeq += str(entry.seq)
	# Remove ambig
	bigSeq = bigSeq.translate(None, '-?N') #remove all gaps and ambiguous
	
	stateFreqs = np.zeros(61)
	for i in range(0, len(bigSeq),3):
		codon = bigSeq[i:i+3]
		ind = molecules.codons.index(codon)
		stateFreqs[ind]+=1
	stateFreqs = np.divide(stateFreqs, len(bigSeq)/3)
	
	return stateFreqs

################################################################################
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

#stateFreqs=readStateFreqs('aln.phy', 'phylip-relaxed')


# read in tree
my_tree  = Newick("trees/10.tre")
my_tree  = my_tree.readTree() # set True to print out the tree

# Build model
mu=1e-4
kappa=5.3
omega=8
stateFreqs = generateStateFreqs()
np.savetxt('stateFreqs.txt', stateFreqs)
myModel = modelBuilder.GY94Model(stateFreqs, kappa, omega)
Q = myModel.buildQ()

# Evolve along the tree
outfile="stephanie.fasta"
sequenceLength=3 # number of codons
myEvolver = Evolver(sequenceLength, stateFreqs, Q, my_tree, outfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeAlignment()




# Super linear runtime as you introduce more sequences
## 0.46 s for 10seq, 500codon
## 1.3 s for 50seq, 500codon
## 2.3 s for 100seq, 500 codon
## 11.2 s for 500seq,500 codon
## 22 s for 1000seq, 500codon



				
				
