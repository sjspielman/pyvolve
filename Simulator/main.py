import numpy as np
import random as rn
from Bio import AlignIO
import sys
sys.path.append("src/")

import misc
import getState
import newick
from modelBuilder import Modeler
from evolver import Evolver

# get dna/amino/codon lists
molecules = misc.Genetics()

# state frequencies
#state = getState.

# read in tree
my_tree  = newick.readTree(file="trees/10.tre", show=False) # set True to print out the tree

# Build model
f = getState.EqualFreqs(type='codon', save='stateFreqs.txt')
stateFreqs = f.setFreqs()
np.savetxt('stateFreqs.txt', stateFreqs)


mu=1e-4
kappa=5.3
omega=8
myModel = GY94Model(stateFreqs, kappa, omega)
Q = myModel.buildQ()
assert 1==0
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



				
				
