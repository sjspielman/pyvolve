import numpy as np
import random as rn
from Bio import AlignIO
import sys
sys.path.append("src/")

import misc
import newick
import stateFreqs
import matrixBuilder
from evolver import Evolver

# get dna/amino/codon lists
molecules = misc.Genetics()

# read in tree
print "reading tree"
my_tree  = newick.readTree(file="trees/10.tre", show=False) # set True to print out the tree

# Build model
print "building model"
print "frequencies"
fgen = stateFreqs.ReadFreqs(type='codon', alnfile = 'refaln99.fasta', save='stateFreqs.txt')
codonFreqs = fgen.getCodonFreqs()
fgen.save2file()


k=5.3
w=8.
print "matrix"
matBuilder = matrixBuilder.GY94Matrix(codonFreqs, k, w)
Q = matBuilder.buildQ()
print Q
print "evolving"
# Evolve along the tree
outfile="stephanie.fasta"
sequenceLength=3 # number of codons
myEvolver = Evolver(sequenceLength, codonFreqs, Q, my_tree, outfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeAlignment()




# Super linear runtime as you introduce more sequences
## 0.46 s for 10seq, 500codon
## 1.3 s for 50seq, 500codon
## 2.3 s for 100seq, 500 codon
## 11.2 s for 500seq,500 codon
## 22 s for 1000seq, 500codon



				
				
