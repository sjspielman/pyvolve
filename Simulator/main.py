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
my_tree  = newick.readTree(file="trees/20.tre", show=False) # set True to print out the tree

# Build model
print "building model"
print "frequencies"
fgen = stateFreqs.EqualFreqs(by='amino', alnfile='hrh1.fasta', save='stateFreqs.txt')
codonFreqs = fgen.getCodonFreqs()
fgen.save2file()


k=.5
w=8
matBuilder = matrixBuilder.GY94Matrix(codonFreqs, k, w)
Q = matBuilder.buildQ()
print "evolving"
# Evolve along the tree
outfile="stephanie.fasta"
sequenceLength=100 # number of codons
myEvolver = Evolver(sequenceLength, codonFreqs, Q, my_tree, outfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeAlignment()
