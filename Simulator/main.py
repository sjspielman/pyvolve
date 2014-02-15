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
my_tree, flag_list  = newick.readTree(file="trees/20.tre", show=False) # set True to print out the tree

# If there are flags, go ahead and assign as needed. Ideally we'll have some kind of config file for the models
#if flag_list:
	# assign a root model first (give a model flag to my_tree)
	# assign models to each flag	
	




# Build model
print "building model"
model = misc.Model()

# state frequencies
fgen = stateFreqs.EqualFreqs(by='amino', alnfile='hrh1.fasta', save='stateFreqs.txt')
model.stateFreqs = fgen.getCodonFreqs()
fgen.save2file()

# parameters. Be sure you give it the correct parameters given the model you want to do
model.params={"kappa": 0.5, "omega": 2.5}

# matrix
mgen = matrixBuilder.GY94Matrix(model)
model.Q = mgen.buildQ()

print "evolving"
# Evolve along the tree
outfile="stephanie.fasta"
sequenceLength=100 # number of codons
myEvolver = Evolver(sequenceLength, model, my_tree, outfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeAlignment()


#######################################
### To incorporate:
# Vary Q over positions. We can define several Q's and randomly assign sites a rate matrix.
# Vary Q over branches. Claus seems to think temporal, but I don't know if that makes any sense. Yes we'd expect changes over time, but why should this be consistent across divergent lineages? That suggests an environmental effect

# Have a function to do something or other with lambda.s























