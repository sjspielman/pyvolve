import numpy as np
import random as rn
from Bio import AlignIO
import sys
import time
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
model1 = misc.Model()
model2 = misc.Model()

# state frequencies
fgen = stateFreqs.EqualFreqs(by='amino', alnfile='hrh1.fasta', save='stateFreqs.txt')
commonFreqs = fgen.getCodonFreqs()
model1.stateFreqs = commonFreqs
model2.stateFreqs = commonFreqs
fgen.save2file()

# parameters
model1.params={"kappa": 3.5, "omega": 2.5}
model2.params={"kappa": 5.2, "omega":0.75}
# matrix
mgen1 = matrixBuilder.GY94Matrix(model1)
mgen2 = matrixBuilder.GY94Matrix(model2)
model1.Q = mgen1.buildQ()
model2.Q = mgen2.buildQ()

# List of partitions. Each partition is a tuple: (size (int, number of codons), model (a model object)). Can be length of 1, leading to NO site heterogeneity.
partitions = [ (100, model1), (100, model2) ] # We can shuffle sites later if we really feel like it.

# Evolve
outfile = time.strftime("%m.%d_%H;%M;%S")+".fasta" # So, : is an illegal filename character
myEvolver = Evolver(partitions, my_tree, outfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeAlignment()


#######################################
### To incorporate:
# Vary Q over branches. Claus seems to think temporal, but I don't know if that makes any sense. Yes we'd expect changes over time, but why should this be consistent across divergent lineages? That suggests an environmental effect

# Have a function to do something or other with lambda.s























