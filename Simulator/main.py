import numpy as np
import random as rn
from Bio import AlignIO
import sys
import time
sys.path.append("src/")

from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

# get dna/amino/codon lists
molecules = Genetics()

# read in tree
print "reading tree"
my_tree, flag_list  = readTree(file="trees/10.tre", show=False) # set True to print out the tree


# state frequencies
print "collecting state frequencies"
fgen = ReadFreqs(by='amino', alnfile='hrh1.fasta', save='stateFreqs.txt')
commonFreqs = fgen.getCodonFreqs()
fgen.save2file()

## temporary code for constructing multiple GY94 models. Will formalize in the future
numPart = 15 # number of partitions
kappa  = 4.5
omegas = np.linspace(0.01, 0.5, num=numPart, dtype=float) # 15 rate categories, inclusive
partLen = np.tile(20, numPart) # for now, equal size partitions
partitions = []
print "constructing models for", numPart, "partitions"
for i in range(numPart):
	# Define model object and its parameters. Build model matrix. Add tuple (partition length, model) to partitions list
	model = misc.Model()
	model.params = { "kappa": kappa, "omega": omegas[i], "stateFreqs": commonFreqs }
	m = GY94(model)
	model.Q = m.buildQ()
	partitions.append( (partLen[i], model) )

# Evolve
outfile = time.strftime("%m.%d_%H;%M;%S")+".fasta" # So, : is an illegal filename character
myEvolver = StaticEvolver(partitions, my_tree, outfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeAlignment()










######################### STUFF ###################
# parameters for rodrigue model
#model1.params={"nucMut":   {"AC":0.16, "AG":0.16, "AT":0.16, "CG":0.16, "CT":0.16, "GT":0.17}, 
#			   "nucFreqs": {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
#		       "aaVector": [0.1, 0.05, 0.1, 0.05, 0.025, 0.025, 0.025, 0.025, 0.05, 0.4, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.05, 0.05]	
#			  }


















