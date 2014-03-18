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
my_tree= readTree(file="trees/100.tre") # set True to print out the tree


# state frequencies
print "collecting state frequencies"
#fgen = UserFreqs(by='amino', freqs = {'I':0.2, 'L':0.2, 'V':0.2, 'C':0.2, 'F':0.2} , save='stateFreqs.txt') # these are 5 hydrophobic.
fgen = EqualFreqs(by='codon')#, file = '
commonFreqs = fgen.getCodonFreqs()
fgen.save2file()


## temporary code for constructing multiple GY94 models. Will formalize in the future
numPart =  1 # number of partitions
kappa  = 2.5
omegas = [0.25]


partLen = 1000 #partLen = np.tile(100, numPart) # for now, equal size partitions
partitions = []
print "constructing models for", numPart, "partitions"
for i in range(numPart):
	# Define model object and its parameters. Build model matrix. Add tuple (partition length, model) to partitions list
	model = misc.Model()
	model.params = { "kappa": kappa, "omega": omegas[i], "stateFreqs": commonFreqs }
	m = GY94(model)
	model.Q = m.buildQ()
	partitions.append( (partLen, model) )

# Evolve
print "evolving"
outfile = time.strftime("%m.%d_%H;%M;%S")+".fasta" # So, : is an illegal filename character
myEvolver = StaticEvolver(partitions = partitions, tree = my_tree)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeSequences()


# Write true rates
truefile = "truerates.txt"
truef = open(truefile, 'w')
truef.write("position\tomega\n")
position = 1
for i in range(numPart):
	for l in range(partLen):
		truef.write(str(position)+'\t'+str(omegas[i])+'\n')
		position+=1
truef.close()










######################### STUFF ###################
# parameters for rodrigue model
#model1.params={"nucMut":   {"AC":0.16, "AG":0.16, "AT":0.16, "CG":0.16, "CT":0.16, "GT":0.17}, 
#			   "nucFreqs": {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
#		       "aaVector": [0.1, 0.05, 0.1, 0.05, 0.025, 0.025, 0.025, 0.025, 0.05, 0.4, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.05, 0.05]	
#			  }


















