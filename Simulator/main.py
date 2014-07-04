# Run the simulator pipeline. This will eventually be totally overhauled when config, user file, sanity checks, etc. introduced.

# IMPORTANT IMPORTANT!!!!  FOR NOW IT IS ALL HARD-CODED AS A CODON MODEL


import numpy as np
import random as rn
import sys
import time
sys.path.append("src/")

from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

genetics = Genetics()

print "Reading tree"
my_tree= readTree(file="trees/10.tre")      #, show=True) # set show=True to print out the tree
printTree(my_tree)
assert 1==6

freqObject = BoltzmannFreqs(by = 'amino', factor = 50, rank = ['D', 'A', 'C', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
myFrequencies = freqObject.calcFreqs(type = 'amino')#, savefile = 'bob.txt')
print myFrequencies
assert 1==5












###### PIPELINE ########
# Read in the tree
# Select a category of molecular evolution model (nucleotide, amino acid, codon, empirical codon, mutation-selection) 
# Set parameters for the model. Includes, for each partition, equilibrium frequencies and required model parameters. 
# Each partition needs a length, too
# Evolve

##########################################################################################
####### READ IN THE TREE #########
print "Reading tree"
my_tree= readTree(file="trees/10.tre")      #, show=True) # set show=True to print out the tree

##########################################################################################
############################## SELECT A MODEL OF EVOLUTION ###############################
##########################################################################################

evoModel = 'codon'
# All models consider equilibrium frequencies
# amino acid params -> empirical matrix
# nucleotide params -> mutation rates
# codon params      -> alpha, beta, mutation rates [note for omega, set alpha=1]
# mutsel params     -> mutation rates, amino acid propensities   
# LATER THIS WILL BECOME IMPORTANT FOR ESTABLISHING MODEL PARAMETERS. FOR NOW, ALL BELOW IS HARD-CODED FOR A CODON MODEL.


##########################################################################################
############################## ESTABLISH EVOLUTIONARY MODEL ##############################
##########################################################################################
print "Constructing evolutionary model(s)"
# Partition specification
partitions = []
numPart =  1
partLen = 100

# Equilbrium frequencies
myFreqs = {'I': 0.33, 'L':0.33, 'V':0.34}
freqObject = UserFreqs(type = 'codon', by = 'amino', freqs = myFreqs)
myFrequencies = freqObject.calcFreqs()
print myFrequencies


# Param dictionary, including the equilibrium frequencies
muCodonParams = {'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}
kappa = 1.0
muCodonParams['AG'] = muCodonParams['AG'] * kappa
muCodonParams['CT'] = muCodonParams['CT'] * kappa
codonParams = {'stateFreqs': myFrequencies, 'mu': muCodonParams, 'alpha': 0.5, 'beta': 0.25}

#muMutSelParams = {'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}
#mutSelParams = {'stateFreqs': myFrequencies, 'mu': muMutSelParams}




# Construct model
for i in range(numPart):
    model = misc.Model()
    #model.params = mutSelParams
    model.params = codonParams
    m = codon_MatrixBuilder(model)
    #m = mutSel_MatrixBuilder(model)
    model.Q = m.buildQ()
    partitions.append( (partLen, model) )
#np.savetxt('matrix.txt', model.Q)
##########################################################################################
###################################### EVOLVE ############################################
##########################################################################################
print "Evolving"
#outfile = time.strftime("%m.%d_%H;%M;%S")+".fasta" # So, : is an illegal filename character
myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = 'temp.fasta')
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeSequences()


# Write true rates
#truefile = "truerates.txt"
#truef = open(truefile, 'w')
#truef.write("position\tomega\n")
#position = 1
#for i in range(numPart):
#    for l in range(partLen):
#        truef.write(str(position)+'\t'+str(omegas[i])+'\n')
#        position+=1
#truef.close()










######################### STUFF ###################
# parameters for rodrigue model
#model1.params={"nucMut":   {"AC":0.16, "AG":0.16, "AT":0.16, "CG":0.16, "CT":0.16, "GT":0.17}, 
#               "nucFreqs": {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
#               "aaVector": [0.1, 0.05, 0.1, 0.05, 0.025, 0.025, 0.025, 0.025, 0.05, 0.4, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.05, 0.05]    
#              }


















