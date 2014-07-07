# Run the simulator pipeline. This will eventually be totally overhauled when config, user file, sanity checks, etc. introduced.

import numpy as np
import random as rn
import sys
sys.path.append("src/")
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


print "Reading tree"
my_tree,flags = readTree(file="trees/10.tre", flags=True) 

# USE THESE LINES ONLY WHEN flag = True FOR readTree 
rootModelName = 'rootModel'
flags.append(rootModelName)


# for now, shared by all models/partitions
freqObject = EqualFreqs(by = 'codon')
myFrequencies = freqObject.calcFreqs()
mu = {'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}
kappa = 3.5
mu['AG'] = mu['AG'] * kappa
mu['CT'] = mu['CT'] * kappa

#### MODELS DEFINED ####
rootModel = misc.Model()
rootModel.substParams = {'stateFreqs': myFrequencies, 'mu': mu, 'alpha': 1.0, 'beta': 0.05}
m = mechCodon_MatrixBuilder(rootModel)
rootModel.Q = m.buildQ()


m1 = misc.Model()
m1.substParams = {'stateFreqs': myFrequencies, 'mu': mu, 'alpha': 1.0, 'beta': 3.5}
m = mechCodon_MatrixBuilder(m1)
m1.Q = m.buildQ()

m2 = misc.Model()
m2.substParams = {'stateFreqs': myFrequencies, 'mu': mu, 'alpha': 1.0, 'beta': 1.5}
m = mechCodon_MatrixBuilder(m2)
m2.Q = m.buildQ()


partitions = []
numPart =  1
partLen = 10000
for n in range(numPart):
    temp = {}
    for flag in flags: 
        temp[flag] = eval(flag) # requires that models are named same as flags corresponding to their introduction. this can probably be relaxed later.
    partitions.append( (partLen, temp ) )
# print partitions 
# [(10000, {'m1': <misc.Model instance at 0x10e4e78c0>, 'rootModel': <misc.Model instance at 0x10e4e6bd8>, 'm2': <misc.Model instance at 0x10e4e6b90>})]

print "Evolving"
myEvolver = IndelEvolver(partitions, rootModelName) # first arg is partition list, second arg is the name for the rootModel
myEvolver.simulate(my_tree) # Since this function is recursive, provide tree here [not to Evolver constructor/whatever python calls it.] 
myEvolver.writeSequences(outfile = 'temp.fasta')







