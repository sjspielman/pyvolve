# Run the simulator pipeline. This will eventually be totally overhauled when config, user file, sanity checks, etc. introduced.

import numpy as np
import random as rn
import sys
sys.path.append("src/")
from misc import *
from newick import *
from state_freqs import *
from matrix_builder import *
#from evolver import *


print "Reading tree"
my_tree,flags = read_tree(file="trees/10.tre", flags=True) 

# USE THESE LINES ONLY WHEN flag = True FOR readTree 
root_model = 'root_model'
flags.append(root_model)


# for now, shared by all models/partitions
freqObject = EqualFreqs(by = 'codon')
f = freqObject.calculate_freqs()

mu = {'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}
kappa = 3.5
mu['AG'] = mu['AG'] * kappa
mu['CT'] = mu['CT'] * kappa

#### MODELS DEFINED ####
root_model = Model()
root_model.subst_params = {'stateFreqs': f, 'mu': mu, 'alpha': 1.0, 'beta': 0.05}
m = mechCodon_MatrixBuilder(root_model)
root_model.Q = m.buildQ()


m1 = Model()
m1.subst_params = {'stateFreqs': f, 'mu': mu, 'alpha': 1.0, 'beta': 3.5}
m = mechCodon_MatrixBuilder(m1)
m1.Q = m.buildQ()

m2 = Model()
m2.subst_params = {'stateFreqs': f, 'mu': mu, 'alpha': 1.0, 'beta': 1.5}
m = mechCodon_MatrixBuilder(m2)
m2.Q = m.buildQ()


partitions = []
num_partitions =  1
partition_length = 100
for n in range(num_partitions):
    temp = {}
    for flag in flags: 
        temp[flag] = eval(flag) # requires that models are named same as flags corresponding to their introduction. this can probably be relaxed later.
    partitions.append( (partition_length, temp ) )
# print partitions 
# [(10000, {'m1': <misc.Model instance at 0x10e4e78c0>, 'root_model': <misc.Model instance at 0x10e4e6bd8>, 'm2': <misc.Model instance at 0x10e4e6b90>})]
print partitions
assert 1==5
print "Evolving"
myEvolver = Evolver(partitions, root_model) # first arg is partition list, second arg is the name for the root_model
myEvolver.simulate(my_tree) # Since this function is recursive, provide tree here [not to Evolver constructor/whatever python calls it.] 
myEvolver.writeSequences(outfile = 'temp.fasta')







