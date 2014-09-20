# Run the simulator pipeline. This will eventually be totally overhauled when config, user file, sanity checks, etc. introduced.

import numpy as np
import random as rn
import sys
sys.path.append("src/")
from misc import *
from newick import *
from state_freqs import *
from matrix_builder import *
from evolver import *


print "Reading tree"
my_tree,flags = read_tree(file="trees/10.tre", flags=True) 

# USE THESE LINES ONLY WHEN flag = True FOR readTree 
root_model_name = 'root_model'
flags.append(root_model_name)


# for now, shared by all models/partitions
print "frequencies"
freqObject = Random_Frequencies(by = 'codon')
myFrequencies = freqObject.calculate_freqs(savefile = 'freqs.txt')
mu = {'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}
kappa = 3.5
mu['AG'] = mu['AG'] * kappa
mu['CT'] = mu['CT'] * kappa

#### MODELS DEFINED ####
root_model = Model()
root_model.params = {'state_freqs': myFrequencies, 'mu': mu, 'alpha': 1.0, 'beta': 0.05}
m = codonGY_Matrix(root_model)
root_model.Q = m.buildQ()


m1 = Model()
m1.params = {'state_freqs': myFrequencies, 'mu': mu, 'alpha': 1.0, 'beta': 3.5}
m = codonGY_Matrix(m1)
m1.Q = m.buildQ()

m2 = Model()
m2.params = {'state_freqs': myFrequencies, 'mu': mu, 'alpha': 1.0, 'beta': 1.5}
m = codonGY_Matrix(m2)
m2.Q = m.buildQ()


partitions = []
numPart =  1
partLen = 50
for n in range(numPart):
    temp = {}
    for flag in flags: 
        temp[flag] = eval(flag) # requires that models are named same as flags corresponding to their introduction. this can probably be relaxed later.
    partitions.append( (partLen, temp ) )
# print partitions 
# [(10000, {'m1': <misc.Model instance at 0x10e4e78c0>, 'rootModel': <misc.Model instance at 0x10e4e6bd8>, 'm2': <misc.Model instance at 0x10e4e6b90>})]

print "Evolving"
myEvolver = Evolver(partitions, root_model_name) # first arg is partition list, second arg is the name for the rootModel
myEvolver.simulate(my_tree) # Since this function is recursive, provide tree here [not to Evolver constructor/whatever python calls it.] 
myEvolver.write_sequences(outfile = 'temp.fasta')







