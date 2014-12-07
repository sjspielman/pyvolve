#! /usr/bin/env python

#from pyvolve import *
import sys
sys.path.append("/Users/sjspielman/Research/pyvolve/src/")
from misc import *
from models import *
from partition import *
from newick import *
from state_freqs import *
from matrix_builder import *
from evolver import *
    
    

my_tree = read_tree(file = 'examples/trees/basic.tre')

f = RandomFrequencies(by = 'codon')()

codonmodel_params = {'state_freqs':f, 'kappa':2.75, 'beta':[1, 2.5, 0.5], 'alpha':[1, .9, 1.1]}
cm = CodonModel(codonmodel_params, "codon")
cm.construct_model(rate_probs = [0.25, 0.25, 0.5])
print cm.params
model_params = {'state_freqs':f, 'kappa':4.5, 'beta':1.5, 'alpha':2}
m = Model(model_params, "MG94")
m.construct_model()
part1 = Partition()
part1.size = 50
part1.models = m
#print part1.models.params


Evolver(partitions = part1, tree = my_tree, seqfile = None)
