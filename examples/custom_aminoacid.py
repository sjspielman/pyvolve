# This example script demonstrates how to evolve according to a customized amino-acid model. The only customizable parameter here is equilibrium frequencies.

import pyvolve
import numpy as np # imported only to generate example amino-acid equilibrium frequencies

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define an amino-acid model, as a pyvolve.Model object. The first argument should be either "JTT", "WAG", "LG", "AB", "mtmam", "mtREV24", "DAYHOFF" (available empirical matrices in Pyvolve) Provide a parameters dictionary with key "state_freqs" to customize model. The provided frequency list must sum to 1! 
aa_freqs = [0.06364529, 0.07224862, 0.06057557, 0.0487006 , 0.05994598, 0.01482058, 0.07141465, 0.03867427, 0.04562237, 0.02361059, 0.03897318, 0.04462156, 0.03720072, 0.06157329, 0.0150518, 0.03246979, 0.0617337 , 0.08068358, 0.07970379, 0.04873008]
my_model = pyvolve.Model("LG", {"state_freqs": aa_freqs})


# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()