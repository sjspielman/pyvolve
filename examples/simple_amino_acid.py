# This example script demonstrates how to evolve according to a simple amin-acid model. Customizable model parameters are default: equal equilibrium frequencies.

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define an amino-acid model, as a pyvolve.Model object. The first argument should be either "JTT", "WAG", or "LG" (available empirical matrices in Pyvolve)
my_model = pyvolve.Model("LG")

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()