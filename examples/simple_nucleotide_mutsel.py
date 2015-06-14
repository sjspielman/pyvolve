# This example script demonstrates how to evolve according to a simple nucleotide mutation-selection (MutSel) model. 
# For a MutSel model, you must supply either fitness values or equilibrium frequencies. Mutation rates are set as default (equal).

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a mutation-selection model, specifying a first argument of "MutSel". These models that you specify either a list of *fitness* values or a list of *equilibrium frequencies* in a parameters dictionary. See the user manual for more information here! Below are examples of acceptable dictionaries
parameters_fitness = {"fitness": [1.5, 0.3, -0.5, 1.8]} # List of length 4 defines nucleotide fitness
parameters_freqs    = {"state_freqs": [0.3, 0.2, 0.15, 0.35]} # List of length 4 defines nucleotide frequencies. This list must sum to 1!

my_model = pyvolve.Model("MutSel", parameters_freqs) # Either of the above parameters dictionaries is acceptable as the second argument!

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions
my_partition = pyvolve.Partition(models = my_model, size = 250)


# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()