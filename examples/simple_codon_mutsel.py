# This example script demonstrates how to evolve according to a simple codon mutation-selection (MutSel) model. 
# For a MutSel model, you must supply either fitness values or equilibrium frequencies. Mutation rates are set as default (equal).

import pyvolve
import numpy as np  # imported to generate example mutation-selection model parameters

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a mutation-selection model, specifying a first argument of "MutSel". These models that you specify either a list of *fitness* values or a list of *equilibrium frequencies* in a parameters dictionary. See the user manual for more information here! Below are examples of acceptable dictionaries
parameters_fitness1 = {"fitness": np.random.uniform(-5, 5, size = 61)} # Numpy array of length 61 defines codon fitness
parameters_fitness2 = {"fitness": np.random.uniform(-5, 5, size = 20)} # Numpy array of length 20 defines amino-acid fitness, which are applied to codons such that synonymous codons will have the same fitness 
parameters_freqs    = {"state_freqs": np.repeat(1./61, 61)} # Numpy array of equal frequencies, just as an example! This list must sum to 1!

my_model = pyvolve.Model("MutSel", parameters_fitness1) # Any of the above parameters dictionaries are acceptable as the second argument!

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 codon positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()