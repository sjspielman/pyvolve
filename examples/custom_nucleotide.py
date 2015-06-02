# This example script demonstrates how to evolve according to a customized nucleotide model. Customizable parameters include mutation rates and equilibrium frequencies.

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a nucleotide model, as a pyvolve.Model object. Below are several example dictionaries of customized parameters. Note that any unspecified parameters remain the default value.

# Customize all mutation rates and state frequencies
parameters1 = {"mu": {"AC": 1.5, "AG": 2.5, "AT": 0.5, "CG": 0.8, "CT": 0.99, "GT": 1.56}, "state_freqs": [0.1, 0.4, 0.4, 0.1]}

# Customize mutation rates using only transtion-to-tranversion bias, and state frequencies
parameters2 = {"kappa": 4.5, "state_freqs": [0.2, 0.4, 0.3, 0.1]}

# Customize mutation rates using only transtion-to-tranversion bias, and retain default (equal) state frequencies
parameters3 = {"kappa": 6.5}

# Customize all mutation rates and retain default (equal) state frequencies
parameters4 = {"mu": {"AC": 1.5, "AG": 2.5, "AT": 0.5, "CG": 0.8, "CT": 0.99, "GT": 1.56}}

# Customize state frequencies only and retain default (equal) mutation rates
parameters5 = {"state_freqs": [0.2, 0.4, 0.3, 0.1]}

my_model = pyvolve.Model("nucleotide", parameters3) # Any of the dictionaries shown above is acceptable!

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()