# This example script demonstrates how to evolve according to a customized mutation-selection nucleotide model. Customizable parameters include mutation rates, and either equilibrium frequencies or fitness values.
# Note that, for MutSel models, mutation rates do not have to be symmetric, so you can provide different rates for A->C ("AC") and C->A ("CA")
import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")


# Below are three example customized parameter dictionaries. Note that each of these could have "fitness" rather than "state_freqs" as a key
nuc_freqs = [0.334, 0.12, 0.41, 0.136]
custom_mutation_sym = {"AC": 1.5, "AG": 2.5, "AT": 0.5, "CG": 0.8, "CT": 0.99, "GT": 1.56} # For MutSel models, if you provide only 1 pair for each mutation rate (e.g. only AC and not CA), then Pyvolve will make mutation rates symmetric
custom_mutation_asym = {"AC": 1.5, "CA": 0.8, "AG": 2.5, "GA": 1.2, "AT": 0.5, "TA": 1.1, "CG": 0.8, "GC": 0.9, "CT": 0.99, "TC": 2.3, "GT": 1.56, "TC": 2.56} 

# Customize mutation rates using symmetric mutation rates, and specify frequencies for the MutSel model
parameters1 = {"state_freqs": nuc_freqs, "mu":custom_mutation_sym}

# Customize mutation rates using asymmetric mutation rates, and specify frequencies for the MutSel model
parameters2 = {"state_freqs": nuc_freqs, "mu":custom_mutation_asym}

# Customize mutation rates using kappa, and specify frequencies for the MutSel model
parameters3 = {"state_freqs": nuc_freqs, "kappa":4.25}

my_model = pyvolve.Model("mutsel", parameters3) # Any of the dictionaries shown above is acceptable!

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()