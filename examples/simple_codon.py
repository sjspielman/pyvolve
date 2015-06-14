# This example script demonstrates how to evolve according to a simple codon model. All model parameters (except dN/dS!) are default: equal mutation rates and equal equilibrium frequencies.

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a codon model, as a pyvolve.Model object. The first argument can be either:
## 1) "GY" or "codon" for the GY94-style (uses codon equilibrium frequencies in the matrix)
## 2) "MG" for the MG94-style (uses nucleotide equilibrium frequencies in the matrix)

# Codon models require you to specify a second argument to pyvolve.Model, a dictionary of parameters. You must specify dN/dS using either "omega" (for the full ratio), or "beta" for dN and "alpha" for dS, as shown below. Either dictionary would be acceptable.
parameters_omega = {"omega": 0.65}
parameters_alpha_beta = {"beta": 0.65, "alpha": 0.98} # Corresponds to dN/dS = 0.65 / 0.98
my_model = pyvolve.Model("MG", parameters_alpha_beta)

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions (for a codon alignment, this means 250 codons, i.e. 750 nucleotide sites)
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()