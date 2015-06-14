# This example script demonstrates how to evolve according to the empirical codon model (ECM) model. All model parameters (except dN/dS!) are default: equal mutation rates empirical model frequencies.

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define an ECM model, as a pyvolve.Model object. The first argument can be either:
## 1) "ECMrest", for the restricted version of this model. This version considers only *single nucleotide* instantaneous changes.
## 2) "ECMunrest", for the unrestricted version of this model. This version considers all *single, double, and triple nucleotide* instantaneous changes.
# Note that ECM uses a special type of kappa parameterization. This may be specified in Pyvolve but is not shown here. Further, dN/dS may be used, although it will NOT have the same interpretation as regular!

my_model = pyvolve.Model("ECMrest")

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions (for a codon alignment, this means 250 codons, i.e. 750 nucleotide sites)
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()