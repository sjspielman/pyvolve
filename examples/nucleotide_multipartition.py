# This example script demonstrates how to evolve according to a nucleotide model with several partitions.
# In this example, the first partition has gamma-distributedsitewise rate heterogeneity, the second partition is homogenous, and the third partition has custom sitewise rate heterogeneity.
# All models use default mutation-rate parameters

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")


# Define first model and partition. This partition has a length of 50 positions
model1 = pyvolve.Model("nucleotide", alpha = 0.7, num_categories = 4 )
part1 = pyvolve.Partition(models = model1, size = 50)

# Define second model and partition. This partition has a length of 20 positions
model2 = pyvolve.Model("nucleotide")
part2 = pyvolve.Partition(models = model2, size = 20)

# Define second model and partition. This partition has a length of 100 positions
model3 = pyvolve.Model("nucleotide", rate_factors = [0.5, 1.6, 4.1], rate_probs = [0.75, 0.2, 0.05])
part3 = pyvolve.Partition(models = model3, size = 100)




# Provide all partitions *in the order in which they should be evolved* to Evolver and evolve
my_evolver = pyvolve.Evolver(partitions = [part1, part2, part3], tree = my_tree)
my_evolver()