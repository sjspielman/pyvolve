# This example script demonstrates how to evolve according to custom model with custom code

import pyvolve
import numpy as np

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a custom model with custom matrix and custom code (states). The matrix must be square and have the same dimension (in 1D) as the provided code. Note that code is a list because, in theory, you can specify multi-character (as in letters) states.
matrix = np.array([ [-0.5, 0.25, 0.25], [0.25, -0.5, 0.25], [0.25, 0.25, -0.5] ]) 
code = ["0", "1", "2"]
my_model = pyvolve.Model("custom", {"matrix":matrix, "code":code} )

my_partition = pyvolve.Partition(models = my_model, size = 1)

my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()