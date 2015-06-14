# This example script demonstrates how to evolve according to a codon model with site heterogeneity. 

import pyvolve
import numpy as np 

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a codon model, as a pyvolve.Model object. The first argument can be either:
## 1) "GY" or "codon" for the GY94-style (uses codon equilibrium frequencies in the matrix)
## 2) "MG" for the MG94-style (uses nucleotide equilibrium frequencies in the matrix)

# To implement dN/dS heterogeneity, provide a list/numpy array of values to the dN and/or dS keys. Note that if you specify both dN and dS heterogeneity, lists must be the same length; dN/dS values will be assigned 1:1 (NOT all combinations!)
# By default, each category is equally likely. To change this, provide probabilities with the keyword argument `rate_probs` when defining a Model object.
# Several possible dictionaries to use are shown below.

parameters1 = {"omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5. 
parameters2 = {"beta": [0.5, 0.4, 0.3], "alpha": [0.98, 1.05, 0.75] } # dN/dS values are 0.5/0.98, 0.4/1.05, and 0.3/0.75.
parameters3 = {"beta": [0.5, 0.4, 0.3], "alpha": [0.98, 1.05, 0.75] } 


# Define model with equal probabilities for all dN/dS categories
my_model1 = pyvolve.Model("GY", parameters3)

# Define model with custom probabilities for each dN/dS category
my_model2 = pyvolve.Model("GY", parameters3, rate_probs = [0.5, 0.2, 0.3])# dN/dS= 0.5/0.98 now occurs at 50% of sites, 0.4/1.05 at 20% of sites, and 0.3/0.75 at 30% of sites.


# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 codon positions
my_partition = pyvolve.Partition(models = my_model1, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()