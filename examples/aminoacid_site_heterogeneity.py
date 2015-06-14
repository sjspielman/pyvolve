# This example script demonstrates how to evolve according to an amino-acid model with sitewise rate heterogeneity. 

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a nucleotide model, as a pyvolve.Model object. For this example, we'll use default parameters, but see the example script custom_aminoacid.py for other options

# To implement rate heterogeneity, do either of these:
## 1) Custom rates: Provide a list of rate_factors when defining a Model object. These rate factors will be assigned to sites with equal probability by default. To change this, provide probabilities with the argument `rate_probs`.
## 2) Gamma rates: Provide the keyword arguments num_categories and alpha when defining a Model object. <num_categories> rates will be drawn from a gamma distribution with shape and scale parameter each equal to <alpha>. These rates will be equiprobable, unless overridden by `rate_probs`.

# Several model definitions are shown below (first argument can be a different model, as desired).

# custom rates
my_model1 = pyvolve.Model("WAG", rate_factors = [0.3, 0.8, 1.5, 2.45] ) # 25% of sites will have each factor.
my_model2 = pyvolve.Model("WAG", rate_factors = [0.3, 0.8, 1.5, 2.45], rate_probs = [0.7, 0.2, 0.05, 0.05] ) # 70% of sites evolve with rate of 0.3, 20% with a rate of 0.8, 5% with a rate of 1.5, and 5% with a rate of 2.45

# gamma rates
my_model3 = pyvolve.Model("WAG", alpha = 0.6, num_categories = 5) 


# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 positions
my_partition = pyvolve.Partition(models = my_model2, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()