# This example script demonstrates how to evolve according to a nucleotide model with *branch* rate heterogeneity. The approach is the same for non-nucleotide models.

import pyvolve

# Define a phylogeny. For clarity, we define this tree with a string. The tree contains model flags for branches which should evolve according to new models. Flags are represented as _name_, where underscores surround the name.
my_tree = pyvolve.read_tree(tree = "((t1:0.5, t2:0.5):0.5_m1_,(t3:0.5, t4:0.5):0.5_m2_));")

# Define a model for each flag. Models should be given names with the keyword argument `name`. These names *MUST* have correspondingly named flags in the tree!
model1 = pyvolve.Model("nucleotide", {"kappa":3.5}, name = "m1")
model2 = pyvolve.Model("nucleotide", {"kappa":4.75}, name = "m2")
rootmodel = pyvolve.Model("nucleotide", name = "root") # We can also define, if we want, a model for the ROOT of the tree that is separate from either of these models.

# Define partition will all models as a list. Include the argument `root_model_name` to indicate the NAME ATTRIBUTE of the model that should be used at the root of the tree. This name's corresponding object must be in the `models` list. Note that a separate root model is not needed - you could easily just start with _m1_ at the root, but you'd still need to give "m1" to `root_model_name`.
my_partition = pyvolve.Partition(models = [model1, model2, rootmodel], size = 250, root_model_name = "root") 

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()