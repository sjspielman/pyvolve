#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    This example script shows how to simulate amino-acid sequences.
    Available amino-acid models include:
        1. JTT
        2. WAG
        3. LG
'''


from pyvolve import *


############## EXAMPLE USING DEFAULT PYVOLVE SETTINGS #############
## This example simulates a 50-column amino-acid alignment according to the LG model. There is no heterogeneity in the evolutionary process.


#### Read in tree from file ####
mytree = read_tree(file = "trees/basic.tre")

#### Set up evolutionary model ####
model = Model("LG") # Other options include "WAG" and "JTT"
model.construct_model()

#### Define partition ####
part = Partition(size = 50, models = model)

#### Evolve ####
Evolver(tree = mytree, partitions = part)() 

####################################################################


############## EXAMPLE USING CUSTOM SETTINGS #############
## This example simulates a 50-column alignment with gamma rate heterogeneity, again according to the LG model.

#### Read in tree from file ####
mytree = read_tree(file = "trees/basic.tre")

#### Set up evolutionary model ####
model = Model("LG")
model.construct_model(alpha = 0.5, num_categories = 5) # site-rate heterogeneity with 5 categories drawn from a discrete gamma distribution with shape = 0.5

#### Define partition ####
part = Partition(size = 50, models = model)

#### Evolve ####
Evolver(tree = mytree, partitions = part)() 



############## EXAMPLE USING CUSTOM SETTINGS #############
## This example simulates a 50-column alignment with custom rate heterogeneity, again according to the LG model.

#### Read in tree from file ####
mytree = read_tree(file = "trees/basic.tre")

#### Set up evolutionary model ####
model = Model("LG")
model.construct_model(rate_factors =[0.2, 1.2, 2.7], rate_probs = [0.25, 0.5, 0.25]) # Custom rate scalars with associated probabilities. (25% of sites will evolve with 0.2, 50% with 1.2, and finally 25% with 2.7). Alternatively, if you
# model.construct_model(rate_factors =[0.2, 1.2, 2.7]) # By not specifying `rate_probs`, rate_factors are given equal probabilities (in this case, 1/3 each) 

#### Define partition ####
part = Partition(size = 50, models = model)

#### Evolve ####
Evolver(tree = mytree, partitions = part)() 


############## EXAMPLE USING CUSTOM SETTINGS #############
## This example simulates a 50-column alignment without rate heterogeneity, but with custom state frequencies (note that default frequencies are equal).

#### Read in tree from file ####
mytree = read_tree(file = "trees/basic.tre")

#### Set up evolutionary model ####
frequencies = RandomFrequencies("amino")() # See the example file setup_frequencies.py for more ways to specify state frequencies
model = Model("LG", params = {'state_freqs':frequencies})
model.construct_model()

#### Define partition ####
part = Partition(size = 50, models = model)

#### Evolve ####
Evolver(tree = mytree, partitions = part)() 



