#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################
## 
## This script provides a basic outline of how to run pyvolve with site-heterogeneity for nucleotide and/or amino acid models.
## The example given here uses an empirical amino acid model.
## It is assumed that pyvolve has been properly installed (see README for installation instructions).
##############################################################################



# Import pyvolve
try:
    from pyvolve import *
except:
    try:
        import sys
        sys.path.append("../src/")
        from misc import *
        from newick import *
        from state_freqs import *
        from matrix_builder import *
        from evolver import *
    except:
        raise AssertionError("\nWhere's pyvolve!!")



##### Read in a newick phylogeny #####
my_tree = read_tree(file = "trees/basic.tre") # This file contains a 10-taxon phylogeny, *with branch lengths*!

##### Set up evolutionary model #####
freq = EqualFrequencies(by = 'amino')()
my_model = Model()
my_model.params = {'state_freqs': freq, 'aa_model': 'LG'} # Use the key 'aa_model' to select an empirical amino acid replacement matrix (JTT, WAG, or LG)
my_model.matrix = aminoAcid_Matrix(my_model.params)()
# Here's the heterogeneity part! Site-wise heterogeneity is modeled with discrete rate categories. You can either provide a list of rates and an associated list of probabilities, or you can draw from a gamma distribution. Let's draw!
rates, probs = draw_gamma_rates(0.45, 5) # Tell the function to draw *5* rate categories from a gamma distribution with a shape parameter (alpha) = 0.45. A final optional argument of a list of associated probabilities may also be provided. If not provided, random probabilities are drawn. 
my_model.rates = rates
my_model.probs = probs

##### Define a partition which uses this model #####
my_partition = Partition()
my_partition.size = 55
my_partition.models = my_model


##### Evolve! #####
my_evolver = Evolver(partitions = my_partition, tree = my_tree, seqfile = "site_heterogeneity_amino.fasta", ratefile = "my_rates.txt", infofile = "rate_info.txt") 
my_evolver() 
# Two files my_rates.txt and rate_info.txt will be output. (Note that analogous files under the names site_rates.txt and site_rate_info.txt, respectively, are these default file names!)
# my_rates.txt tells, for each site in the simulated alignment, which partition it belongs to and which rate category it belongs to.
# rate_info.txt gives the rate factors for each rate category, across partitions.






