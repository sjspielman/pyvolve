#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################
## 
## This script provides a basic outline of how to run pyvolve. (Note that documentation, examples, etc. will be greatly expanded in the coming months. This is a temporary example file!!) 
## It is assumed that pyvolve has been properly installed (see README for installation instructions).
## To run pyvolve, you must:
##     1. Specify a phylogeny.
##     2. Specify an evolutionary model, including a) relevant model parameters, and  b) equilibrium frequencies
##     3. Specify evolutionary partitions (site-heterogeneity).
##     4. Finally, evolve!
##     
##
## In this example script, we will evolve two partitions. Both partitions will evolve according to the GY94 model (Goldman and Yang 1994), with few added flairs!
## The first partition will have kappa = 2.5 and omega (dN/dS) = 0.5. These parameters correspond to the standard parameters used in GY94.
## We will set up the second partition a little differently! We will use a GTR mutation model (rather than HKY85) and we will specify different rates for dN and dS. 
## We will provide custom equilibrium frequencies for the first partition, but equal frequencies for the second partition.
## Note that this example DOES NOT incorporate branch heterogeneity (but such an example is coming!)
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

# Read in a newick phylogeny using the function "read_tree". You may specify either a file containing the tree (flag "file"), or the newick string itself (flag "tstring").
my_tree = read_tree(file = "trees/basic.tre") # This file contains a 10-taxon phylogeny, *with branch lengths*!

############# Set up the evolutionary models. ###################

# FIRST, we will define equilibrium frequencies. pyvolve has extensive flexibility in defining equilibrium frequencies.
# In particular, the frequency *calculations* may use a different alphabet from the frequencies you ultimately want. For example, let's say you want to evolve amino acids, but you have information about the distribution of codons. pyvolve can compute codon steady state frequencies, given certain options, but then ultimately return amino acid frequencies for your evolutionary model.
# There are several python classes for defining equilibrium frequencies (see their docstrings for more info), including - 
## 1. EqualFrequencies           - Sets frequencies as equal (i.e. 1/4 for all nucleotides if by='nuc', and so on.) 
## 2. RandomFrequencies          - Computes (semi-)random frequency values for a given alphabet. Basically flat distributions but with some amount of noise.
## 3. CustomFrequencies          - Computes frequencies based on a user-provided dictionary of frequencies.
## 4. ReadFrequencies            - Computes frequencies from a sequence file. Contains an option to select specific columns from sequence file only, but this requires that the file is an alignemnt.
## 5. EmpiricalModelFrequencies  - Suitable for empirical models (AA models JTT, WAG and LG, and codon models ECMrest/ECMunrest). Will use the default state frequencies of the empirical model. 

# For the first partition, we will use the CustomFrequencies class. We will provide a dictionary of amino acid frequencies we are interested in, and then obtain the corresponding codon frequencies.
# The next line sets up an instance of StateFrequencies (but performs no calculations!). The "by" argument (can be amino/codon/nuc) tells pyvolve how we want to perform the calculations. The "freq_dict" argument, specific to CustomFrequencies, takes a *python dictionary* of frequencies.
freq_calculator_model1 = CustomFrequencies(by = 'amino', freq_dict = {'A':0.25, 'G':0.25, 'W':0.5})
# Compute the frequencies. We use the argument "type" to indicate that we want *codon* frequencies (not amino acid, as originally provided with "by"!)
frequencies_model1  = freq_calculator_model1(type = 'codon') 

# For the second partition, we will use the EqualFrequencies class to simply obtain equal codon frequencies. Nothing fancy! We will save these frequencies to a file with the argument "savefile", in the line that calls the calculate_freqs() function.
frequencies_model2 = EqualFrequencies(by = 'codon')(savefile = "partition2_codon_frequencies.txt") # all in 1 line! Note that if "type" is same as "by", no need to provide.

# SECOND, we will define our models for each partition. For this, we'll create a Model() object which will contain attributes representing model parameters and the instantaneous rate matrix.
# Model parameters are stored in a python dictionaries, with specific keys as parameter names (as shown below!). 
# The model instantaneous rate matrix is created using the MatrixBuilder children classes, which consist of - 
## 1. aminoAcid_Matrix         - Empirical amino acid models (JTT, WAG, LG)
## 2. nucleotide_Matrix        - GTR and all nested nucleotide models (incl. JC, K2P, TN93, etc.)
## 3. ECM_Matrix               - Empirical codon model by Kosiol 2007
## 4. mechCodon_Matrix         - So-called mechanistic codon matrices (GY94 and MG94)
## 5. mutSel_Matrix            - mutation-selection balance model proposed by Halpern and Bruno 1998. In pyvolve, mutSel models can also be used for DNA sequence evolution.
# Note that model parameters for the GY94 matrix include mutation rates, dN and/or dS, and state frequencies. 
model1 = Model()
model1.params = {'omega': 0.5, 'kappa': 2.5, 'state_freqs': frequencies_model1} # Define GY94 model parameters as attributes of the model1 object
build_matrix1 = mechCodon_Matrix(model1.params, type = "MG") # Define MatrixBuilder object. Any matrix builder takes a parameter dictionary as its argumentg
model1.matrix  = build_matrix1() # Build the matrix by calling this class and assign as attribute to model1 object


# Parameters for the second partition are a little different. To define different rates for dN and dS, we use the keys 'beta' and 'alpha', respectively, following the notation of MuseGaut1994.
# In addition, we can define arbitrary (symmetric!!) mutational parameters with the key 'mu', which itself is a dictionary.
model2 = Model()
part_mu = {'AC': 1.5, 'AG': 2.6, 'AT': 0.4, 'CG': 1.0, 'CT': 0.004, 'GT': 1.34} # Here, the rate of A->C and C->A is therefore 1.5, and so on. The keys for this list must be alphabetically ordered, as in "AG" should be used and not "GA".  
model2.params = {'beta': 2.5, 'alpha': 0.75, 'mu': part_mu, 'state_freqs': frequencies_model2}
model2.matrix = mechCodon_Matrix(model2.params)() # Define MatrixBuilder object and build matrix in 1 line!


# THIRD, we will define our partitions which make use of the models defined above. Each partition is again a Partition() object with several attributes.
part1 = Partition()
part1.size = 850     # partition1 should contain 850 positions (in this case, 850 codon positions = 2550 nucleotide positions)
part1.models = model1 
part2 = Partition()
part2.size = 1000
part2.models = model2





# FOURTH, we will evolve sequences according to our partitions along our provided phylogeny!
# First, set up an instance of the Evolver class, which does the evolving. We'll initiate an Evolver instance with several arguments, as follows...
part_list = [part1, part2]
my_evolver = Evolver(partitions = part_list, tree = my_tree, seqfile = "my_simulated_alignment.phy", seqfmt = "phylip")
# Now we actually call the Evolver instance. Your sequences will appear in phylip format in the file my_simulated_alignment.phy.
my_evolver() 





























