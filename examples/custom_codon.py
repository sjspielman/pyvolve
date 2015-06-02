# This example script demonstrates how to evolve according to a customized codon model. Customizable parameters include dN/dS, mutation rates, and equilibrium frequencies.

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a codon model, as a pyvolve.Model object. Below are several example dictionaries of customized parameters. Note that any unspecified parameters remain the default value. Further, for codon models, dN/dS is always required!

codon_freqs = [0.02792, 0.01502, 0.01755, 0.01635, 0.02512, 0.01505, 0.01997, 0.00779, 0.02043, 0.01193, 0.01404, 0.01176, 0.01449, 0.01774, 0.00658, 0.02969, 0.02937, 0.00726, 0.01316, 0.00458, 0.02227, 0.00045, 0.00697, 0.00368, 0.01169, 0.01274, 0.01866, 0.0125, 0.00914, 0.00119, 0.02332, 0.02301, 0.00315, 0.02554, 0.02328, 0.01468, 0.02868, 0.02669, 0.00417, 0.01947, 0.0145, 0.01586, 0.02783, 0.01179, 0.006, 6e-05, 0.00549, 0.02555, 0.03147, 0.03111, 0.02524, 0.00276, 0.02051, 0.01129, 0.02267, 0.02258, 0.00012, 0.03009, 0.02104, 0.02865, 0.0283]

# Customize all mutation rates, state frequencies, and provide dN/dS = 1.5
parameters1 = {"mu": {"AC": 1.5, "AG": 2.5, "AT": 0.5, "CG": 0.8, "CT": 0.99, "GT": 1.56}, "state_freqs": codon_freqs, "omega": 1.5}

# Customize mutation rates using only transtion-to-tranversion bias, state frequencies, and provide dN/dS = 0.75
parameters2 = {"kappa": 4.5, "state_freqs": codon_freqs, "omega":0.75}

# Customize mutation rates using only transtion-to-tranversion bias, retain default (equal) state frequencies, and provide separate dN (0.75) and dS (0.8)
parameters3 = {"kappa": 6.5, "beta": 0.75, "alpha": 0.8}

# Customize all mutation rates and retain default (equal) state frequencies. Use dN = 0.95, dS = 0.97
parameters4 = {"mu": {"AC": 1.5, "AG": 2.5, "AT": 0.5, "CG": 0.8, "CT": 0.99, "GT": 1.56}, "beta":0.95, "alpha":0.97}

# Customize state frequencies only and retain default (equal) mutation rates. Use dN/dS = 0.35
parameters5 = {"state_freqs": codon_freqs, "omega":0.35}

my_model = pyvolve.Model("nucleotide", parameters3) # Any of the dictionaries shown above is acceptable!

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 codon positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()