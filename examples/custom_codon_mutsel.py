# This example script demonstrates how to evolve according to a customized mutation-selection codon model. Customizable parameters include mutation rates, and either equilibrium frequencies or fitness values.
# Note that, for MutSel models, mutation rates do not have to be symmetric, so you can provide different rates for A->C ("AC") and C->A ("CA")
import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")


# Below are three example customized parameter dictionaries. Note that each of these could have "fitness" rather than "state_freqs" as a key
codon_freqs = [0.02792, 0.01502, 0.01755, 0.01635, 0.02512, 0.01505, 0.01997, 0.00779, 0.02043, 0.01193, 0.01404, 0.01176, 0.01449, 0.01774, 0.00658, 0.02969, 0.02937, 0.00726, 0.01316, 0.00458, 0.02227, 0.00045, 0.00697, 0.00368, 0.01169, 0.01274, 0.01866, 0.0125, 0.00914, 0.00119, 0.02332, 0.02301, 0.00315, 0.02554, 0.02328, 0.01468, 0.02868, 0.02669, 0.00417, 0.01947, 0.0145, 0.01586, 0.02783, 0.01179, 0.006, 6e-05, 0.00549, 0.02555, 0.03147, 0.03111, 0.02524, 0.00276, 0.02051, 0.01129, 0.02267, 0.02258, 0.00012, 0.03009, 0.02104, 0.02865, 0.0283]
custom_mutation_sym = {"AC": 1.5, "AG": 2.5, "AT": 0.5, "CG": 0.8, "CT": 0.99, "GT": 1.56} # For MutSel models, if you provide only 1 pair for each mutation rate (e.g. only AC and not CA), then Pyvolve will make mutation rates symmetric
custom_mutation_asym = {"AC": 1.5, "CA": 0.8, "AG": 2.5, "GA": 1.2, "AT": 0.5, "TA": 1.1, "CG": 0.8, "GC": 0.9, "CT": 0.99, "TC": 2.3, "GT": 1.56, "TC": 2.56} 

# Customize mutation rates using symmetric mutation rates, and specify frequencies for the MutSel model
parameters1 = {"state_freqs": codon_freqs, "mu":custom_mutation_syn}

# Customize mutation rates using asymmetric mutation rates, and specify frequencies for the MutSel model
parameters2 = {"state_freqs": codon_freqs, "mu":custom_mutation_asyn}

# Customize mutation rates using kappa, and specify frequencies for the MutSel model
parameters3 = {"state_freqs": codon_freqs, "kappa":4.25}

my_model = pyvolve.Model("mutsel", parameters3) # Any of the dictionaries shown above is acceptable!

# Assign the model to a pyvolve.Partition. The size argument indicates to evolve 250 codon positions
my_partition = pyvolve.Partition(models = my_model, size = 250)

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver()