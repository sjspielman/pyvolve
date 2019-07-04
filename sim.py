from pyvolve import *

#(t1:0.0121481549,t2:0.0161978254,t3:0.0080811336);

t = read_tree(tree = "(t1:0.05, t2:0.07, t3:0.020);")
m = Model("MG", {"omega":1.0})
p = Partition(size=100000, models = m)
e = Evolver(partitions = p, tree = t)
e(algorithm = 1, ratefile = False, infofile = False)
print(e.branch_substitutions_codon)
print(e.branch_substitutions_nuc)
print(e.branch_substitutions_syn)
print(e.branch_substitutions_nonsyn)




######## now smaller
# alg 0 - iqtree good
# {'t1': 4839, 't2': 6779, 't3': 2001}
# {'t1': 4907, 't2': 6953, 't3': 2018}
# {'t1': 1269, 't2': 1647, 't3': 485}
# {'t1': 3570, 't2': 5132, 't3': 1516}

# alg 1
# {'t1': 5007, 't2': 6981, 't3': 2012}
# {'t1': 5093, 't2': 7145, 't3': 2027}
# {'t1': 1212, 't2': 1780, 't3': 524}
# {'t1': 3795, 't2': 5201, 't3': 1488}



###### confirmed ok for large branch lengths ##########
# iqtree got it totally fucking spot on
# {'t1': 48951, 't2': 68141, 't3': 19965}
# {'t1': 56317, 't2': 82091, 't3': 21265}
# {'t1': 10866, 't2': 14142, 't3': 4760}
# {'t1': 38085, 't2': 53999, 't3': 15205}


# alg 0
# {'t1': 38529, 't2': 48731, 't3': 17870}
# {'t1': 44895, 't2': 59913, 't3': 19078}
# {'t1': 8435, 't2': 10121, 't3': 4221}
# {'t1': 30094, 't2': 38610, 't3': 13649}