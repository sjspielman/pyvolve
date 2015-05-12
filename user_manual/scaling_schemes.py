# SJS
# Script to simulate with dN/dS model to demonstrate scaling issues with Yang vs. neutral approach
# The R script plot_scaling_results.R will perform post-processing and plot the results

from pyvolve import *
from dendropy import *
import numpy as np
import subprocess

tree_string = "((((t11:0.006869548629,((t14:0.02744615343,(((((t17:0.009474494869,t13:0.02755469065):0.03020011152,t23:0.01637109739):0.06635743803,t18:0.05894008901):0.03145557342,t12:0.05107510961):0.003171503637,t20:0.04171180078):0.09332439501):0.006789946905,((t9:0.06393237635,(t7:0.04177429336,t25:0.07494378996):0.09668286787):0.06210535988,(t4:0.08427849947,(t8:0.01699725552,(t2:0.09613203502,t21:0.08250621709):0.0515646199):0.08814909002):0.001062202132):0.01712481421):0.07721304006):0.08939141796,t6:0.03179403555):0.06586158879,((t19:0.01458299537,t22:0.08076849533):0.04711015159,((t5:0.03288242628,(t1:0.09949872966,t24:0.007778408006):0.08514633912):0.006163714104,((t16:0.03605371944,t10:0.09133298337):0.09026077781,t3:0.0303836077):0.08432583029):0.03588969493):0.06904304423):0.09048717315,t15:0.002226415281);"
simtree = read_tree(tree = tree_string)
omegas =  np.arange(0.1, 2.1, 0.1) 

outstring = "scaling,rep,omega,treelen\n"

for dn in omegas:
    yang_seqfile = "yang.fasta"
    neutral_seqfile = "neutral.fasta"
    
    yang_model = Model("codon", {'omega': dn})
    yang_model.construct_model()
    neutral_model = Model("codon", {'omega': dn}, "neutral")
    neutral_model.construct_model()
    
    yang_partition = Partition(models = yang_model, size = 200)
    neutral_partition = Partition(models = neutral_model, size = 200)
    
    
    for rep in range(50):
        Evolver(partitions = yang_partition, tree = simtree, infofile = None, ratefile = None, seqfile = yang_seqfile)()
        Evolver(partitions = neutral_partition, tree = simtree, infofile = None, ratefile = None, seqfile = neutral_seqfile)()
    
        build_yang_tree = subprocess.call("raxmlHPC-AVX -m GTRGAMMA -p 12345 -s " + yang_seqfile + " -n yang ", shell = True)
        build_neutral_tree = subprocess.call("raxmlHPC-AVX -m GTRGAMMA -p 12345 -s " + neutral_seqfile + " -n neutral ", shell = True)
    
        yang_tree = Tree.get_from_path("RAxML_bestTree.yang", 'newick')
        yang_length = float(yang_tree.length())
        outstring += "yang," + str(rep) + "," + str(dn) + "," + str(yang_length) + "\n"
    
        neutral_tree = Tree.get_from_path("RAxML_bestTree.neutral", 'newick')
        neutral_length = float(neutral_tree.length())
        outstring += "neutral," + str(rep) + "," + str(dn) + "," + str(neutral_length) + "\n"
    
        subprocess.call("rm RAxML*", shell=True)
        subprocess.call("rm *.fasta", shell=True)
        subprocess.call("rm *.reduced", shell=True)
        
with open("compare_scaling_results.csv", "w") as outf:
    outf.write(outstring)
    
    

    

    
    
    