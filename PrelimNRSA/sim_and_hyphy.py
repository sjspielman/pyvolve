# SJS. Started 3/15/14.

# 3/17/14: rewritten for cluster. 
# Conduct 100 simulations each for varying kappa and omega. For varying omega, kappa constant, and vice versa.
# Tree has 200 taxa. Same state frequencies for all. 


import os
import re
import sys
import shutil
import subprocess
import numpy as np


# Use this line on my computer only, not on the cluster!
#sys.path.append("/Users/sjspielman/Research/MutSel/Simulator/src/")

from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


#################################################################################################################################
def ensure_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
#################################################################################################################################


#################################################################################################################################
def getTree(tree_file):
	tree = open(tree_file, 'rU')
	tree_string = tree.read()
	tree.close()
	tree_string = tree_string.rstrip()
	return tree_string
#################################################################################################################################


#################################################################################################################################
def buildSingleModel(freqs, kappa, omega, size):
	partitions = []
	model = misc.Model()
	model.params = { "kappa": kappa, "omega": omega, "stateFreqs": freqs }
	m = GY94(model)
	model.Q = m.buildQ()
	partitions.append( (size, model) )
	return partitions
#################################################################################################################################


#################################################################################################################################
def callSim(partitions, my_tree, seqfile):

	print "evolving"
	myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
	myEvolver.sim_sub_tree(my_tree)
	myEvolver.writeSequences()
#################################################################################################################################


#################################################################################################################################
def prepHyPhy(alnfile, tree):
	aln = open(alnfile, 'r')
	alnlines = aln.read()
	aln.close()
	hyin = open('temp.fasta', 'w')
	hyin.write(alnlines+'\n'+tree)
	hyin.close()	
#################################################################################################################################


#################################################################################################################################
def runHyPhy(hyphy_exec, param):

	if param == 'omega':
		param = 'w'
	elif param == 'kappa':
		param = 'k'

	outfile = "hyout.txt"
	callHyphy = hyphy_exec+" run.bf > "+outfile
	print "running hyphy as called:",callHyphy
	check = subprocess.call(callHyphy, shell=True)
	assert (check==0), "HyPhy didn't run properly."
	
	# Parse hyphy output
	file = open(outfile, 'rU')
	hyout=file.readlines()
	file.close()
	for line in hyout:
		find = re.search("^"+param+"=(\d+\.*\d*)$", line)
		if find:
			final_param = float(find.group(1))
			break
	return final_param
#################################################################################################################################
	
	
	
	
############################################### PARAMETERS, 3/18/14.  ###########################################################
kappas = np.linspace(2.50, 6.50, num = 100)
omegas = np.linspace(0.05, 0.9, num = 100)
omega_fixed = 1.0
kappa_fixed = 4.0



tree_string = getTree("100tree_samebl.tre")


size = 1000   #1000 codons per alignment

# NOTE:Setting equal freqs for run on 3/18
freqaln = 'aln_aa_H1N1_HA.fasta'

hydir = 'HyPhyMaterials/'
#hyphy_exec = "/home/sjs3495/bin/bin/HYPHYMP"
hyphy_exec = "HYPHYMP"
#################################################################################################################################


############################################ CLUSTER-SPECIFIC FOR AN ARRAY JOB ##################################################
param = sys.argv[1]
rep = int(sys.argv[2]) - 1

param_of_interest = 0
if param == "kappa":
	kappa = kappas[rep]
	omega = omega_fixed
	param_of_interest = kappa
elif param == "omega":
	kappa = kappa_fixed
	omega = omegas[rep]
	param_of_interest = omega
else:
	sys.exit("big fail! you have no param!")
assert(param_of_interest != 0), "Param not defined properly"



param_dir = "params/"
seq_dir = "seqs/"
ensure_dir(param_dir)
ensure_dir(seq_dir)

name = param+str(rep)
seqfile = seq_dir+name+".fasta"
paramfile = param_dir+name+".txt"
#################################################################################################################################

# Convert tree to sim-usable format
my_tree = readTree(tree = tree_string)

# Collect state frequencies
fgen = EqualFreqs(by='codon') #, alnfile=freqaln) ## by needs to be the type of data in the file
freqs = fgen.getCodonFreqs()

# Build the model
print "building model"
partition = buildSingleModel(freqs, kappa, omega, size)
	
# Simulate
print "simulating"
callSim(partition, my_tree, "seq.fasta")

# Save sequences
shutil.copy("seq.fasta", seqfile)

# Copy sequences into hydir
shutil.copy("seq.fasta", hydir)
	
os.chdir(hydir)

# Call HyPhy
print "hyphying"
prepHyPhy("seq.fasta", tree_string)
hyparam = runHyPhy(hyphy_exec, param)

# Save param results
file = open("../"+paramfile, 'w')
file.write(str(param_of_interest)+'\t'+str(hyparam))
file.close()



















