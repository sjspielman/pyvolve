######## TO BE RUN ON MAC!!!! ###########

import os
import re
import sys
import shutil
import subprocess
import numpy as np


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
def indelibleSim(kappa, omega):
	cfile=open('control.txt', 'r')
	control=cfile.read()
	cfile.close()
	
	## line is \s+[submodel]\s+kappa\s+omega
	control=re.sub(r'\n\s+\[submodel\]\s+\d\.\d*\s+\d\.\d*', r'\n    [submodel]    '+str(kappa)+'    '+str(omega), control)
	
	cfile=open('control.txt', 'w')
	cfile.write(control)
	cfile.close()
		
	runIndelible = subprocess.call('indelible control.txt', shell=True)
	assert (runIndelible == 0), "Indelible did not run properly."
#################################################################################################################################

	
##################################################### PARAMETERS ################################################################
kappas = np.linspace(2.50, 6.50, num = 100)
omegas = np.linspace(0.05, 1.10, num = 100)
omega_fixed = 1.0
kappa_fixed = 2.5


param = sys.argv[1]

if param == "kappa":
	omega = omega_fixed
elif param == "omega":
	kappa = kappa_fixed
else:
	sys.exit("big fail! you have no param to vary!")


seq_dir = param+"seqs/"
ensure_dir(seq_dir)

tree_string = getTree("200rand.tre")

count=0
for omega in omegas:
	indelibleSim(kappa, omega)

	# Save sequences
	out = "seqs"+str(count)+".fasta"
	shutil.copy('results.fas', seq_dir+out)
	count+=1


















