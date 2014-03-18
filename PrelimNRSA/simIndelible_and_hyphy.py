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
	

	
############################################### PARAMETERS, 3/17/14.  ###########################################################
kappas = np.linspace(2.50, 6.50, num = 100)
omegas = np.linspace(0.05, 1.10, num = 100)
omega_fixed = 1.0
kappa_fixed = 2.5


# NOTE:Setting equal freqs for run on 3/18
freqaln = 'aln_aa_H1N1_HA.fasta'

hydir = 'HyPhyMaterials/'
#hyphy_exec = "/home/sjs3495/bin/bin/HYPHYMP"
hyphy_exec = "HYPHYMP"
#################################################################################################################################


############################################ CLUSTER-SPECIFIC FOR AN ARRAY JOB ##################################################
param = sys.argv[1]
rep = int(sys.argv[2])
if rep == 100:
	rep = 0

param_of_interest = 0
if param == "kappa":
	omega = omega_fixed
	kappa = kappas[rep]
	param_of_interest=kappa
elif param == "omega":
	kappa = kappa_fixed
	omega = omegas[rep]
	param_of_interest=omega
else:
	sys.exit("big fail! you have no param!")
assert(param_of_interest != 0), "Param not defined properly"


param_dir = "params/"
seq_dir = "indelibleOmega_Seqs/"
ensure_dir(param_dir)
ensure_dir(seq_dir)

tree_string = getTree("200rand.tre")

#name = param+str(rep)
#seqfile = seq_dir+name+".fasta"
#paramfile = param_dir+name+".txt"
#################################################################################################################################

####### TO SIMULATE ON MAC SINCE INDELIBLE WONT WORK ON CLUSTER #################################################################
#count=0
#for omega in omegas:
	#indelibleSim(kappa, omega)
	
	# Save sequences
	#out = "seq"+str(count)+".fasta"
	#shutil.copy('results.fas', out)
	#count+=1

#assert 1==0	
#################################################################################################################################

seqfile = seq_dir+"seq"+str(rep)+".fasta"
shutil.copy(seqfile, hydir+"seqs.fasta")
os.chdir(hydir)

# Call HyPhy
print "hyphying"
prepHyPhy('seqs.fasta', tree_string)
hyparam = runHyPhy(hyphy_exec, param)

# Save param results
file = open("../"+paramfile, 'w')
file.write(str(param_of_interest)+'\t'+str(hyparam))
file.close()




















