# SJS. Started 3/15/14.

# Conduct 100 simulations, 100 taxa each. Same state frequencies and kappa for all. Vary omega evenly (0.01-3). 
# Infer omega with hyphy, save resulting true sim and inferred. Plot to demonstrate that our simulation software works.


import os
import re
import sys
import subprocess
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *

sys.path.append("/Users/sjspielman/Research/MutSel/Simulator/src/")

from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


#################################################################################################################################
def ensure_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
    return 0
#################################################################################################################################

def prepSim(freq_aln, kappa, omegas, numPart, partLen):
	''' Use when frequencies are global and not site-specific. All sites share same freqs. ''' 
	partitions = []
	
	# global param
	fgen = ReadFreqs(by='amino', alnfile=freq_aln) ## by needs to be the type of data in the file
	
	## UNIFORM DISTRIBUTION
	#fgen = EqualFreqs(by='amino')
	
	freqs = fgen.getCodonFreqs()
	
	print "constructing models for", numPart, "partitions"
	for i in range(numPart):
		# Define model object and its parameters. Build model matrix. Add tuple (partition length, model) to partitions list
		model = misc.Model()
		model.params = { "kappa": kappa, "omega": omegas[i], "stateFreqs": freqs }
		m = GY94(model)
		model.Q = m.buildQ()
		partitions.append( (partLen, model) )

	return partitions

def buildSingleModel(freqs, kappa, omega, size):
	partitions = []
	model = misc.Model()
	model.params = { "kappa": kappa, "omega": omega, "stateFreqs": freqs }
	m = GY94(model)
	model.Q = m.buildQ()
	partitions.append( (size, model) )
	return partitions
	
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
def runHyPhy():
	outfile = "hyout.txt"
	callHyphy = "HYPHYMP run.bf > "+outfile
	print "running hyphy as called:",callHyphy
	check = subprocess.call(callHyphy, shell=True)
	assert (check==0), "HyPhy didn't run properly."
	
	# Parse hyphy output
	file = open(outfile, 'rU')
	hyout=file.readlines()
	file.close()
	for line in hyout:
		find = re.search("^w=(\d+\.\d*)$", line)
		if find:
			omega = float(find.group(1))
			break
	return omega
#################################################################################################################################
	
	
results_dir='/Users/sjspielman/Dropbox/MutSelProject/PreliminaryData/SimSeqs/'
ensure_dir(results_dir)


zero=1e-10

tree_string = "(t26:0.8157528578,((((t19:0.130143242,t80:0.7591279587):0.9975617055,t86:0.2094392162):0.7539197295,((t99:0.714059568,t38:0.3790428748):0.9462954316,((t8:0.02399082482,t69:0.4831460086):0.5862778549,(((t9:0.2352395647,t64:0.7096769114):0.8309416587,t50:0.2220890515):0.7175409112,t92:0.7234821764):0.4819605483):0.2681629495):0.02097934368):0.2734834186,((((t83:0.7843712003,t52:0.2901142205):0.1374024423,((t20:0.02574405214,((t74:0.3084802043,(t47:0.1138049611,t14:0.6317933276):0.3454655285):0.96972045,t24:0.3594481982):0.5955672206):0.7026010049,((t27:0.8782107995,t58:0.8896448328):0.7314315287,t34:0.1242245506):0.5693691911):0.976664278):0.4722047253,((t68:0.3034443925,t42:0.09149145428):0.05196091859,((t96:0.3599379426,t95:0.8641043194):0.9317530163,(t66:0.07659937069,((t7:0.5128562001,t55:0.06975672254):0.3449339604,t16:0.6597705341):0.7747370843):0.7912949338):0.8788079387):0.6016570281):0.2477935245,(((((t72:0.3799303034,t75:0.3532710788):0.07015336538,((t36:0.02703211713,t90:0.4426369618):0.3160473928,((t98:0.8061463048,t89:0.144463154):0.1857714909,(((t62:0.5307802558,t94:0.8982971627):0.7463883571,t63:0.7711612736):0.4759638396,(t53:0.7423221956,t17:0.686133717):0.132485118):0.9145437798):0.3240567176):0.7411617341):0.1697024996,(t59:0.8547446851,t76:0.3970330856):0.6356236765):0.3735436248,t78:0.3690058561):0.2423712793,((((t4:0.9409420351,t25:0.1384483045):0.6021607046,(t41:0.5857174119,t77:0.830866555):0.3623668503):0.3275486603,(((((t87:0.5843094843,(t29:0.9604068338,(t35:0.6937076449,t30:0.5904012059):0.230183077):0.8039164867):0.3730038516,t79:0.9049554598):0.008393579628,(t46:0.05780046433,(t6:0.1631849981,t81:0.5384514888):0.8815424331):0.7129145358):0.9547338288,((t37:0.6507847607,((t31:0.02607493685,t11:0.04472297546):0.3350441111,t13:0.3689027205):0.2309699834):0.2834663277,(t45:0.5464109634,t56:0.623756083):0.8225920612):0.3664859182):0.7387945433,(((((t57:0.7937401866,t97:0.5238280646):0.8410244284,t88:0.7108538796):0.9742555241,((t70:0.8903303037,t100:0.7963381845):0.04238839727,((t3:0.3148138926,(((t32:0.6826657334,t85:0.5152513909):0.2917491719,t5:0.8117499282):0.1037160209,t82:0.737436238):0.3962228657):0.6861416439,((t93:0.9698872003,t23:0.3531935748):0.9532342199,t61:0.03933861828):0.6078367378):0.8907218107):0.7618982717):0.2459336291,(((t2:0.6646096616,t43:0.2235449431):0.1258856647,t49:0.6103038758):0.708156846,t91:0.4830829434):0.8970741369):0.8757817054,t28:0.02927036583):0.586141187):0.4437033541):0.7749263747,(((t10:0.8176077947,t40:0.5686039282):0.2708641936,((t84:0.7844968771,t67:0.2490973365):0.1329914862,(t54:0.9972389175,(((t39:0.8772778097,t22:0.6142720587):0.4587370832,t21:0.6421265274):0.3914671526,(t71:0.9076540405,t51:0.5560081138):0.2632204019):0.835884884):0.784519562):0.4417822429):0.9865452512,(t44:0.5330150193,(((t73:0.05110091181,t15:0.3320527144):0.3521587751,(t65:0.05273675127,t33:0.9678427358):0.2134595928):0.9132453285,((t48:0.2262219382,t18:0.2089061344):0.8540851672,((t1:0.7213522177,t60:0.5620779982):0.6258146379,t12:0.05284262123):0.8915971264):0.8343131964):0.8695364469):0.5049135908):0.8888139946):0.7708704474):0.6049848683):0.04934649402):0.2062134217);"
my_tree = readTree(tree = tree_string)

freqaln = 'aln_aa_H1N1_HA.fasta'
size = 200
kappa = 4.5
fgen = ReadFreqs(by='amino', alnfile=freqaln) ## by needs to be the type of data in the file
freqs = fgen.getCodonFreqs()
omegas = np.linspace( 0.01, 3.0, num = 100 )


hydir = '/Users/sjspielman/Research/MutSel/PrelimNRSA/HyPhyMaterials/'
os.chdir(hydir)


# Save omegas
omegafile = results_dir + "omegas.txt"
file = open(omegafile, 'w')
file.write("true\thyphy\n")

for n in range(100):
	print n

	print "building model"
	partition = buildSingleModel(freqs, kappa, omegas[n], size)
	
	# sim output file
	seqfile = results_dir+"seqs"+str(n)+".fasta"
	
	# Simulate
	callSim(partition, my_tree, seqfile)
	
	# Call hyphy
	os.chdir(hydir)
	prepHyPhy(seqfile, tree_string)
	hyomega = runHyPhy()
	file.write(str(omegas[n])+'\t'+str(hyomega)+'\n')
	
file.close()