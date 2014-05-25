## SJS 5/24/14. 
## Script to run various simulations for NIH


import os
import re
import sys
import subprocess
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *

sys.path.append('Simulator/src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

from functions_simandinf import *

################################ GLOBAL #######################################
molecules = Genetics()
treefile = "100.tre" # in here, all branch lengths are 0.5
my_tree = readTree(file=treefile)
length = 100
codonParams = {'alpha': 1.0, 'mu': {'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}}
mutSelParams = {'mu': {'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
###############################################################################

hrh1_columns = [134, 155, 157, 158, 160, 161, 162, 164, 165, 166, 167, 168, 169, 172, 175, 179, 182, 183, 185, 186, 188, 200, 207, 208, 212, 213, 214, 215, 220, 221, 222, 231, 232, 233, 235, 239, 252, 257, 264, 268, 271, 275, 276, 277, 279, 280, 281, 283, 284, 285, 287, 288, 290, 291, 293, 297, 298, 300, 318, 320, 322, 324, 326, 327, 329, 330, 334, 337, 338, 342, 345, 346, 349, 353, 356, 360, 362, 363, 364, 369, 370, 371, 372, 373, 374, 384, 386, 387, 389, 408, 419, 420, 441, 442, 466, 529, 563, 566, 567, 573, 574, 576, 593, 607, 610, 612, 613, 616, 617, 618, 620, 636, 637, 638, 641, 645, 649, 652, 653, 655, 656, 659, 664, 666, 675, 678, 679, 683, 684, 685, 686, 693, 696, 700, 707, 711, 712, 715, 716, 717, 718, 719, 720]
numaa        = [6, 5, 3, 4, 7, 4, 4, 3, 3, 5, 5, 4, 4, 3, 7, 3, 4, 3, 4, 5, 5, 4, 3, 3, 5, 4, 3, 7, 4, 4, 3, 3, 4, 3, 6, 3, 3, 6, 3, 3, 5, 5, 6, 5, 5, 3, 4, 4, 4, 5, 3, 4, 3, 5, 6, 4, 3, 6, 5, 4, 4, 6, 3, 5, 5, 7, 6, 3, 3, 5, 6, 6, 5, 5, 7, 7, 6, 5, 6, 5, 7, 7, 5, 4, 6, 4, 7, 7, 6, 6, 3, 3, 5, 5, 7, 4, 5, 3, 3, 4, 3, 4, 4, 5, 5, 4, 5, 6, 5, 4, 3, 5, 5, 3, 3, 4, 3, 4, 3, 3, 3, 5, 5]

run = int(sys.argv[1])
final_outfile = sys.argv[2]

outfile = "mutsel_"+str(run)+".txt"

# Set up model for simulation	
f = prepFreq(1.0, i)
mutSelParams['stateFreqs'] = f
model = Model()
model.params = mutSelParams
mat = mutSel_MatrixBuilder(model)
model.Q = mat.buildQ()
partitions = [(length, model)]		

# Evolve
print "Evolving"
seqfile = 'temp.fasta'	
myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeSequences()


# Use math to derive an omega for the simulated sequences
print "Deriving"
derived_w = deriveAverageOmegaAlignment(seqfile, f)


# Send to HyPhy
print "hyphy"
setuphyphy = "cat 100.tre >> temp.fasta"
setup = subprocess.call(setuphyphy, shell = True)
assert(setup == 0), "couldn't prep file for hyphy"
hyphy = "./HYPHYMP globalGY94.bf > hyout.txt"
runhyphy = subprocess.call(hyphy, shell = True)
assert (runhyphy == 0), "hyphy fail"

# grab hyphy output
hyout = open('hyout.txt', 'r')
hylines = hyout.readlines()
hyout.close()
for line in hylines:
	findw = re.search("^w=(\d+\.*\d*)", line)
	if findw:
		hyphy_w = findw.group(1)
		break
			
# Now save everything to file
outf = open(outfile, 'w')
outf.write(str(numaa[run]) + '\t' + str(hrh1_columns[run]) + '\t' + str(hyphy_w) + '\t' + str(derived_w) + '\n')
outf.close()

# And now send to the final outfile
save = "cat "+outfile+" >> "+final_outfile
saverun = subprocess.call(save, shell=True)
assert(saverun == 0), "couldn't save final file"

