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

sys.path.append('../Simulator/src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *
from mutation_counter import *
from site_counter import *
from predict_omega import *
from functions_simandinf import *

################################ GLOBAL #######################################
molecules = Genetics()
treefile = "2.tre" # in here, all branch lengths are 0.5
my_tree = readTree(file=treefile)
length = 100000
codonParams = {'alpha': 1.0, 'mu': {'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}}
mutSelParams = {'mu': {'AC': 0.001, 'CA':0.001, 'AG': 0.001, 'GA':0.001, 'AT': 0.001, 'TA':0.001, 'CG': 0.001, 'GC':0.001, 'CT': 0.001, 'TC':0.001, 'GT': 0.001, 'TG':0.001}}
M = MutationCounter()
S = SiteCounter()
###############################################################################

hrh1_columns = [134, 155, 157, 158, 160, 161, 162, 164, 165, 166, 167, 168, 169, 172, 175, 179, 182, 183, 185, 186, 188, 200, 207, 208, 212, 213, 214, 215, 220, 221, 222, 231, 232, 233, 235, 239, 252, 257, 264, 268, 271, 275, 276, 277, 279, 280, 281, 283, 284, 285, 287, 288, 290, 291, 293, 297, 298, 300, 318, 320, 322, 324, 326, 327, 329, 330, 334, 337, 338, 342, 345, 346, 349, 353, 356, 360, 362, 363, 364, 369, 370, 371, 372, 373, 374, 384, 386, 387, 389, 408, 419, 420, 441, 442, 466, 529, 563, 566, 567, 573, 574, 576, 593, 607, 610, 612, 613, 616, 617, 618, 620, 636, 637, 638, 641, 645, 649, 652, 653, 655, 656, 659, 664, 666, 675, 678, 679, 683, 684, 685, 686, 693, 696, 700, 707, 711, 712, 715, 716, 717, 718, 719, 720]
numaa        = [6, 5, 3, 4, 7, 4, 4, 3, 3, 5, 5, 4, 4, 3, 7, 3, 4, 3, 4, 5, 5, 4, 3, 3, 5, 4, 3, 7, 4, 4, 3, 3, 4, 3, 6, 3, 3, 6, 3, 3, 5, 5, 6, 5, 5, 3, 4, 4, 4, 5, 3, 4, 3, 5, 6, 4, 3, 6, 5, 4, 4, 6, 3, 5, 5, 7, 6, 3, 3, 5, 6, 6, 5, 5, 7, 7, 6, 5, 6, 5, 7, 7, 5, 4, 6, 4, 7, 7, 6, 6, 3, 3, 5, 5, 7, 4, 5, 3, 3, 4, 3, 4, 4, 5, 5, 4, 5, 6, 5, 4, 3, 5, 5, 3, 3, 4, 3, 4, 3, 3, 3, 5, 5]
'''
cpu = sys.argv[1]
run = int(sys.argv[2])
rdir = sys.argv[3]
if rdir[-1] != '/':
	rdir += '/'
final_outfile = rdir + sys.argv[4]



outfile = rdir+"mutsel_"+str(run)+".txt"
seqfile = str(run)+'.fasta'
'''
# Set up model for simulation	
f = prepFreq(1.0, hrh1_columns[0])
mutSelParams['stateFreqs'] = f
model = Model()
model.params = mutSelParams
mat = mutSel_MatrixBuilder(model)
model.Q = mat.buildQ()
partitions = [(length, model)]		

# Evolve
myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
myEvolver.sim_sub_tree(my_tree)
myEvolver.writeSequences()


# Use math to derive an omega for the simulated sequences
derived_w = deriveAverageOmegaAlignment(seqfile, f)

# Get omega using counting method
records = list(SeqIO.parse(seqfile, 'fasta'))
s1 = records[0].seq
s2 = records[1].seq
( ns_mut, s_mut ) = M.countMutations( s1, s2 )
( ns_sites1, s_sites1 ) = S.countSites( s1 )
( ns_sites2, s_sites2 ) = S.countSites( s2 )
dS = 2*sum( s_mut )/(sum( s_sites1 ) + sum( s_sites2 ))
dN = 2*sum( ns_mut )/(sum( ns_sites2 ) + sum( ns_sites2 ))
count_w = dN/dS



# Send to HyPhy
setuphyphy1 = "cp "+seqfile+" temp.fasta"
setup1 = subprocess.call(setuphyphy1, shell = True)
assert(setup1 == 0), "couldn't create temp.fasta"

setuphyphy2 = "cat 2.tre >> temp.fasta"
setup2 = subprocess.call(setuphyphy2, shell = True)
assert(setup2 == 0), "couldn't add tree to hyphy infile"

hyf = freq2Hyphy(f)
setuphyphy3 = "sed 's/PLACEHOLDER/"+hyf+"/g' globalGY94.bf > run.bf"
setup3 = subprocess.call(setuphyphy3, shell = True)
assert(setup3 == 0), "couldn't properly add in frequencies"

hyphy = "./HYPHYMP globalGY94.bf CPU="+cpu+" > hyout.txt"
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
			
# Now results everything to file
outf = open(outfile, 'w')
outf.write(str(numaa[run]) + '\t' + str(hrh1_columns[run]) + '\t' + str(count_w) + '\t' + str(hyphy_w) + '\t' + str(derived_w) + '\n')
outf.close()

# And now send to the final outfile
save = "cat "+outfile+" >> "+final_outfile
saverun = subprocess.call(save, shell=True)
assert(saverun == 0), "couldn't save final file"
save2 = "cp "+seqfile+" "+rdir
saverun2 = subprocess.call(save2, shell=True)
assert(saverun2 == 0), "couldn't save seqfile"



