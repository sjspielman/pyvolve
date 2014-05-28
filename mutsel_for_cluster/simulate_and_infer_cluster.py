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
from mutation_counter import *
from site_counter import *
from predict_omega import *
from functions_simandinf import *

hrh1_columns = [134, 155, 157, 158, 160, 161, 162, 164, 165, 166, 167, 168, 169, 172, 175, 179, 182, 183, 185, 186, 188, 200, 207, 208, 212, 213, 214, 215, 220, 221, 222, 231, 232, 233, 235, 239, 252, 257, 264, 268, 271, 275, 276, 277, 279, 280, 281, 283, 284, 285, 287, 288, 290, 291, 293, 297, 298, 300, 318, 320, 322, 324, 326, 327, 329, 330, 334, 337, 338, 342, 345, 346, 349, 353, 356, 360, 362, 363, 364, 369, 370, 371, 372, 373, 374, 384, 386, 387, 389, 408, 419, 420, 441, 442, 466, 529, 563, 566, 567, 573, 574, 576, 593, 607, 610, 612, 613, 616, 617, 618, 620, 636, 637, 638, 641, 645, 649, 652, 653, 655, 656, 659, 664, 666, 675, 678, 679, 683, 684, 685, 686, 693, 696, 700, 707, 711, 712, 715, 716, 717, 718, 719, 720]
numaa        = [6, 5, 3, 4, 7, 4, 4, 3, 3, 5, 5, 4, 4, 3, 7, 3, 4, 3, 4, 5, 5, 4, 3, 3, 5, 4, 3, 7, 4, 4, 3, 3, 4, 3, 6, 3, 3, 6, 3, 3, 5, 5, 6, 5, 5, 3, 4, 4, 4, 5, 3, 4, 3, 5, 6, 4, 3, 6, 5, 4, 4, 6, 3, 5, 5, 7, 6, 3, 3, 5, 6, 6, 5, 5, 7, 7, 6, 5, 6, 5, 7, 7, 5, 4, 6, 4, 7, 7, 6, 6, 3, 3, 5, 5, 7, 4, 5, 3, 3, 4, 3, 4, 4, 5, 5, 4, 5, 6, 5, 4, 3, 5, 5, 3, 3, 4, 3, 4, 3, 3, 3, 5, 5]

cpu = sys.argv[1]
run = int(sys.argv[2])
treefile = sys.argv[3]
rdir = sys.argv[4]
if rdir[-1] != '/':
	rdir += '/'
final_outfile = rdir + sys.argv[5]

outfile = rdir+"mutsel_"+str(run)+".txt"
seqfile = str(run)+'.fasta'

# Simulate
print "simulating"
f = simulate(treefile, seqfile, 10000, 'mutsel', 'random', 'amino', 'codon')
#simulate(treefile, seqfile, length, modeltype, freq, by, type, columns = None, column_index = None, constraint = None)

# Use math to derive an omega for the simulated sequences
print "deriving"
derived_w = deriveAverageOmegaAlignment(seqfile, f)

# Nei-Gojobori Method
#print "nei-gojobori"
#nei_w = run_neigojo(seqfile)

# Send to PAML
#print "paml"
#paml_w = runpaml(seqfile)

# Send to HyPhy
print "hyphy"
hyphy_w = runhyphy("globalGY94.bf", seqfile, treefile, cpu, f)
		
# Now results everything to file
outf = open(outfile, 'w')
#outf.write(str(numaa[run]) + '\t' + str(hrh1_columns[run]) + '\t' + str(count_w) + '\t' + str(paml_w) + '\t' + str(derived_w) + '\n')
outf.write(str(derived_w) + '\t' + str(hyphy_w) + '\n')
outf.close()

# And now send to the final outfile
save = "cat "+outfile+" >> "+final_outfile
saverun = subprocess.call(save, shell=True)
assert(saverun == 0), "couldn't save final file"
save2 = "cp "+seqfile+" "+rdir
saverun2 = subprocess.call(save2, shell=True)
assert(saverun2 == 0), "couldn't save seqfile"



