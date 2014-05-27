## SJS. Functions for simulate_and_infer.py

import os
import re
import sys
import subprocess
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *

# FIX THIS
path_to_sim = '/Users/sjspielman/Research/MutSel/Simulator/'
sys.path.append(path_to_sim + 'src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


#number of amino acids per column:  4aa, 5aa, 6aa, 7aa, 8aa, 9aa, 10aa (from hrh1_aa.fasta). so, numaa = index + 4
constrained_propensities_columns = [175, 169, 326, 129, 327, 154, 121]
omegas_purifying = np.linspace(0.05, 1.00, num=20)
omegas_positive  = np.linspace(1.25, 4.25, num=20)

zero = 1e-8
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nslist = [['CAA', 'GAA', 'ACA', 'ATA', 'AGA', 'AAC', 'AAT'], ['CAC', 'TAC', 'GAC', 'ACC', 'ATC', 'AGC', 'AAA', 'AAG'], ['CAG', 'GAG', 'ACG', 'ATG', 'AGG', 'AAC', 'AAT'], ['CAT', 'TAT', 'GAT', 'ACT', 'ATT', 'AGT', 'AAA', 'AAG'], ['CCA', 'TCA', 'GCA', 'AAA', 'ATA', 'AGA'], ['CCC', 'TCC', 'GCC', 'AAC', 'ATC', 'AGC'], ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG'], ['CCT', 'TCT', 'GCT', 'AAT', 'ATT', 'AGT'], ['GGA', 'AAA', 'ACA', 'ATA', 'AGC', 'AGT'], ['CGC', 'TGC', 'GGC', 'AAC', 'ACC', 'ATC', 'AGA', 'AGG'], ['TGG', 'GGG', 'AAG', 'ACG', 'ATG', 'AGC', 'AGT'], ['CGT', 'TGT', 'GGT', 'AAT', 'ACT', 'ATT', 'AGA', 'AGG'], ['CTA', 'TTA', 'GTA', 'AAA', 'ACA', 'AGA', 'ATG'], ['CTC', 'TTC', 'GTC', 'AAC', 'ACC', 'AGC', 'ATG'], ['CTG', 'TTG', 'GTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'], ['CTT', 'TTT', 'GTT', 'AAT', 'ACT', 'AGT', 'ATG'], ['AAA', 'GAA', 'CCA', 'CTA', 'CGA', 'CAC', 'CAT'], ['AAC', 'TAC', 'GAC', 'CCC', 'CTC', 'CGC', 'CAA', 'CAG'], ['AAG', 'GAG', 'CCG', 'CTG', 'CGG', 'CAC', 'CAT'], ['AAT', 'TAT', 'GAT', 'CCT', 'CTT', 'CGT', 'CAA', 'CAG'], ['ACA', 'TCA', 'GCA', 'CAA', 'CTA', 'CGA'], ['ACC', 'TCC', 'GCC', 'CAC', 'CTC', 'CGC'], ['ACG', 'TCG', 'GCG', 'CAG', 'CTG', 'CGG'], ['ACT', 'TCT', 'GCT', 'CAT', 'CTT', 'CGT'], ['GGA', 'CAA', 'CCA', 'CTA'], ['AGC', 'TGC', 'GGC', 'CAC', 'CCC', 'CTC'], ['TGG', 'GGG', 'CAG', 'CCG', 'CTG'], ['AGT', 'TGT', 'GGT', 'CAT', 'CCT', 'CTT'], ['ATA', 'GTA', 'CAA', 'CCA', 'CGA'], ['ATC', 'TTC', 'GTC', 'CAC', 'CCC', 'CGC'], ['ATG', 'GTG', 'CAG', 'CCG', 'CGG'], ['ATT', 'TTT', 'GTT', 'CAT', 'CCT', 'CGT'], ['AAA', 'CAA', 'GCA', 'GTA', 'GGA', 'GAC', 'GAT'], ['AAC', 'CAC', 'TAC', 'GCC', 'GTC', 'GGC', 'GAA', 'GAG'], ['AAG', 'CAG', 'GCG', 'GTG', 'GGG', 'GAC', 'GAT'], ['AAT', 'CAT', 'TAT', 'GCT', 'GTT', 'GGT', 'GAA', 'GAG'], ['ACA', 'CCA', 'TCA', 'GAA', 'GTA', 'GGA'], ['ACC', 'CCC', 'TCC', 'GAC', 'GTC', 'GGC'], ['ACG', 'CCG', 'TCG', 'GAG', 'GTG', 'GGG'], ['ACT', 'CCT', 'TCT', 'GAT', 'GTT', 'GGT'], ['AGA', 'CGA', 'GAA', 'GCA', 'GTA'], ['AGC', 'CGC', 'TGC', 'GAC', 'GCC', 'GTC'], ['AGG', 'CGG', 'TGG', 'GAG', 'GCG', 'GTG'], ['AGT', 'CGT', 'TGT', 'GAT', 'GCT', 'GTT'], ['ATA', 'CTA', 'TTA', 'GAA', 'GCA', 'GGA'], ['ATC', 'CTC', 'TTC', 'GAC', 'GCC', 'GGC'], ['ATG', 'CTG', 'TTG', 'GAG', 'GCG', 'GGG'], ['ATT', 'CTT', 'TTT', 'GAT', 'GCT', 'GGT'], ['AAC', 'CAC', 'GAC', 'TCC', 'TTC', 'TGC'], ['AAT', 'CAT', 'GAT', 'TCT', 'TTT', 'TGT'], ['ACA', 'CCA', 'GCA', 'TTA'], ['ACC', 'CCC', 'GCC', 'TAC', 'TTC', 'TGC'], ['ACG', 'CCG', 'GCG', 'TTG', 'TGG'], ['ACT', 'CCT', 'GCT', 'TAT', 'TTT', 'TGT'], ['AGC', 'CGC', 'GGC', 'TAC', 'TCC', 'TTC', 'TGG'], ['AGG', 'CGG', 'GGG', 'TCG', 'TTG', 'TGC', 'TGT'], ['AGT', 'CGT', 'GGT', 'TAT', 'TCT', 'TTT', 'TGG'], ['ATA', 'GTA', 'TCA', 'TTC', 'TTT'], ['ATC', 'CTC', 'GTC', 'TAC', 'TCC', 'TGC', 'TTA', 'TTG'], ['ATG', 'GTG', 'TCG', 'TGG', 'TTC', 'TTT'], ['ATT', 'CTT', 'GTT', 'TAT', 'TCT', 'TGT', 'TTA', 'TTG']]



def freq2Hyphy(f):
	''' Convert codon frequencies to a form hyphy can use. '''
	hyphy_f = "{"
	for freq in f:
		hyphy_f += "{"
		hyphy_f += str(freq)
		hyphy_f += "},"
	hyphy_f = hyphy_f[:-1]
	hyphy_f += "}"
	return hyphy_f

#################################################################################################################################
## Given F(i) and F(j), where F() is the frequency of the given codon in that column, return fix_(i->j).
def fix(fi, fj):
	if fi == fj:
		return 1.
	elif fi == 0  or fj == 0:
		return 0.
	else:
		return (np.log(fj) - np.log(fi)) / (1 - fi/fj)
#################################################################################################################################


#################################################################################################################################
def deriveOmegaColumn(codonFreq, nonZero):
	''' Use the math to derive an omega at a given position ''' 
	kN=0 #dN numerator
	nN=0 #dN denominator. Does not consider number of nonsyn options
	fix_sum=0

	# Calculations
	for i in nonZero:
		fix_sum=0
		
		### Nonsynonymous.
		for nscodon in nslist[i]:
			nscodon_freq = codonFreq[codons.index(nscodon)]
			fix_sum += fix(float(codonFreq[i]), float(nscodon_freq))					
			nN += codonFreq[i]
		kN += fix_sum*codonFreq[i]

	# Final dN/dS
	if kN < zero:
		dNdS = 0
	else:
		dNdS=kN/nN
	
	return dNdS
#################################################################################################################################

	
#################################################################################################################################
def readAln(seqfile):
	aln=[] #list of lists wherein each nested list is a row
	aln_raw = list(SeqIO.parse(seqfile, 'fasta'))
	for record in aln_raw:
		aln.append(str(record.seq))
	
	alnlen=len(aln[0])
	numseq=len(aln)
	return (aln, alnlen, numseq)
#################################################################################################################################


#################################################################################################################################
def getNonZeroFreqs(freq):
	''' Return indices whose frequencies are not 0.'''
	nonZero = [] # indices
	for i in range(len(freq)):
		if freq[i] > zero:
			nonZero.append(i)
	return nonZero
	
def getFreqs(numseq, aln, col):	
	''' Return indices whose frequencies are not 0, based on codon frequencies seen in the simulated data - NOT the true frequencies.'''
	
	codonFreq=np.zeros(len(codons)) # will contain frequencies for all codons in a given column		
	#Find codon frequencies
	for row in range(numseq):
		codon=aln[row][col:col+3]
		codonFreq[codons.index(codon)] += 1	
	codonFreq/=float(numseq)
		
	# Fill nonZero with codonFreq indices whose values are not 0
	nonZero = getNonZeroFreqs(codonFreq)
	
	return codonFreq,nonZero
#################################################################################################################################

#################################################################################################################################
   
#################################################################################################################################
def deriveAverageOmegaAlignment(alnfile, globalFreqs = None):
	''' Calculate and average derived omegas for each column in an alignment. Return a single global derived omega for that alignment.'''
	aln, alnlen, numseq = readAln(alnfile)
	
	
	# THIS LINE IS FOR GLOBAL FREQUENCIES
	if globalFreqs is not None:
		nonZero = getNonZeroFreqs(globalFreqs)
		dNdS = deriveOmegaColumn(globalFreqs, nonZero)
		return dNdS
	
	else:	
		dNdS_values = np.zeros(alnlen)
		count = 0
		for col in range(0,alnlen,3):
			columnFreqs, nonZero = getFreqs(numseq, aln, col)
			dNdS = deriveOmegaColumn(columnFreqs, nonZero)
			dNdS_values[count] = float(dNdS)
			count += 1
		mean_dNdS = np.mean(dNdS_values)
		return mean_dNdS
#################################################################################################################################


#################################################################################################################################
def prepFreq(freq, column):
	''' Prepare equilibrium frequency object'''
	
	if freq == 'equal':
		freqObject = EqualFreqs(by = 'amino', type = 'codon')
	else:
		assert(type(freq) is float and freq > 0. and freq <=1.), "Booo learn to type."
		freqObject = ReadFreqs(by = 'amino', type = 'codon', columns=[column], file = 'hrh1_aa.fasta', constraint = freq)	
	f = freqObject.calcFreqs()
	return f
#################################################################################################################################




#################################################################################################################################
def runOmegaModel(omegas, cols, outfile, freq):
	''' Simulate and infer for omega models.'''
	
	outf = open(outfile, 'w')
	outf.write('num_aa\tsimulated_w\tderived_w\tkN\n')
	
	for i in cols:
		numaa = str(constrained_propensities_columns.index(i) + 4)
		f = prepFreq(freq, i)
		
		for w in omegas:
		
			print "running", numaa, "amino acids, with omega = ", w
			# Construct model for these conditions
			codonParams['stateFreqs'] = f
			codonParams['beta'] = w
			model = Model()
			model.params = codonParams
			mat = codon_MatrixBuilder(model)
			model.Q = mat.buildQ()
			partitions = [(length, model)]
		
			# Evolve
			print "Evolving"
			seqfile = "out.fasta"	
			myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
			myEvolver.sim_sub_tree(my_tree)
			myEvolver.writeSequences()
	
			# Use math to derive an omega for the simulated sequences
			print "Deriving"
			derivedw, kN = deriveAverageOmegaAlignment(seqfile, f)
	
			# And save it all
			outf.write(numaa+'\t'+str(w)+'\t'+str(derivedw)+'\t'+str(kN)+'\n')
	outf.close()
#################################################################################################################################
