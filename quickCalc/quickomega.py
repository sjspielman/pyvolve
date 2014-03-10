###### Rewritten for my simulator.

import os
import re
import sys
import subprocess
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *

sys.path.append("/Users/sjspielman/Omega_MutSel/Simulator/src/")

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

def prepSimColumns(freq_aln, whichcol, kappa, omegas, numPart, partLen):
	''' Use when frequencies are different (column-specific) for each partition ''' 
	
	partitions = []
	
	print "constructing models for", numPart, "partitions"
	for i in range(numPart):
		# Define model object and its parameters. Build model matrix. Add tuple (partition length, model) to partitions list
		model = misc.Model()
		
		# For site-specific distribution
		fgen = ReadFreqs(by='amino', alnfile=freq_aln, which = whichcol[i]) ## by needs to be the type of data in the file
		freqs = fgen.getCodonFreqs()
		
		model.params = { "kappa": kappa, "omega": omegas[i], "stateFreqs": freqs }
		m = GY94(model)
		model.Q = m.buildQ()
		partitions.append( (partLen, model) )

	return partitions
	
	
	
	
def prepSim(freq_aln, kappa, omegas, numPart, partLen):
	''' Use when frequencies are global and not site-specific. All sites share same freqs. ''' 
	partitions = []
	
	# global param
	#fgen = ReadFreqs(by='amino', alnfile=freq_aln) ## by needs to be the type of data in the file
	
	## UNIFORM DISTRIBUTION
	fgen = EqualFreqs(by='amino')
	
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
#################################################################################################################################
def callSim(partitions, my_tree, seqfile):

	print "evolving"
	myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
	myEvolver.sim_sub_tree(my_tree)
	myEvolver.writeSequences()
#################################################################################################################################

#################################################################################################################################
## For a given codon, returns a list of its synonymous changes and of its nonsynonymous changes, with only one nucleotide change permitted.
def findSynNonsyn(codon_raw):
	syn=[]
	nonsyn=[]
	codon=Seq(codon_raw, generic_dna)
	aa_source=str(codon.translate())
	for n in range(0,3):
		for nuc in ['A', 'C', 'T', 'G']:
			#Don't count the changes to the same nucleotide. upper as in upper case.
			if codon[n].upper()==nuc:
				continue
			target = codon[0:n]+nuc+codon[n+1:3]
			aa_target = str(target.translate())
			if aa_target=='*':	#disregard nonsense mutations
				continue
			if aa_source==aa_target:
				syn.append(str(target))
			elif aa_source != target:
				nonsyn.append(str(target))
	return(syn, nonsyn)
#################################################################################################################################


#################################################################################################################################
## Given F(i) and F(j), where F() is the frequency of the given codon in that column, return fix_(i->j).
def fix(fi, fj):
	return (np.log(fj) - np.log(fi)) / (1 - fi/fj)
#################################################################################################################################

#################################################################################################################################
### Given a source codon and target codon, return True if the change is nonsynonymous. probably not needed.
def isNonsyn(codon_source_raw, codon_target_raw):

	# Convert to biopython object so can use the translate() function
	codon_source = Seq(codon_source_raw, generic_dna)
	codon_target = Seq(codon_target_raw, generic_dna)

	if (codon_source.translate() == codon_target.translate()):
		return False
	else:
		return True
#################################################################################################################################

#################################################################################################################################
def deriveOmegaColumn(seqfile, 

	codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
	slist=[['AAG'], ['AAT'], ['AAA'], ['AAC'], ['ACC', 'ACT', 'ACG'], ['ACA', 'ACT', 'ACG'], ['ACA', 'ACC', 'ACT'], ['ACA', 'ACC', 'ACG'], ['CGA', 'AGG'], ['AGT'], ['CGG', 'AGA'], ['AGC'], ['ATC', 'ATT'], ['ATA', 'ATT'], [], ['ATA', 'ATC'], ['CAG'], ['CAT'], ['CAA'], ['CAC'], ['CCC', 'CCT', 'CCG'], ['CCA', 'CCT', 'CCG'], ['CCA', 'CCC', 'CCT'], ['CCA', 'CCC', 'CCG'], ['AGA', 'CGC', 'CGT', 'CGG'], ['CGA', 'CGT', 'CGG'], ['AGG', 'CGA', 'CGC', 'CGT'], ['CGA', 'CGC', 'CGG'], ['TTA', 'CTC', 'CTT', 'CTG'], ['CTA', 'CTT', 'CTG'], ['TTG', 'CTA', 'CTC', 'CTT'], ['CTA', 'CTC', 'CTG'], ['GAG'], ['GAT'], ['GAA'], ['GAC'], ['GCC', 'GCT', 'GCG'], ['GCA', 'GCT', 'GCG'], ['GCA', 'GCC', 'GCT'], ['GCA', 'GCC', 'GCG'], ['GGC', 'GGT', 'GGG'], ['GGA', 'GGT', 'GGG'], ['GGA', 'GGC', 'GGT'], ['GGA', 'GGC', 'GGG'], ['GTC', 'GTT', 'GTG'], ['GTA', 'GTT', 'GTG'], ['GTA', 'GTC', 'GTT'], ['GTA', 'GTC', 'GTG'], ['TAT'], ['TAC'], ['TCC', 'TCT', 'TCG'], ['TCA', 'TCT', 'TCG'], ['TCA', 'TCC', 'TCT'], ['TCA', 'TCC', 'TCG'], ['TGT'], [], ['TGC'], ['CTA', 'TTG'], ['TTT'], ['CTG', 'TTA'], ['TTC']]
	nslist = [['CAA', 'GAA', 'ACA', 'ATA', 'AGA', 'AAC', 'AAT'], ['CAC', 'TAC', 'GAC', 'ACC', 'ATC', 'AGC', 'AAA', 'AAG'], ['CAG', 'GAG', 'ACG', 'ATG', 'AGG', 'AAC', 'AAT'], ['CAT', 'TAT', 'GAT', 'ACT', 'ATT', 'AGT', 'AAA', 'AAG'], ['CCA', 'TCA', 'GCA', 'AAA', 'ATA', 'AGA'], ['CCC', 'TCC', 'GCC', 'AAC', 'ATC', 'AGC'], ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG'], ['CCT', 'TCT', 'GCT', 'AAT', 'ATT', 'AGT'], ['GGA', 'AAA', 'ACA', 'ATA', 'AGC', 'AGT'], ['CGC', 'TGC', 'GGC', 'AAC', 'ACC', 'ATC', 'AGA', 'AGG'], ['TGG', 'GGG', 'AAG', 'ACG', 'ATG', 'AGC', 'AGT'], ['CGT', 'TGT', 'GGT', 'AAT', 'ACT', 'ATT', 'AGA', 'AGG'], ['CTA', 'TTA', 'GTA', 'AAA', 'ACA', 'AGA', 'ATG'], ['CTC', 'TTC', 'GTC', 'AAC', 'ACC', 'AGC', 'ATG'], ['CTG', 'TTG', 'GTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'], ['CTT', 'TTT', 'GTT', 'AAT', 'ACT', 'AGT', 'ATG'], ['AAA', 'GAA', 'CCA', 'CTA', 'CGA', 'CAC', 'CAT'], ['AAC', 'TAC', 'GAC', 'CCC', 'CTC', 'CGC', 'CAA', 'CAG'], ['AAG', 'GAG', 'CCG', 'CTG', 'CGG', 'CAC', 'CAT'], ['AAT', 'TAT', 'GAT', 'CCT', 'CTT', 'CGT', 'CAA', 'CAG'], ['ACA', 'TCA', 'GCA', 'CAA', 'CTA', 'CGA'], ['ACC', 'TCC', 'GCC', 'CAC', 'CTC', 'CGC'], ['ACG', 'TCG', 'GCG', 'CAG', 'CTG', 'CGG'], ['ACT', 'TCT', 'GCT', 'CAT', 'CTT', 'CGT'], ['GGA', 'CAA', 'CCA', 'CTA'], ['AGC', 'TGC', 'GGC', 'CAC', 'CCC', 'CTC'], ['TGG', 'GGG', 'CAG', 'CCG', 'CTG'], ['AGT', 'TGT', 'GGT', 'CAT', 'CCT', 'CTT'], ['ATA', 'GTA', 'CAA', 'CCA', 'CGA'], ['ATC', 'TTC', 'GTC', 'CAC', 'CCC', 'CGC'], ['ATG', 'GTG', 'CAG', 'CCG', 'CGG'], ['ATT', 'TTT', 'GTT', 'CAT', 'CCT', 'CGT'], ['AAA', 'CAA', 'GCA', 'GTA', 'GGA', 'GAC', 'GAT'], ['AAC', 'CAC', 'TAC', 'GCC', 'GTC', 'GGC', 'GAA', 'GAG'], ['AAG', 'CAG', 'GCG', 'GTG', 'GGG', 'GAC', 'GAT'], ['AAT', 'CAT', 'TAT', 'GCT', 'GTT', 'GGT', 'GAA', 'GAG'], ['ACA', 'CCA', 'TCA', 'GAA', 'GTA', 'GGA'], ['ACC', 'CCC', 'TCC', 'GAC', 'GTC', 'GGC'], ['ACG', 'CCG', 'TCG', 'GAG', 'GTG', 'GGG'], ['ACT', 'CCT', 'TCT', 'GAT', 'GTT', 'GGT'], ['AGA', 'CGA', 'GAA', 'GCA', 'GTA'], ['AGC', 'CGC', 'TGC', 'GAC', 'GCC', 'GTC'], ['AGG', 'CGG', 'TGG', 'GAG', 'GCG', 'GTG'], ['AGT', 'CGT', 'TGT', 'GAT', 'GCT', 'GTT'], ['ATA', 'CTA', 'TTA', 'GAA', 'GCA', 'GGA'], ['ATC', 'CTC', 'TTC', 'GAC', 'GCC', 'GGC'], ['ATG', 'CTG', 'TTG', 'GAG', 'GCG', 'GGG'], ['ATT', 'CTT', 'TTT', 'GAT', 'GCT', 'GGT'], ['AAC', 'CAC', 'GAC', 'TCC', 'TTC', 'TGC'], ['AAT', 'CAT', 'GAT', 'TCT', 'TTT', 'TGT'], ['ACA', 'CCA', 'GCA', 'TTA'], ['ACC', 'CCC', 'GCC', 'TAC', 'TTC', 'TGC'], ['ACG', 'CCG', 'GCG', 'TTG', 'TGG'], ['ACT', 'CCT', 'GCT', 'TAT', 'TTT', 'TGT'], ['AGC', 'CGC', 'GGC', 'TAC', 'TCC', 'TTC', 'TGG'], ['AGG', 'CGG', 'GGG', 'TCG', 'TTG', 'TGC', 'TGT'], ['AGT', 'CGT', 'GGT', 'TAT', 'TCT', 'TTT', 'TGG'], ['ATA', 'GTA', 'TCA', 'TTC', 'TTT'], ['ATC', 'CTC', 'GTC', 'TAC', 'TCC', 'TGC', 'TTA', 'TTG'], ['ATG', 'GTG', 'TCG', 'TGG', 'TTC', 'TTT'], ['ATT', 'CTT', 'GTT', 'TAT', 'TCT', 'TGT', 'TTA', 'TTG']]
	kN=0 #dN numerator
	nN=0 #dN denominator. Does not consider number of nonsyn options
	fix_sum=0
	
	codonFreq, nonZero = getFreq(codons, numseq, aln)

	# Calculations
	for i in nonZero:
		fix_sum=0
		
		### Nonsynonymous.
		for nscodon in nslist[i]:
			nscodon_freq = codonFreq[codons.index(nscodon)]
			if nscodon_freq==0 or codonFreq[i]==nscodon_freq:
				continue
			else:
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
def deriveOmegaAlignment(seqfile, outfile):
	aln, alnlen, numseq = readAln(seqfile)
	
	out = open(outfile, 'w')
	position = 1
	for i in range(0,alnlen,3):
		dNdS = deriveOmegaColumn(aln, numseq, alnlen)
		out.write(str(position)+'\t'+str(dNdS)+'\n')
		position+=1	
	out.close()
#################################################################################################################################

#################################################################################################################################
#################################################################################################################################



#################################################################################################################################
def getFreq(codons, numseq, aln):	
	codonFreq=np.zeros(len(codons)) # will contain frequencies for all codons in a given column
	nonZero=[] # will contain the nonzero indices for codonFreq
		
	#Find codon frequencies
	for row in range(numseq):
		codon=aln[row][col:col+3]
		codonFreq[codons.index(codon)] += 1	
	codonFreq/=float(numseq)
		
	# Fill nonZero with codonFreq indices whose values are not 0
	for i in range(len(codonFreq)):
		if codonFreq[i] > zero:
			nonZero.append(i)
	
	return codonFreq,nonZero
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
#################################################################################################################################	
# List of all codons. Note that the slist and nslist lists are also in this order.
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
slist=[['AAG'], ['AAT'], ['AAA'], ['AAC'], ['ACC', 'ACT', 'ACG'], ['ACA', 'ACT', 'ACG'], ['ACA', 'ACC', 'ACT'], ['ACA', 'ACC', 'ACG'], ['CGA', 'AGG'], ['AGT'], ['CGG', 'AGA'], ['AGC'], ['ATC', 'ATT'], ['ATA', 'ATT'], [], ['ATA', 'ATC'], ['CAG'], ['CAT'], ['CAA'], ['CAC'], ['CCC', 'CCT', 'CCG'], ['CCA', 'CCT', 'CCG'], ['CCA', 'CCC', 'CCT'], ['CCA', 'CCC', 'CCG'], ['AGA', 'CGC', 'CGT', 'CGG'], ['CGA', 'CGT', 'CGG'], ['AGG', 'CGA', 'CGC', 'CGT'], ['CGA', 'CGC', 'CGG'], ['TTA', 'CTC', 'CTT', 'CTG'], ['CTA', 'CTT', 'CTG'], ['TTG', 'CTA', 'CTC', 'CTT'], ['CTA', 'CTC', 'CTG'], ['GAG'], ['GAT'], ['GAA'], ['GAC'], ['GCC', 'GCT', 'GCG'], ['GCA', 'GCT', 'GCG'], ['GCA', 'GCC', 'GCT'], ['GCA', 'GCC', 'GCG'], ['GGC', 'GGT', 'GGG'], ['GGA', 'GGT', 'GGG'], ['GGA', 'GGC', 'GGT'], ['GGA', 'GGC', 'GGG'], ['GTC', 'GTT', 'GTG'], ['GTA', 'GTT', 'GTG'], ['GTA', 'GTC', 'GTT'], ['GTA', 'GTC', 'GTG'], ['TAT'], ['TAC'], ['TCC', 'TCT', 'TCG'], ['TCA', 'TCT', 'TCG'], ['TCA', 'TCC', 'TCT'], ['TCA', 'TCC', 'TCG'], ['TGT'], [], ['TGC'], ['CTA', 'TTG'], ['TTT'], ['CTG', 'TTA'], ['TTC']]
nslist = [['CAA', 'GAA', 'ACA', 'ATA', 'AGA', 'AAC', 'AAT'], ['CAC', 'TAC', 'GAC', 'ACC', 'ATC', 'AGC', 'AAA', 'AAG'], ['CAG', 'GAG', 'ACG', 'ATG', 'AGG', 'AAC', 'AAT'], ['CAT', 'TAT', 'GAT', 'ACT', 'ATT', 'AGT', 'AAA', 'AAG'], ['CCA', 'TCA', 'GCA', 'AAA', 'ATA', 'AGA'], ['CCC', 'TCC', 'GCC', 'AAC', 'ATC', 'AGC'], ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG'], ['CCT', 'TCT', 'GCT', 'AAT', 'ATT', 'AGT'], ['GGA', 'AAA', 'ACA', 'ATA', 'AGC', 'AGT'], ['CGC', 'TGC', 'GGC', 'AAC', 'ACC', 'ATC', 'AGA', 'AGG'], ['TGG', 'GGG', 'AAG', 'ACG', 'ATG', 'AGC', 'AGT'], ['CGT', 'TGT', 'GGT', 'AAT', 'ACT', 'ATT', 'AGA', 'AGG'], ['CTA', 'TTA', 'GTA', 'AAA', 'ACA', 'AGA', 'ATG'], ['CTC', 'TTC', 'GTC', 'AAC', 'ACC', 'AGC', 'ATG'], ['CTG', 'TTG', 'GTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'], ['CTT', 'TTT', 'GTT', 'AAT', 'ACT', 'AGT', 'ATG'], ['AAA', 'GAA', 'CCA', 'CTA', 'CGA', 'CAC', 'CAT'], ['AAC', 'TAC', 'GAC', 'CCC', 'CTC', 'CGC', 'CAA', 'CAG'], ['AAG', 'GAG', 'CCG', 'CTG', 'CGG', 'CAC', 'CAT'], ['AAT', 'TAT', 'GAT', 'CCT', 'CTT', 'CGT', 'CAA', 'CAG'], ['ACA', 'TCA', 'GCA', 'CAA', 'CTA', 'CGA'], ['ACC', 'TCC', 'GCC', 'CAC', 'CTC', 'CGC'], ['ACG', 'TCG', 'GCG', 'CAG', 'CTG', 'CGG'], ['ACT', 'TCT', 'GCT', 'CAT', 'CTT', 'CGT'], ['GGA', 'CAA', 'CCA', 'CTA'], ['AGC', 'TGC', 'GGC', 'CAC', 'CCC', 'CTC'], ['TGG', 'GGG', 'CAG', 'CCG', 'CTG'], ['AGT', 'TGT', 'GGT', 'CAT', 'CCT', 'CTT'], ['ATA', 'GTA', 'CAA', 'CCA', 'CGA'], ['ATC', 'TTC', 'GTC', 'CAC', 'CCC', 'CGC'], ['ATG', 'GTG', 'CAG', 'CCG', 'CGG'], ['ATT', 'TTT', 'GTT', 'CAT', 'CCT', 'CGT'], ['AAA', 'CAA', 'GCA', 'GTA', 'GGA', 'GAC', 'GAT'], ['AAC', 'CAC', 'TAC', 'GCC', 'GTC', 'GGC', 'GAA', 'GAG'], ['AAG', 'CAG', 'GCG', 'GTG', 'GGG', 'GAC', 'GAT'], ['AAT', 'CAT', 'TAT', 'GCT', 'GTT', 'GGT', 'GAA', 'GAG'], ['ACA', 'CCA', 'TCA', 'GAA', 'GTA', 'GGA'], ['ACC', 'CCC', 'TCC', 'GAC', 'GTC', 'GGC'], ['ACG', 'CCG', 'TCG', 'GAG', 'GTG', 'GGG'], ['ACT', 'CCT', 'TCT', 'GAT', 'GTT', 'GGT'], ['AGA', 'CGA', 'GAA', 'GCA', 'GTA'], ['AGC', 'CGC', 'TGC', 'GAC', 'GCC', 'GTC'], ['AGG', 'CGG', 'TGG', 'GAG', 'GCG', 'GTG'], ['AGT', 'CGT', 'TGT', 'GAT', 'GCT', 'GTT'], ['ATA', 'CTA', 'TTA', 'GAA', 'GCA', 'GGA'], ['ATC', 'CTC', 'TTC', 'GAC', 'GCC', 'GGC'], ['ATG', 'CTG', 'TTG', 'GAG', 'GCG', 'GGG'], ['ATT', 'CTT', 'TTT', 'GAT', 'GCT', 'GGT'], ['AAC', 'CAC', 'GAC', 'TCC', 'TTC', 'TGC'], ['AAT', 'CAT', 'GAT', 'TCT', 'TTT', 'TGT'], ['ACA', 'CCA', 'GCA', 'TTA'], ['ACC', 'CCC', 'GCC', 'TAC', 'TTC', 'TGC'], ['ACG', 'CCG', 'GCG', 'TTG', 'TGG'], ['ACT', 'CCT', 'GCT', 'TAT', 'TTT', 'TGT'], ['AGC', 'CGC', 'GGC', 'TAC', 'TCC', 'TTC', 'TGG'], ['AGG', 'CGG', 'GGG', 'TCG', 'TTG', 'TGC', 'TGT'], ['AGT', 'CGT', 'GGT', 'TAT', 'TCT', 'TTT', 'TGG'], ['ATA', 'GTA', 'TCA', 'TTC', 'TTT'], ['ATC', 'CTC', 'GTC', 'TAC', 'TCC', 'TGC', 'TTA', 'TTG'], ['ATG', 'GTG', 'TCG', 'TGG', 'TTC', 'TTT'], ['ATT', 'CTT', 'GTT', 'TAT', 'TCT', 'TGT', 'TTA', 'TTG']]

######################## code to generate slist, nslist ################################
#slist=[] #Contains syn codons for the given index. Indices as in list, codons
#nslist=[] #Contains nonsyn codons for the given index. Indices as in list, codons
#for codon in codons:
#	(syn, nonsyn) = findSynNonsyn(codon)
#	slist.append(syn)
#	nslist.append(nonsyn)

#################################################################################################################################
#################################################################################################################################	


## Conduct 100 simulations to confirm that the correlation between ML dN/dS and our derived dN/dS really exists

zero=1e-10


freq_aln = '/Users/sjspielman/structural_prediction_of_ER/evolutionary_rates/alignments/aminoacid/aln_aa_H1N1_HA.fasta'
######## THESE VALUES ARE FOR SITE-SPECIFIC STUFF!! ############
### Taken from structural_prediction_of_ER/FINAL_evolutionary_rates/siterates_REL/H1N1_HA.txt. All sites in these dicts are indexed at 0
# dict = {position:omega}...
#pos = {2: 1.62, 110:1.61, 97: 1.597, 239: 1.15, 146: 1.48, 278: 1.247, 88: 1.61, 157: 1.61, 158:1.29, 581: 3.27 } # sites indexed at 0. EXCLUDE ANY WITH DNDS=2.21. those sites are crazy.
#pur = {1: 0.18, 3: 0.24, 4: 0.4, 165: 0.66, 73: 0.308, 11: 0.57, 178: 0.725, 214: 0.47, 282: 0.89, 5: 0.817} # sites indexed at 0
#omegas = []
#columns = []
#for entry in pos:
#	omegas.append(pos[entry])
#	columns.append(entry)
################################################################



############# USE THIS FOR GLOBAL EQUILIBRIUM FREQUENCIES ##############
numPart = 15
#omegas = np.linspace(0.05, 0.85, num=numPart)  ####### PURIFYING SELECTION
omegas = np.linspace(1.1, 2.5, num=numPart)  ####### POSITIVE  SELECTION

#numPart = len(omegas)
partLen = 20
kappa = 4.5


results_dir='/Users/sjspielman/Dropbox/MutSelProject/quickCalc/SimSeqs_pos_unif/'
ensure_dir(results_dir)

#### Write a truerates files. Note that this will apply to all simulations.
truefile = results_dir+"truerates.txt"
truef = open(truefile, 'w')
truef.write("position\tomega\n")
position = 1
for i in range(numPart):
	for l in range(partLen):
		truef.write(str(position)+'\t'+str(omegas[i])+'\n')
		position+=1
truef.close()


#partitions = prepSimColumns(freq_aln, columns, kappa, omegas, numPart, partLen) # SITE-SPECIFIC
partitions = prepSim(freq_aln, kappa, omegas, numPart, partLen) # GLOBAL FREQS

my_tree, flag_list  = readTree(file="/Users/sjspielman/Omega_MutSel/Simulator/trees/100.tre", show=False) # set True to print out the tree

for n in range(100):
	print n
	
	# sim output
	seqfile = results_dir+"seqs"+str(n)+".fasta"
		
	# Outfile
	ratename = results_dir+'rates_codonfreq'+str(n)+'.txt'
	
	# Simulate
	callSim(partitions, my_tree, seqfile)
	
	# Read in sequences to process
	aln, alnlen, numseq = readAln(seqfile)
	
	# Process
	deriveOmegaAlignment(seqfile, outfile)
