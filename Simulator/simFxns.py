from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *
import re
import random as rn

##### Day 1 of Stephanie Writes a Simulator #####
#################################################

#Alphabetical list of sense codons
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]

# Alphabetical genetic code
gencode = [['GCA', 'GCC', 'GCG', 'GCT'], ['TGC','TGT'], ['GAC', 'GAT'], ['GAA', 'GAG'], ['TTC', 'TTT'], ['GGA', 'GGC', 'GGG', 'GGT'], ['CAC', 'CAT'], ['ATA', 'ATC', 'ATT'], ['AAA', 'AAG'], ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], ['ATG'], ['AAC', 'AAT'], ['CCA', 'CCC', 'CCG', 'CCT'], ['CAA', 'CAG'], ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'] , ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], ['ACA', 'ACC', 'ACG', 'ACT'], ['GTA', 'GTC', 'GTG', 'GTT'], ['TGG'], ['TAC', 'TAT']]

amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N' 'P', 'Q', 'R' , 'S', 'T', 'V', 'W', 'Y']

def generateStateFreqs():
	''' List of randomly-generated reasonable amino acid frequencies.'''
	stateFreqs=[]
	freq=float(1)
	max=0.1
	sum=float(0)
	for i in range(19):
		freq = rn.uniform(0,max)
		while ((freq!=0) & (sum + freq > 1)):
			freq = rn.uniform(0,max)
		sum += freq
		stateFreqs.append(freq)
	stateFreqs.append(1-sum)
	return stateFreqs

# If we want to calculate lambda from the stateFreqs
#def assignFitness(stateFreqs):
# return 0

def generateRoot(stateFreqs, gencode):
	''' Generate a root codon sequence given equilibrium frequencies.'''
	''' Given this root amino acid, select one of its codons at random.''' 
	
	r = rn.uniform(0,1)
	i=0
	sum=stateFreqs[i]
	while sum <= r:
		i+=1
		sum+=stateFreqs[i]	
	codNum = rn.randint(0,len(gencode[i])-1)
	return gencode[i][codNum]
	

stateFreqs=generateStateFreqs()
print stateFreqs
rootSeq = generateRoot(stateFreqs, gencode)
print rootSeq



def sim_subtree(baseSeq, node):
	newSeq = evolve_branch(baseSeq, node) # Provide starting sequence, branch length along which to evolve
	
# Traverse and simulate along tree.	
def traverseTree(baseSeq, tree):
	
	evolve_branch(baseSeq, tree)

	if len(tree.children)>0:
		print len(tree.childre)
		for node in tree.children:
			newSeq = evolve_branch(baseSeq, node)
	
	
	
	
	
	
	print indent, tree.name, tree.BL, num
	num+=1
	if len(tree.children)>0:
		for node in tree.children:
			num=printTree(node, level+1, num)
	return num
































	
	
	
	
	
	
	
	
	
	
	
