from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *
import re
import random as rn
import numpy as np


#Alphabetical list of sense codons
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
#
# Alphabetical genetic code
gencode = [['GCA', 'GCC', 'GCG', 'GCT'], ['TGC','TGT'], ['GAC', 'GAT'], ['GAA', 'GAG'], ['TTC', 'TTT'], ['GGA', 'GGC', 'GGG', 'GGT'], ['CAC', 'CAT'], ['ATA', 'ATC', 'ATT'], ['AAA', 'AAG'], ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], ['ATG'], ['AAC', 'AAT'], ['CCA', 'CCC', 'CCG', 'CCT'], ['CAA', 'CAG'], ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'] , ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], ['ACA', 'ACC', 'ACG', 'ACT'], ['GTA', 'GTC', 'GTG', 'GTT'], ['TGG'], ['TAC', 'TAT']]

amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N' 'P', 'Q', 'R' , 'S', 'T', 'V', 'W', 'Y']


slist=[['AAG'], ['AAT'], ['AAA'], ['AAC'], ['ACC', 'ACT', 'ACG'], ['ACA', 'ACT', 'ACG'], ['ACA', 'ACC', 'ACT'], ['ACA', 'ACC', 'ACG'], ['CGA', 'AGG'], ['AGT'], ['CGG', 'AGA'], ['AGC'], ['ATC', 'ATT'], ['ATA', 'ATT'], [], ['ATA', 'ATC'], ['CAG'], ['CAT'], ['CAA'], ['CAC'], ['CCC', 'CCT', 'CCG'], ['CCA', 'CCT', 'CCG'], ['CCA', 'CCC', 'CCT'], ['CCA', 'CCC', 'CCG'], ['AGA', 'CGC', 'CGT', 'CGG'], ['CGA', 'CGT', 'CGG'], ['AGG', 'CGA', 'CGC', 'CGT'], ['CGA', 'CGC', 'CGG'], ['TTA', 'CTC', 'CTT', 'CTG'], ['CTA', 'CTT', 'CTG'], ['TTG', 'CTA', 'CTC', 'CTT'], ['CTA', 'CTC', 'CTG'], ['GAG'], ['GAT'], ['GAA'], ['GAC'], ['GCC', 'GCT', 'GCG'], ['GCA', 'GCT', 'GCG'], ['GCA', 'GCC', 'GCT'], ['GCA', 'GCC', 'GCG'], ['GGC', 'GGT', 'GGG'], ['GGA', 'GGT', 'GGG'], ['GGA', 'GGC', 'GGT'], ['GGA', 'GGC', 'GGG'], ['GTC', 'GTT', 'GTG'], ['GTA', 'GTT', 'GTG'], ['GTA', 'GTC', 'GTT'], ['GTA', 'GTC', 'GTG'], ['TAT'], ['TAC'], ['TCC', 'TCT', 'TCG'], ['TCA', 'TCT', 'TCG'], ['TCA', 'TCC', 'TCT'], ['TCA', 'TCC', 'TCG'], ['TGT'], [], ['TGC'], ['CTA', 'TTG'], ['TTT'], ['CTG', 'TTA'], ['TTC']]

def generateStateFreqs():
	''' List of randomly-generated reasonable codon frequencies.'''
	aaFreqs=[]
	freq=float(1)
	max=0.1
	sum=float(0)
	for i in range(19):
		freq = rn.uniform(0,max)
		while ((freq!=0) & (sum + freq > 1)):
			freq = rn.uniform(0,max)
		sum += freq
		aaFreqs.append(freq)
	aaFreqs.append(1-sum)
	print aaFreqs
	# Now we need to convert this to codons by dividing based on the number of synonymous codons that amino acid has
	stateFreqs=[]
	for i in range(20):		
		num = len(gencode[i])
		new = aaFreqs[i] / num
		for j in range(num):
			stateFreqs.append(new)
		print new	
	#print stateFreqs
	
	
	

def generateCodon(freqList, gencode):
	''' Generate a codon sequence given a set of probabilities. For generating a root sequence, these are equilibrium state frequencies. For an evolving sequence, these are the transition probabilities calc'd for that branch.'''
	''' Given this root amino acid, select one of its codons at random as all codons for a given amino acid all currently have the same frequencies (synonymous changes are therefore neutral).''' 	
	r = rn.uniform(0,1)
	i=0
	sum=stateFreqs[i]
	while sum <= r:
		i+=1
		sum+=stateFreqs[i]	
	codNum = rn.randint(0,len(gencode[i])-1)
	
	return gencode[i][codNum]
	
generateStateFreqs()




































	
	
	
	
	
	
	
	
	
	
	
