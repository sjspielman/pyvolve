import os, re, subprocess, shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *
import numpy as np


#################################################################################################################################
def ensure_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
    return 0
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
### Parse indelible rates
def parseIndelible(results_dir, truerates, ratefile, n):

	indelible_rate = open(ratefile, 'r')
	truelines = indelible_rate.readlines()
	indelible_rate.close()
	
	truelines=truelines[10:] ## only keep these lines since before that it's all header crap.	
	
	outfile=results_dir+'truerates'+str(n)+'.txt'
	outrates=open(outfile, 'w')
	outrates.write("position\tomega\n")
	for line in truelines:
		
		find=re.search('^(\d+)\t(\d+)\t', line)
		if find:
			position=find.group(1)
			rate=int(find.group(2))
			outrates.write(position+'\t'+str(truerates[rate])+'\n')
	outrates.close()
				
	return 0	
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
for codon in codons:
	(syn, nonsyn) = findSynNonsyn(codon)
	slist.append(syn)
	nslist.append(nonsyn)

#################################################################################################################################
#################################################################################################################################	


## Conduct 100 simulations to confirm that the correlation between ML dN/dS and our derived dN/dS really exists

# These are the rate categories according to which we simulate on 1/31/14.
# 50 equally spaced omega categories from [0.5,0.95]
truerates = [0.05000000, 0.06836735, 0.08673469, 0.10510204, 0.12346939, 0.14183673, 0.16020408, 0.17857143, 0.19693878, 0.21530612, 0.23367347, 0.25204082, 0.27040816, 0.28877551, 0.30714286, 0.32551020, 0.34387755, 0.36224490, 0.38061224, 0.39897959, 0.41734694, 0.43571429, 0.45408163, 0.47244898, 0.49081633, 0.50918367, 0.52755102, 0.54591837, 0.56428571, 0.58265306, 0.60102041, 0.61938776, 0.63775510, 0.65612245, 0.67448980, 0.69285714, 0.71122449, 0.72959184, 0.74795918, 0.76632653, 0.78469388, 0.80306122, 0.82142857, 0.83979592, 0.85816327, 0.87653061, 0.89489796, 0.91326531, 0.93163265, 0.95]


date='2.2.14'
home='/Users/sjspielman/' # Change if on MacMini or MacBook

sim_dir=home+'Dropbox/MutSelProject/quickCalc/SimMaterials/'
results_dir=home+'Dropbox/MutSelProject/quickCalc/SimSeqs_'+date+'/'
ensure_dir(results_dir)

os.chdir(sim_dir) # Stay here so all the files Indelible generates stay there.
for n in range(50):
	print n
	
	# Simulate, get true dN/dS, save alignment and truerates file
	simCall = 'indelible control.txt'
	subprocess.call(simCall, shell=True)
	
	parseIndelible(results_dir, truerates, 'results_RATES.txt', n)
	newAln = results_dir+'aln'+str(n)+'.fasta'
	shutil.copy('results.fas', newAln)
	
	# Calculate derived dN/dS from the alignment. Save those values to rates_codonfreq(n).txt
	# aln will contain the alignment. 
	aln=[] #list of lists wherein each nested list is a row
	aln_raw = list(SeqIO.parse(newAln, 'fasta'))
	for record in aln_raw:
		aln.append(str(record.seq))
	
	alnlen=len(aln[0])
	numseq=len(aln)

	ratename = results_dir+'rates_codonfreq'+str(n)+'.txt'
	ratefile=open(ratename, 'w')
	ratefile.write("position\tomega_simple\tomega_count\n")
	
	position=1
	for col in range(0,alnlen,3):
	
		kN=0 #dN numerator
		nN_simple=0 #dN denominator. Does not consider number of nonsyn options
		nN_count=0 #dN denominator. DOES consider number of nonsyn options.
		
		fix_sum=0
		
		codonFreq=np.zeros(len(codons)) # will contain frequencies for all codons in a given column
		nonZero=[] # will contain the nonzero indices for codonFreq
		
		#Find codon frequencies
		for row in range(numseq):
			codon=aln[row][col:col+3]
			codonFreq[codons.index(codon)] += 1	
		codonFreq/=float(numseq)

		# Fill nonZero with codonFreq indices whose values are not 0
		for i in range(len(codonFreq)):
			if codonFreq[i] > 0:
				nonZero.append(i)
	
		# Calculations
		for i in nonZero:
			fix_sum=0
			
			### Nonsynonymous.  and BH methods here
			for nscodon in nslist[i]:
				nscodon_freq = codonFreq[codons.index(nscodon)]
				if nscodon_freq==0 or codonFreq[i]==nscodon_freq:
					continue
				else:
					fix_sum += fix(float(codonFreq[i]), float(nscodon_freq))					
					nN_simple += codonFreq[i]
					nN_count += codonFreq[i] * len(nslist[i])
			kN += fix_sum*codonFreq[i]
	
		# Final dN/dS
		dNdS_simple=kN/nN_simple
		dNdS_count=kN/nN_count
		
		ratefile.write(str(position)+'\t'+str(dNdS_simple)+'\t'+str(dNdS_count)+'\n')
		position+=1
			
	ratefile.close()