import os
import re
import numpy as np
import random as rn
from Bio import SeqIO

from misc import Genetics


class StateFreqs(object):
	'''Will return frequencies. '''
	def __init__(self, **kwargs):
		self.type  = kwargs.get('type') # Type of frequencies to RETURN to user. Either amino, codon, nuc.
		self.by    = kwargs.get('by', self.type) # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. If nucleotide, well, yeah.
	
		self.savefile   = kwargs.get('savefile', 'stateFreqs.txt') #default file for saving 
		self.molecules  = Genetics()
		self.aminoFreqs = np.zeros(20)
		self.codonFreqs = np.zeros(61)
		self.nucFreqs   = np.zeros(4)
		self.zero = 1e-10

	def generate(self):
		''' BASE CLASS. NOT IMPLEMENTED. ''' 
	def sanityCheck(self):
		'''BASE CLASS. NOT IMPLEMENTED. ''' 
		
	def setCodeLength(self):
		''' Set the codes and lengths once all, if any, "by" issues are resolved ''' 
		if self.by == 'amino':
			self.code = self.molecules.amino_acids
		elif self.by == 'codon':
			self.code = self.molecules.codons
		elif self.by == 'nuc':
			self.code = self.molecules.nucleotides
		self.length = len(self.code)
			
	def sanityByType(self):
		# Are type and by compatible?
		if self.by == 'nuc' and self.type != 'nuc' and self.type is not None:
			print "CAUTION: If calculations are performed with nucleotides, you can only retrieve nucleotide frequencies."
			print "I'm going to calculate nucleotide frequencies for you."
			self.type = 'nuc'
		if self.by == 'nuc' and self.type == 'amino':
			print "CAUTION: Amino acid frequencies cannot be calculated from nucleotide frequencies."
			print "I'm going to calculate your frequencies using amino acid frequencies."
			self.by = 'amino'
		if self.by == 'amino' and self.type == 'nuc':
			print "CAUTION: Nucleotide frequencies cannot be calculated from amino acid frequencies."
			print "I'm going to calculate nucleotide frequencies for you."
			self.by = 'nuc'

	def amino2codon(self):
		''' Calculate codon frequencies from amino acid frequencies. CAUTION: assumes equal frequencies for synonymous codons!! '''
		count = 0
		for codon in self.molecules.codons:
			ind = self.molecules.amino_acids.index(self.molecules.codon_dict[codon])	
			if codon in self.molecules.genetic_code[ind]:
				numsyn = float(len(self.molecules.genetic_code[ind]))
				self.codonFreqs[count] = self.aminoFreqs[ind]/numsyn
			count += 1
		assert( abs(np.sum(self.codonFreqs) - 1.) < self.zero), "Codon state frequencies improperly generated from amino acid frequencies. Do not sum to 1." 				
				
				
	def codon2amino(self):
		''' Calculate amino acid frequencies from codon frequencies. ''' 
		for a in range(len(self.molecules.amino_acids)):
			codons1 = self.molecules.genetic_code[a]
			for c in codons1:
				ind = self.molecules.codons.index(c)
				self.aminoFreqs[a] += self.codonFreqs[ind]
		assert( abs(np.sum(self.aminoFreqs) - 1.) < self.zero), "Amino acid state frequencies improperly generated from codon frequencies. Do not sum to 1." 
	
	def codon2nuc(self):
		''' Calculate the nucleotide frequencies from the codon frequencies. ''' 
		self.generate() # This will get us the codon frequencies. Now convert those to nucleotide
		self.nucFreqs = np.zeros(4) ## ACGT
		for i in range(61):
			codon_freq = self.codonFreqs[i]
			codon = self.molecules.codons[i]
			for n in range(4):
				nuc =  self.molecules.nucleotides[n]
				nuc_freq = float(codon.count(nuc))/3. # number of that nucleotide in the codon
				if nuc_freq > 0 :
					self.nucFreqs[n] += codon_freq * nuc_freq
		assert( abs(np.sum(self.nucFreqs) - 1.) < self.zero), "Nucleotide state frequencies improperly generated. Do not sum to 1." 


	def assignFreqs(self, freqs):
		''' For generate() functions when frequencies are created generally, assign to a specific type with this function. '''
		if self.by == 'codon':
			self.codonFreqs = freqs
		elif self.by == 'amino':
			self.aminoFreqs = freqs
		elif self.by == 'nuc':
			self.nucFreqs = freqs

	def calcFreqs(self):
		''' Calculate and return state frequencies. type = the type of frequency to return (nuc, codon, amino). '''
		
		# Some sanity checking
		self.sanityCheck()
		assert(self.type is not None), "I don't know what type of frequencies to calculate! I'm quitting."
		self.setCodeLength()
		# This function will generate frequencies for whatever the 'by' is. If the 'type' is different, convert below before returning.
		self.generate() 
		
		if self.type == 'codon':
			if self.by == 'amino':
				self.amino2codon()
			return self.codonFreqs	
		elif self.type == 'amino':
			if self.by == 'codon':
				self.codon2amino()
			return self.aminoFreqs	
		elif self.type == 'nuc':
			if self.by == 'codon':
				self.codon2nuc()
			return self.nucFreqs
		else:
			raise AssertionError("Type of frequencies must be either amino, codon, or nuc.")
			
	######## UNTESTED #########
	def save2file(self, type):
		if type == 'codon':
			np.savetxt(self.savefile, self.codonFreqs)
		elif type == 'amino':
			np.savetxt(self.savefile, self.aminoFreqs)
		elif type == 'nuc':
			np.savetxt(self.savefile, self.nucFreqs)










class EqualFreqs(StateFreqs):
	def __init__(self, 	**kwargs):
		super(EqualFreqs, self).__init__(**kwargs)

	def sanityCheck(self):
		''' Check that internals are all ok before performing calculations. ''' 
		self.sanityByType()

	def generate(self):
		freqs = np.array(np.repeat(1./float(self.length), self.length))
		assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		self.assignFreqs(freqs)
					
					
					
					
					
					
class RandFreqs(StateFreqs):
	def __init__(self, **kwargs):
		super(RandFreqs, self).__init__(**kwargs)

	def sanityCheck(self):
		''' Check that internals are all ok before performing calculations. ''' 
		self.sanityByType()
		# Notify users if they have strange options, but proceed anyways.
		if self.type != self.by and self.type is not None:
			print "You have specified", self.type, "random state frequencies to be calculated as random", self.by, "frequencies, and then converted."
			print "This is a strange choice, but I will proceed anyways."
	
	def generate(self):
		freqs = np.zeros(self.length)
		max = 2. / (self.length) # times 2 for ease/speed/slightly less ridiculousness.
		sum = 0.
		for i in range(int(self.length) - 1):
			freq = rn.uniform(0,max)
			while ((freq==0) or (sum + freq > 1)):
				freq = rn.uniform(0,max)
			sum += freq
			freqs[i] = freq
		freqs[-1] = (1.-sum)	
		assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly generated. Do not sum to 1." 
		self.assignFreqs(freqs)
		
	
	
	

class UserFreqs(StateFreqs):
	''' Assign frequencies based on user input. Assume that if not specified, the frequency is zero. 
		For this subclass, 'by' should correspond to the type that they've entered. 'type' should correspond to what they want.
		If 'by' is not provided, I can try to guess what they want, but I won't always succeed.	
	'''
	def __init__(self, **kwargs):
		super(UserFreqs, self).__init__(**kwargs)	
		self.givenFreqs = kwargs.get('freqs', {}) # Dictionary of desired frequencies.	
	
	
	
	#self.type  = kwargs.get('type') # Type of frequencies to RETURN to user. Either amino, codon, nuc.
	#self.by    = kwargs.get('by', self.type) # Type of frequencies to base generation on. If amino, get amino acid freqs and convert to codon freqs, with all synonymous having same frequency. If codon, simply calculate codon frequencies independent of their amino acid. If nucleotide, well, yeah.
	
	def guessBy(self):
		''' If neither by='' nor type='' is not provided, we can use the keys to try and ascertain. Note that we will not always be able to (e.g. if frequencies are provided for ACGT, we're in trouble.)'''
		keys = self.givenFreqs.keys() 
		assert(type(keys[0]) is str), "Your keys must be strings, not any other type (eg int, float, list, etc."
		keysize = len(keys[0]) # Size of first key. All other keys should be the same size as this one.
		# For keysize of 3.
		if keysize == 3:
			if self.by is None:
				print "I think you have codon data!"
				self.by = 'codon'
			elif self.by != 'codon':
				raise AssertionError("You told me you were providing codon frequencies, but your dictionary keys aren't codons. I'm quitting!")
		# If the keysize is 1, we need to guess if nucleotide or amino acid.
		elif keysize == 1:
			notNuc=False
			for key in keys:
				# All nucleotides are also aa codes, so can just check this list.
				assert(key in self.molecules.amino_acids), "Your keys don't correspond to genetics. I'm quitting!"
				if key not in self.molecules.nucleotides:
					notNuc = True
			if notNuc:
				print "I think you have amino acid data!"
				self.by = 'amino'
			else:
				raise AssertionError("It is unclear if you have amino acid or nucleotide data. Please specify using the 'by=' flag.")
		else:
			raise AssertionError("Your keys need to be of length 1 (amino acids and nucleotides) or 3 (codons).")	


	def sanityCheck(self):
		''' Perform sanity checks on user-provided frequencies. ''' 
	
		# Did the user actually provide a dictionary?
		assert(type(self.givenFreqs) is dict), "You must provide a dictionary of frequencies. Keys should be single letter amino acid symbols, single letter nucleotides, or three-letter codons."
		# Try to guess the type of frequencies if they didn't specify. This can happen only if ACGT not provided alone, because this (although unlikely!!!!) could be amino acid data.
		if self.by is None:
			print "I'm going to try and guess your data"
			self.guessBy()
		
		# Is the dictionary correct given the type of frequencies they are providing?
		self.setCodeLength()
		if self.by == 'codon':
			keylen = 3
		else:
			keylen = 1
		sum = 0
		keysize = len( str(self.givenFreqs.keys()[0]) ) # Size of first key. All other keys should be the same size as this one. NOTE THAT IF THIS IS REALLY NOT A STRING, IT WILL BE CAUGHT LATER!! Perhaps/definitely this is inelegant, but I'll deal w/ it later.
		for key in self.givenFreqs:
			assert( type(key) is str), "Your keys must be strings, not any other type (eg int, float, list, etc."
			assert( len(key) == keysize), "The keys for your frequency dictionary do not have the same length. All keys should be ONE of the following: single letter amino acid symbols, single letter nucleotides, or three-letter codons."
			self.givenFreqs[key.upper()] = self.givenFreqs.pop(key) # Ensure upper-case key
			assert ( len(key) == keylen ), ("\n\nThis key,", key, "is an unaccepted format. Please use three-letter codes for codons and one-letter codes for amino acids or nucleotides.")
			assert ( key.upper() in self.code ), ("\n\nThis key,", key, "is not part of the genetic code. Remember, no ambiguities or stop codons allowed.")
			sum += self.givenFreqs[key.upper()]
		assert ( abs(sum - 1.) < self.zero), ("\n\nIf you provide frequencies, they must sum to 1. The provided frequencies sum to",sum,".")
		
		if self.type is None:
			print "You did not specify the final type of frequencies to return. I'm going to return exactly what you provided."
			self.type = self.by
		# Can now check by/type compatibility now that by and type are assigned.
		self.sanityByType()
	
	def generate(self):
		freqs = np.zeros(self.length)
		for i in range(self.length):
			element = self.code[i]
			if element in self.givenFreqs:
			 	freqs[i] = self.givenFreqs[element]
		assert( abs(np.sum(freqs) - 1.) < self.zero), "State frequencies improperly converted from provided to internal object. Do not sum to 1." 
		self.assignFreqs(freqs)
		









########## NEEDS EXTRAORDINARY AMOUNTS OF DEBUGGING ############## 5/10/14.
class ReadFreqs(StateFreqs):
	''' Retrieve frequencies from a file. Can either do global or specify a particular column/group of columns ''' 
	def __init__(self, **kwargs):
		super(ReadFreqs, self).__init__(**kwargs)
		self.seqfile     = kwargs.get('file', None)   # Can also read frequencies from a sequence file
		self.format      = kwargs.get('format', 'fasta') # Default for that file is fasta
		self.whichCol    = kwargs.get('columns', None)     # Which columns we are collecting frequencies from. Default is all columns combined. IF YOU GIVE IT A NUMBER, INDEX AT 0!!!!
	
	def setUpSeqs(self):
		''' Set up sequences and relevent variables for frequency collection. '''
		self.seqs = []
		self.numseq = len(self.rawrecords)
		self.alnlen = len(self.rawrecords[0]) # This will only come into play if we're collecting columns.
		for entry in self.rawrecords:
			self.seqs.append(str(entry.seq))
		print self.seqs
			
	def checkAlphabet(self):	
		''' Ensure sequences in file provided are the correct alphabet. ''' 
		fullSeq = "".join(self.seqs)
		fullSeq = fullSeq.translate(None, '-?.*BJOUXZ') # Remove ambiguities
		fullLength = len(fullSeq)
		DNA = re.compile(r"[ACGT]")
		PROT = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]")
		# Double check that everything in this culled dataset is a correct character given intended alphabets
		dnaChar = len(DNA.findall(fullSeq))
		protChar = len(PROT.findall(fullSeq))
		#print dnaChar, protChar, fullLength
		if self.by == 'amino':
			assert(dnaChar != protChar), "Are you sure this is a protein sequence file?"
			assert(protChar == fullLength), "Your sequences do not appear to be the correct alphabet for your frequency calculations, or you have bizarre characters in the sequences."
		elif self.by == 'codon' or self.by == 'nuc':
			assert(dnaChar == fullLength), "Your sequences do not appear to be the correct alphabet for your frequency calculations, or you have bizarre characters in the sequences."
	
	def checkColumns(self):
		''' Sanity checking for proper column specification, ONLY if user has requested frequencies come from columns.
			1. Columns must be a list. We can work with a single int or tuple, but that's it.
			2. Sequences must be an alignment if columns are requested.
			3. There actually need to be that many columns in the alignment. As in, if user requests column #20, there better be a column #20.
		'''
		# Ensure self.whichCol is a list.
		if type(self.whichCol) is not list:
			if type(self.whichCol) is int:
				self.whichCol = [self.whichCol]
			if type(self.whichCol) is tuple:
				self.whichCol = list(self.whichCol)
			else:
				raise AssertionError("If you'd like to read frequencies by column, you must provide a *list* of which columns (indexed at 0) you want.")	
		
		# Ensure sequences were an alignment
		if self.whichCol is not None:
			for i in range(1, self.numseq):
				seqlen = len(self.rawrecords[i].seq)
				assert(seqlen == self.alnlen), "Your provided file does not appear to be an alignment. If you want to collect frequency information from columns, it *must* be an alignment!!"	
		
		# Ensure that indexing is within range. 
		assert(max(self.whichCol) < self.alnlen), "Your provided alignment does not contain the columns you requested."	
			
			
			
	def sanityCheck(self):
		''' Sanity checking for ReadFreqs.''' 
		# Make sure the file we need to read from exists and has correct format
		assert (os.path.exists(self.seqfile)), "The file you would like to use for state frequencies doesn't exist."
		try:
			self.rawrecords = list(SeqIO.parse(self.seqfile, self.format))
			self.setUpSeqs()
		except:
			raise AssertionError("Your sequence/alignment file is not formatted the way you specified (as", str(self.format),"). Note that default format is FASTA.")
		
		# Check columns, as needed
		if self.whichCol is not None:
			self.checkColumns()
		
		# Make sure alphabet is correct, given the 'by'
		if self.by == 'codon':
			for seq in self.seqs:
				assert( len(self.seqs) % 3 == 0), "Your file does not appear to contain codon sequences (sequence lengths are not all multiples of three)."
		self.checkAlphabet()
		
		
		# Good to go!
		self.sanityByType()
		self.setCodeLength()
		
		
	'''
	def getSeq(self):
		'''''' If we want columns, we must get a string of the specific columns we're collecting from.
			Otherwise, we can just turn the whole alignment into a single string.
		''''''
		seq = ''
		if self.by == "codon":
			print "i'm a codon"

		elif self.by == "amino":
			if self.whichCol:
				for col in self.whichCol:
					for row in self.seqs:
						seq += row[col]
			else:
				for entry in self.seqs:
					seq += entry
			#Remove ambig, nonstandard, gaps
			seq = seq.upper()
			seq = seq.translate(None, '-?.*BJOUXZ')
		return seq
		

	def generate(self):
	
		seq = self.getSeq()	
		
		if self.by == 'codon':
			for i in range(0, len(seq),3):
				codon = seq[i:i+3]
				## check for stop codons.
				ind = self.code.index(codon)
				self.codonFreqs[ind]+=1
			self.codonFreqs = np.divide(self.codonFreqs, len(seq)/3)
		
		elif self.by == 'amino':
			for i in range(0, len(seq)):
				ind = self.code.index(seq[i])
				self.aminoFreqs[ind]+=1
			self.aminoFreqs = np.divide(self.aminoFreqs, len(seq))		
		
	'''
		
		
		
		
		
		
		
		
		
		
		
		
		
		
