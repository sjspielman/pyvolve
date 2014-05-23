## Parser for temp_config.txt. Note that this is almost definitely not a final parser to be released, but for current use to become comfortable with ConfigParser module, etc.

import os
import sys
import re
from string import whitespace
from misc import Model, Genetics
from stateFreqs import *

class parseConfig:
	''' Parse temp_config.txt file '''
	def __init__(self, **kwargs):
		self.configFile = kwargs.get('file', '../temp_config.txt')
		self.options = {} # nested dictionaries of options.
		self.model = Model()
		# Parse the config file into the dictionary, self.options
		confile = open(self.configFile, 'rU')
		parsed = confile.readlines()
		confile.close()
		for line in parsed:
			newline = line.translate(None, whitespace)
			if newline != '':
				self.options[newline.split('=')[0]] = newline.split('=')[1]
		
		# Following two variables are assigned in the function self.configModel() based on the modelClass parameter
		# They may be overwritten based on user specifications in ...
		self.freqType = None
		self.freqBy   = None
		
		# TO DO: will need to check this in some way. 
		self.saveFreq = None
	
		self.molecules = Genetics()
	
	#########################################################################################################################
	################################### CONFIGURE THE EVOLUTIONARY MODEL PARAMETERS #########################################
	#########################################################################################################################

	def configModel(self):
		''' Configure the default evolutionary model specifications.
			self.options['modelClass'] should be either mutsel, codon, amino, nucleotide.
			Depending on which model type, construct default model specifications. 
			Defaults are also be overwritten by any provided parameters .
		'''

		if self.options['modelClass'] == 'codon':
			self.model.params = {'alpha': 1., 'beta': 1., 'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.} }
			self.parseCodonParams()
			self.freqType = 'codon'
			self.freqBy = 'amino'
		
		elif self.options['modelClass'] == 'mutsel':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
			self.parseMuKappa()
			self.freqType = 'codon'
			self.freqBy = 'amino'

		elif self.options['modelClass'] == 'nucleotide':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
			self.parseMuKappa()
			self.freqType = 'nuc'
			self.freqBy = 'nuc'
			
		elif self.options['modelClass'] == 'amino':
			self.model.params = {'myparams': 'fill in later once this class is constructed.'}
			self.freqType = 'amino'
			self.freqBy = 'amino'
			
		else:
			raise AssertionError("You must provide a model class (codon, mutsel, nucleotide, or amino).")




	def parseAminoParam(self):
		''' Write this later etc.'''
		stephanie = 'awesome'


	def parseMuKappa(self):
		''' codon, mutsel, and nucleotide all use mu and kappa. Use this function for all three.'''
	
		if "mu" in self.options:
			try:
				self.options_mu = eval(self.options['mu'])
			except:
				raise ValueError("\nYour mutational parameters do not appear to be numerical.\n")
			for k in self.options_mu:
				try:
					provided = float(self.options_mu[k])
					self.model.params['mu'][k] = provided
				except:
					raise ValueError("\nYou must provide numeric values for any mutational parameters.\n")		
		## NOT YET CLEAR TO ME IF THIS SHOULD BE AN IF OR ELSE!!!!!!!! ARE THESE MUTUALLY EXCLUSIVE? NEED TO THINK ON IT.
		if "kappa" in self.options:
			try:
				kappa = float(self.options['kappa'])
			except:
				raise ValueError("\nIf you wish to provide a kappa (TI/TV) value, it must be numeric.\n")
			self.model.params['mu']['AG'] = self.model.params['mu']['AG'] * kappa
			self.model.params['mu']['GA'] = self.model.params['mu']['GA'] * kappa
			self.model.params['mu']['CT'] = self.model.params['mu']['CT'] * kappa
			self.model.params['mu']['TC'] = self.model.params['mu']['TC'] * kappa
				
				
	def parseCodonParams(self):
		''' Parse the self.options codon model parameters. Use these to overwrite any defaults.
			Allowed input parameters: mu, kappa, omega OR alpha/beta.		
		'''
		
		# Check that reasonable combinations are self.options. Expand on this part later.
		if "alpha" in self.options and "omega" in self.options:
			print "\n\nDANGER, WILL ROBINSON!: You specified value(s) for alpha (dS) and/or beta(dN) in addition to omega (dN/dS). Therefore, I will only pay attention to your self.options omega value.\n\n"
		
		# Grab any user-provided mu and/or kappa parameters
		self.parseMuKappa()
		
		if "omega" in self.options:
			try:
				omega = float(self.options['omega'])
			except:
				raise ValueError("\nIf you wish to provide an omega (equivalent to dN/dS) value, it must be numeric.\n")
			self.model.params['beta'] = omega
			self.model.params['alpha'] = 1.0	
		
		else:
			if "alpha" in self.options:
				try:
					alpha = float(self.options['alpha'])
				except:
					raise ValueError("\nIf you wish to provide an alpha (equivalent to dS) value, it must be numeric.\n")
				self.model.params['alpha'] = alpha
			
			if "beta" in self.options:
				try:
					beta = float(self.options['beta'])
				except:
					raise ValueError("\nIf you wish to provide a beta (equivalent to dN) value, it must be numeric.\n")
				self.model.params['beta'] = beta	
	#########################################################################################################################
	#########################################################################################################################

	#########################################################################################################################
	################################### CONFIGURE THE STATE FREQUENCY OBJECT ################################################

	def configFrequencies(self):
		''' Configure the default evolutionary model specifications.
			self.options['frequencyClass'] should be either equal, random, a user dictionary, or an indication to read from a file (TBD how this will work).
			Construct default frequency specifications (equal) and overwrite as needed.
			NOTE: for codon and mutation-selection models, default codon frequencies are *equal by amino acid*. 
		'''
		
		try:
			self.parseFrequencyClass(self.options['freqClass'])
		except KeyError:
			self.freqObject = EqualFreqs(by = self.freqBy, type = self.freqType)
		
	
	def parseFrequencyClass(self, freqClass):
		''' Define state frequency instance based on user-provided specifications.
			Sanity checking is also performed here.
		'''
		
		# TO DO: self.saveFreqs should be incorporated in some way
		if freqClass == 'equal':
			self.freqObject = EqualFreqs(by = self.freqBy, type = self.freqType)
		elif freqClass == 'random':
			self.freqObject = RandFreqs(by = self.freqBy, type = self.freqType)
		elif type(freqClass) == dict:
			userDict = eval(freqClass)
			userDict = self.sanityUserFreqs(userDict)
			self.freqObject = UserFreqs(by = self.freqBy, type = self.freqType, freqs = userDict)
		else:
			# Check type
			if type(freqClass) == str:
				
				# Here, we are seeing if they have provided a file. If so, sanity check it. Note that they must also provide a format for this file (fasta, phylip, nexus...)
				if os.path.exists(freqClass):
					
					
					
			# could be a file name - check if so here. This should be something dynamic because you can't look for a dictionary (above) and then yell if not a dictionary.
			# DEAL WITH THIS!!
			stephanie = 'really awesome'
			 
	
	def sanityReadFreqs(self, freqClass):
		''' Sanity check the frequency specification if a file was provided.'''
		
	
	
	def sanityReadFreqs_checkAlphabet(self):	
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

	def sanityReadFreqs_checkColumns(self):
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
			
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	def sanityUserFreqs(self, userDict):
		''' Sanity check to ensure that the dictionary provided is compatible with other specifications.'''
		if self.freqBy == 'codon':
			keylen = 3
			self.code = self.molecules.codons
		else:
			keylen = 1
			if self.freqBy == 'nuc':
				self.code = self.molecules.nucleotides
			elif self.freqBy == 'amino':
				self.code = self.molecules.amino_acids	
		sum = 0.
		keysize = len( str(userDict.keys()[0]) ) # Size of first key. All other keys should be the same size as this one. NOTE THAT IF THIS IS REALLY NOT A STRING, IT WILL BE CAUGHT LATER!! Perhaps/definitely this is inelegant, but I'll deal w/ it later.
		for key in userDict:
			assert( type(key) is str), "Your keys must be strings, not any other type (eg int, float, list, etc)."
			assert( len(key) == keysize), "The keys for your frequency dictionary do not have the same length. All keys should be ONE of the following: single letter amino acid symbols, single letter nucleotides, or three-letter codons."
			userDict[key.upper()] = userDict.pop(key) # Ensure upper-case key
			assert ( len(key) == keylen ), ("\n\nThis key,", key, "is an unaccepted format. Please use three-letter codes for codons and one-letter codes for amino acids or nucleotides.")
			assert ( key.upper() in self.code ), ("\n\nThis key,", key, "is not part of the genetic code. Remember, no ambiguities or stop codons allowed.")
			sum += userDict[key.upper()]
		assert ( abs(sum - 1.) < self.zero), ("\n\nIf you provide frequencies, they must sum to 1. The provided frequencies sum to",sum,".")
		return userDict
			
		
	#########################################################################################################################
	#########################################################################################################################
	
		
		
		
		
		
		
					
				
				
			

conf = parseConfig()
conf.configModel()