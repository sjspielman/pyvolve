## Parser for temp_config.txt. Note that this is almost definitely not a final parser to be released, but for current use to become comfortable with ConfigParser module, etc.

import os
import sys
import re
from Bio import SeqIO
from string import whitespace
from misc import Model, Genetics
from stateFreqs import *

class parseConfig:
	''' Parse temp_config.txt file '''
	def __init__(self, **kwargs):
		self.configFile = kwargs.get('file', '../temp_config.txt')
		self.molecules = Genetics()
		self.model = Model()   # defined in self.configModelParams
		self.freqObject = None # defined in self.configModelFrequencies
		self.zero = 1e-10
		
		# Nested dictionary of options. Will all go to None, by default, here, and any provided options will be overridden.
		####### THESE OPTIONS ARE HARD-CODED THROUGHOUT THIS SCRIPT. THEY CANNOT BE CHANGED. ########
		self.options = {'treefile': None, 
						'modelClass': None, 
						'equilibriumFrequencies': None, 
						'freqBy':None, 
						'saveFrequencies': None,
						'freqConstraint': None,
						'mu':None, 
						'kappa': None, 
						'omega': None, 
						'alpha': None, 
						'beta': None, 
						'numberPartitions': None, 
						'partitionSize': None, 
						'outpath': None, 
						'outfile': None
					   }
	
		# Parse the config file into the dictionary, self.options
		confile = open(self.configFile, 'rU')
		parsed = confile.readlines()
		confile.close()
		for line in parsed:
			newline = line.translate(None, whitespace)
			if newline != '':
				self.options[newline.split('=')[0]] = newline.split('=')[1]



	
	
	def configModel(self):
		''' Crux function. Will generate the full model object (incl. params and eqfreqs). 
			Calls other crux-like functions:
				1. configModelParams
				2. configModelFrequencies
				etc?
		'''

		self.configModelParams()
		self.configModelFrequencies()
		#print self.model.params
		#print self.freqObject.by
		
	
	#########################################################################################################################
	################################### CONFIGURE THE EVOLUTIONARY MODEL PARAMETERS #########################################
	#########################################################################################################################

	def configModelParams(self):
		''' Configure the default evolutionary model specifications.
			self.options['modelClass'] should be either mutsel, codon, amino, nucleotide.
			Depending on which model type, construct default model specifications. 
			Defaults are also be overwritten by any provided parameters .
		'''

		if self.options['modelClass'] == 'codon':
			self.model.params = {'alpha': 1., 'beta': 1., 'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.} }
			self.parseCodonParams()
			self.freqType = 'codon'
		
		elif self.options['modelClass'] == 'mutsel':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
			self.parseMuKappa()
			self.freqType = 'codon'

		elif self.options['modelClass'] == 'nucleotide':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
			self.parseMuKappa()
			self.freqType = 'nuc'
			
		elif self.options['modelClass'] == 'amino':
			self.model.params = {'myparams': 'fill in later once this class is constructed.'}
			self.freqType = 'amino'
			
		else:
			raise AssertionError("You must provide a model class (codon, mutsel, nucleotide, or amino).")

		## Configure self.freqBy (how the frequencies are calculated). Must be consistent w/ specified model.
		self.configFreqBy()
		
		
		
		
		
		
	def parseAminoParam(self):
		''' Write this later etc.'''
		stephanie = 'awesome'


	def parseMuKappa(self):
		''' codon, mutsel, and nucleotide all use mu and kappa. Use this function for all three.'''
	
		if self.options['mu'] is not None:
			try:
				self.options_mu = eval(self.options['mu'])
			except:
				raise ValueError("\nYour mutational parameters do not appear to be numerical.\n")
			for key in self.options_mu:
				try:
					provided = float(self.options_mu[key])
					self.model.params['mu'][key] = provided
				except:
					raise ValueError("\nYou must provide numeric values for any mutational parameters.\n")		
		## NOT YET CLEAR TO ME IF THIS SHOULD BE AN IF OR ELSE!!!!!!!! ARE THESE MUTUALLY EXCLUSIVE? NEED TO THINK ON IT.
		if self.options['kappa'] is not None:
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
		if self.options['omega'] is not None and (self.options['alpha'] is not None or self.options['beta'] is not None):
			print "\n\nDANGER, WILL ROBINSON!: You specified value(s) for alpha (dS) and/or beta (dN) in addition to omega (dN/dS). Therefore, I will only pay attention to your specified omega value.\n\n"
		
		# Grab any user-provided mu and/or kappa parameters
		self.parseMuKappa()
		
		if self.options['omega'] is not None:
			try:
				omega = float(self.options['omega'])
			except:
				raise ValueError("\nIf you wish to provide an omega (equivalent to dN/dS) value, it must be numeric.\n")
			self.model.params['beta'] = omega
			self.model.params['alpha'] = 1.0	
		
		else:
			if self.options['alpha'] is not None:
				try:
					alpha = float(self.options['alpha'])
				except:
					raise ValueError("\nIf you wish to provide an alpha (equivalent to dS) value, it must be numeric.\n")
				self.model.params['alpha'] = alpha
			
			if self.options['beta'] is not None:
				try:
					beta = float(self.options['beta'])
				except:
					raise ValueError("\nIf you wish to provide a beta (equivalent to dN) value, it must be numeric.\n")
				self.model.params['beta'] = beta	
	
	def configFreqBy(self):
		''' Ensure that freqBy and model are fully compatible.'''
		m = self.options['modelClass'] # just for readability
		self.freqBy = self.options['freqBy']
		
		# If not specified, assign default
		if self.freqBy is None:
			if self.options['modelClass'] != 'nucleotide':
				self.freqBy = 'amino'
			else:
				self.freqBy = 'nuc'	
		# If specified, ensure it's ok
		else:
			if (m == 'codon' or m == 'mutsel' or m == 'amino') and self.freqBy == 'nuc':
				print "\nCAUTION: Nucleotide frequencies cannot be used for your specified evolutionary model class,", m, ". By default, amino acids will be used for calculations.\n"
				self.freqBy = 'amino'
			elif m == 'nucleotide' and self.freqBy != 'nuc':
				print "\nCAUTION: For nucleotide models, only nucleotide frequencies can be used in internal calculations. I'm defaulting to that setting.\n"
				self.freqBy = 'nuc'

	#########################################################################################################################
	#########################################################################################################################










	#########################################################################################################################
	################################### CONFIGURE THE STATE FREQUENCY OBJECT ################################################

	def configModelFrequencies(self):
		''' Configure the default evolutionary model specifications.
			self.options['frequencyClass'] should be either equal, random, a user dictionary, or an indication to read from a file (TBD how this will work).
			Construct default frequency specifications (equal) and overwrite as needed.
			NOTE: for codon and mutation-selection models, default codon frequencies are *equal by amino acid*. 
		'''
		
		## TO DO: incorporate functionality for savefile, other ....
		#try:
		self.parseFrequencyClass( freqClass = self.options['equilibriumFrequencies'] )
		#except:
		#	print "\n CAUTION: I don't understand your specification for equilibriumFrequencies. So, I'm going to default to using equal frequencies.\n"
		#	self.freqObject = EqualFreqs(by = self.freqBy, type = self.freqType)
		
	
	def parseFrequencyClass(self, **kwargs):
		''' Define state frequency instance based on user-provided specifications.
			Sanity checking is also performed here.
			TO DO: SAVEFREQS FILE SHOULD BE INCORPORATED IN SOME WAY
		'''
		freqClass = kwargs.get('freqClass')
		saveFile  = kwargs.get('saveFile', None)

		# Equal frequencies
		if freqClass.lower() == 'equal':
			self.freqObject = EqualFreqs(by = self.freqBy, type = self.freqType)
		
		# Random frequencies
		elif freqClass.lower() == 'random':
			self.freqObject = RandFreqs(by = self.freqBy, type = self.freqType)
			
		# Read frequencies
		elif os.path.exists(freqClass):
			# Check if columns were also specified
			if self.options['columns']:
				userColumns = self.options['columns']
			else:
				userColumns = None
			seqfile_ext, final_columns = self.sanityReadFreqs(freqClass, userColumns)
			self.freqObject = ReadFreqs(by = self.freqBy, type = self.freqType, file = freqClass, format = seqfile_ext, columns = final_columns)
			if self.options['freqConstraint']:
				constraint = self.sanityReadUserFreqs_checkConstraint()
				self.freqObject.constraint = constraint
		
		# User frequencies, as last parsing attempt
		else:
			try:
				userDict = eval(freqClass)
			except:
				raise AssertionError() # Will ultimately result in using default
			if type(userDict) == dict:
				userDict = self.sanityUserFreqs(userDict)
				self.freqObject = UserFreqs(by = self.freqBy, type = self.freqType, freqs = userDict)
				if self.options['freqConstraint']:
					constraint = self.sanityReadUserFreqs_checkConstraint()
					self.freqObject.constraint = constraint
			else:
				raise AssertionError() # Will ultimately result in using default

	
	
		####### TO DO ##########
		# Any other frequency options?? For instance, savefile????????
		if saveFile:
			self.freqObject.savefile = saveFile

	
	
	def sanityReadUserFreqs_checkConstraint(self):
		'''If a constraint is specified, must be float between 0-1.'''
		try:
			constraint = float(self.options['freqConstraint'])
		except ValueError:
			print "To specify a constraint on equilibrium frequencies, provide a decimal in the range (0,1]."
		assert(constraint > 0. and constraint <= 1.), "\n\nTo specify a constraint on equilibrium frequencies, provide a decimal in the range (0,1]."
		return constraint
		
		
		
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
		keysize = len( str(userDict.keys()[0]) ) # Size of first key. All other keys should be the same size as this one. NOTE THAT IF THIS IS REALLY NOT A STRING, IT WILL BE CAUGHT LATER!! Perhaps/definitely this is inelegant, but I'll deal w/ it later.
		sum = 0.
		for key in userDict:
			assert( type(key) is str), "\n\nYour keys must be strings, not any other type (eg int, float, list, etc)."
			assert( len(key) == keysize), "\n\nThe keys for your frequency dictionary do not have the same length. All keys should be ONE of the following: single letter amino acid symbols, single letter nucleotides, or three-letter codons."
			userDict[key.upper()] = userDict.pop(key) # Ensure upper-case key
			assert ( len(key) == keylen ), ("\n\nThis key,", key, "is an unaccepted format. Please use three-letter codes for codons and one-letter codes for amino acids or nucleotides.")
			assert ( key.upper() in self.code ), ("\n\nThis key,", key, "is not part of the genetic code. Remember, no ambiguities or stop codons allowed.")
			sum += userDict[key.upper()]
		assert ( abs(sum - 1.) < self.zero), ("\n\nIf you provide frequencies, they must sum to 1. The provided frequencies sum to",sum,".")
		return userDict
	
	
	
	
	
	
	
	
	
	
	
	
	
	def sanityReadFreqs(self, seqfile, userColumns):
		''' Sanity check the frequency specification if a file was provided.
			Return a proper column list, if exists, and sequence file extension.		
		'''
		
		# Ensure that the file can be parsed. 
		seqs, ext = self.sanityReadFreqs_readSeqFile(seqfile)

		# Ensure correct alphabet
		self.sanityReadFreqs_checkAlphabet(seqs)
		
		# If columns provided, ensure correct format
		if userColumns:
			columns = self.sanityReadFreqs_checkColumns(seqs, userColumns)
	
		return (ext, userColumns)
	
	
	
	def sanityReadFreqs_readSeqFile(self, seqfile):
		''' Read in a sequence file. Also return the format the file is in, if we can parse it.'''
		# Accepted formats are fasta, phylip, clustal, nexus, stockholm. 
		exts = ['fasta', 'phylip', 'clustal', 'nexus', 'stockholm']
		for ext in exts:
			try:
				records = list(SeqIO.parse(file, ext))
			except:
				pass
			if records:
				break
		if records is None:
			raise AssertionError("\n\n I'm having trouble reading the file from which I need to calculated equilibrium frequencies.\n Please specify a file in either fasta, phylip, nexus, clustal, or stockholm format.")
		seqs = []
		for record in records:
			seqs.append(str(record.seq))
		return (seqs, ext)
	
	
	
	
	
	def sanityReadFreqs_checkAlphabet(self, seqs):	
		''' Ensure sequences in file provided are the correct alphabet. ''' 
		fullSeq = "".join(seqs)
		fullSeq = fullSeq.translate(None, '-?.*BJOUXZ') # Remove ambiguities
		fullLength = len(fullSeq)
		DNA = re.compile(r"[ACGT]")
		PROT = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]")
		dnaChar = len(DNA.findall(fullSeq))
		protChar = len(PROT.findall(fullSeq))
		if self.freqBy == 'amino':
			assert(dnaChar != protChar), "\n\nAre you sure this is a protein sequence file? I'm quitting."
			assert(protChar == fullLength), "\n\nYour sequences do not appear to be the correct alphabet for your frequency calculations, or you have bizarre characters in the sequences. I'm quitting."
		elif self.freqBy == 'codon' or self.freqBy == 'nuc':
			assert(dnaChar == fullLength), "\n\nYour sequences do not appear to be the correct alphabet for your frequency calculations, or you have bizarre characters in the sequences. I'm quitting."




	def sanityReadFreqs_checkColumns(self, seqs, userColumns):
		''' Sanity checking for proper column specification, ONLY if user has requested frequencies come from columns.
			1. Columns must be a list. We can work with a single int or tuple, but that's it.
			2. Sequences must be an alignment if columns are requested.
			3. There actually need to be that many columns in the alignment. As in, if user requests column #20, there better be a column #20.
		'''
		# Ensure columns is a list.
		if type(userColumns) != list:
			if type(userColumns) == tuple:
				columns = list(userColumns)
			elif type(userColumns) == str:
				columns = eval(userColumns)
			else:	
				raise AssertionError("If you'd like to read frequencies by column, you must provide a list of which columns (indexed at 0) you want.")	
		else:
			columns = userColumns
		
		# Ensure columns contains only integers
		for col in columns:	
			assert(type(col) == int), "\n\nYou must provide integer values for any columns you want to specify when collecting equilibrium frequencies."
	
		# Ensure sequences were an alignment
		alnlen = len(seqs[0])
		for i in range(1, len(seqs)):
			seqlen = len(seqs[i])
			assert(seqlen == alnlen), "\n\nIf you want to collect frequency information from columns, your specified sequence file *must* be an alignment. Currently, it is not.\n"	
		
		# Ensure that indexing is within range. 
		assert(min(columns) >= 0 and max(columns) < alnlen), "\n\nYour provided alignment does not contain the columns you requested. Remember that you should index at 0."	
			
			

conf = parseConfig()
conf.configModel()