## Parser for temp_config.txt. Note that this is almost definitely not a final parser to be released, but for current use to become comfortable with ConfigParser module, etc.

import ConfigParser
import os
import sys
from misc import Model

class parseConfig:
	''' Parse temp_config.txt file '''
	def __init__(self, **kwargs):
		self.configFile = kwargs.get('file', '../temp_config.txt')
		print self.configFile
		self.parser = ConfigParser.ConfigParser()	
		self.parser.read(self.configFile)
		self.options = {} # nested dictionaries of options.
		self.model = Model()

	def parseSections(self):
		''' Organize config specifications into a nested dictionary of options.'''
		for section in self.parser.sections():
			self.options[section] = {}
			for (k,v) in self.parser.items(section):
				self.options[section][k] = v
	
	def configModel(self):
		''' Configure the default evolutionary model specifications.
			options['MODEL'] should contain either mutsel, codon, amino, nucleotide.
			Depending on which model type, construct default model specifications. 
			These will automatically be overwritten by anything different provided in options['PARAMETERS'].
		'''
		modelClass = self.options['MODEL']['modelclass'] # this is a bit odd, how configparser makes things lowercase(?)	
		if modelClass == 'codon':
			self.model.params = {'alpha': 1., 'beta': 1., 'mu':{'AC': 1., 'AG': 1., 'AT': 1., 'CG': 1., 'CT': 1., 'GT': 1.}}
			self.parseCodonParams()
			print self.model.params
		elif modelClass == 'mutsel':
			self.model.params = {'alpha': 1., 'beta': 1., }
		elif modelClass == 'nucleotide':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
		elif modelClass == 'amino':
			self.model.params = {'myparams': 'fill in later once this class is constructed.'}
		else:
			raise AssertionError("You must provide a model class (codon, mutsel, nucleotide, or amino).")
	



	def parseMutselParams(self):
		''' Parse the provided mutation-selection model parameters. Use these to overwrite any defaults.
			Allowed input parameters: mu, kappa.		
		'''
		provided = self.options['PARAMETERS']
		# Grab any mu and/or kappa to overwrite
		self.parseCommonMK(provided)
				
	def parseCodonParams(self):
		''' Parse the provided codon model parameters. Use these to overwrite any defaults.
			Allowed input parameters: mu, kappa, omega OR alpha/beta.		
		'''
		provided = self.options['PARAMETERS']
		
		# Check that reasonable combinations are provided. Expand on this part later.
		if "alpha" in provided and "omega" in provided:
			print "DANGER, WILL ROBINSON!: You specified value(s) for alpha (dS) and/or beta(dN) in addition to omega (dN/dS). Therefore, I will only pay attention to your provided omega value."
		
		# Grab any mu and/or kappa to overwrite
		self.parseCommonMK(provided)
		
		if "omega" in provided:
			try:
				omega = float(provided['omega'])
			except TypeError:
				"If you wish to provide an omega (equivalent to dN/dS) value, it must be numeric."
			self.model.params['beta'] = omega
			self.model.params['alpha'] = 1.0	
		
		else:
			if "alpha" in provided:
				try:
					alpha = float(provided['alpha'])
				except TypeError:
					"If you wish to provide an alpha (equivalent to dS) value, it must be numeric."
				self.model.params['alpha'] = alpha
			
			if "beta" in provided:
				try:
					beta = float(provided['beta'])
				except TypeError:
					"If you wish to provide a beta (equivalent to dN) value, it must be numeric."
				self.model.params['beta'] = beta
			
			
			
		
		def parseCommonMK(self, provided):
			''' codon, mutsel, and nucleotide all use mu and kappa. Use this function for all three.'''
		
			if "mu" in provided:
				provided_mu = eval(provided['mu'])
				for k in provided_mu:
					if provided_mu[k] != 1.:
						self.model.params['mu'][k] = provided_mu[k] 
			
			## NOT YET CLEAR TO ME IF THIS SHOULD BE AN IF OR ELSE!!!!!!!! ARE THESE MUTUALLY EXCLUSIVE? NEED TO THINK ON IT.
			if "kappa" in provided:
				try:
					kappa = float(provided['kappa'])
				except TypeError:
					"If you wish to provide a kappa (TI/TV) value, it must be numeric."
				self.model.params['AG'] = self.model.params['AG'] * kappa
				self.model.params['GA'] = self.model.params['GA'] * kappa
				self.model.params['CT'] = self.model.params['CT'] * kappa
				self.model.params['TC'] = self.model.params['TC'] * kappa
				
				
				
				
				
				
				
				
			

conf = parseConfig()
conf.parseSections()
conf.configModel()


#mu = eval(options['PARAMETERS']['muparams']) # this works.