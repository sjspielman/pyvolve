## Parser for temp_config.txt. Note that this is almost definitely not a final parser to be released, but for current use to become comfortable with ConfigParser module, etc.

import os
import sys
from string import whitespace
from misc import Model

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
		
	
	def configModel(self):
		''' Configure the default evolutionary model specifications.
			options['modelClass'] should be either mutsel, codon, amino, nucleotide.
			Depending on which model type, construct default model specifications. 
			These will automatically be overwritten by anything different self.options in options['PARAMETERS'].
		'''

		if self.options['modelClass'] == 'codon':
			self.model.params = {'alpha': 1., 'beta': 1., 'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.} }
			self.parseCodonParams()
			print self.model.params
		
		elif self.options['modelClass'] == 'mutsel':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}

		elif self.options['modelClass'] == 'nucleotide':
			self.model.params = {'mu':{'AC': 1., 'CA':1., 'AG': 1., 'GA':1., 'AT': 1., 'TA':1., 'CG': 1., 'GC':1., 'CT': 1., 'TC':1., 'GT': 1., 'TG':1.}}
		
		elif self.options['modelClass'] == 'amino':
			self.model.params = {'myparams': 'fill in later once this class is constructed.'}
		
		else:
			raise AssertionError("You must provide a model class (codon, mutsel, nucleotide, or amino).")


	### NOT YET CLEAR TO ME HOW TO DEAL WITH NUC AND MUTSEL SINCE HAVE THE SAME PARAMETERS OF MU AND KAPPA.
	def parseMutselParams(self):
		''' Parse the self.options mutation-selection model parameters. Use these to overwrite any defaults.
			Allowed input parameters: mu, kappa.		
		'''
		self.parseCommonMK()
				
	def parseCodonParams(self):
		''' Parse the self.options codon model parameters. Use these to overwrite any defaults.
			Allowed input parameters: mu, kappa, omega OR alpha/beta.		
		'''
		
		# Check that reasonable combinations are self.options. Expand on this part later.
		if "alpha" in self.options and "omega" in self.options:
			print "\n\nDANGER, WILL ROBINSON!: You specified value(s) for alpha (dS) and/or beta(dN) in addition to omega (dN/dS). Therefore, I will only pay attention to your self.options omega value.\n\n"
		
		# Grab any mu and/or kappa to overwrite
		self.parseCommonMK()
		
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
			
			
			
		
	def parseCommonMK(self):
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
				
				
				
				
				
				
				
				
			

conf = parseConfig()
conf.configModel()


#mu = eval(options['PARAMETERS']['muparams']) # this works.