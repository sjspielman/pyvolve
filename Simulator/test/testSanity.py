## SJS 5/11/14.
######## sanitycheck tests. these were written for stateFreqs, but then I realized that all parsing stuff should be its own class/module, clearly, so moving functions here for now until can use them in a better way.
######## these almost definitely will NOT WORK PROPERLY, but will guide future coding when I need to reimplement it properly in its own class.


class stateFreqs_EqualFreqs_Tests(unittest.TestCase):
	''' Set of "unittests" for the EqualFreqs subclass of StateFreqs.'''
	
	def setUp(self):
		self.dec = 8 # For accuracy
	

	############### type is not provided - should raise assertion #################
	def test_EqualFreqs_calcFreqs_bycodon_notype(self):
		self.eqFreqs = EqualFreqs(by = 'codon')
		self.assertRaises(AssertionError, self.eqFreqs.sanityCheck)
	def test_EqualFreqs_calcFreqs_byamino_notype(self):
		self.eqFreqs = EqualFreqs(by = 'amino')
		self.assertRaises(AssertionError, self.eqFreqs.sanityCheck)
	def test_EqualFreqs_calcFreqs_bynuc_notype(self):
		self.eqFreqs = EqualFreqs(by = 'nuc')
		self.assertRaises(AssertionError, self.eqFreqs.sanityCheck)
	
	

class stateFreqs_UserFreqs_Tests(unittest.TestCase):
	''' Set of "unittests" for the UserFreqs subclass of StateFreqs.'''
	
	def setUp(self):
		self.dec = 8 # For accuracy


	########### incorrect dictionaries provided. should raise assertions. ##############
	def test_UserFreqs_calcFreqs_bycodon_badkey_length(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'A':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_byamino_badkey_length(self):
		self.uFreqs = UserFreqs(by = 'amino', freqs = {'AAA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
		
	def test_UserFreqs_calcFreqs_bycodon_keytoolong(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAAAAA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
		
	def test_UserFreqs_calcFreqs_bycodon_keytooshort(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_byamino_keytoolong(self):
		self.uFreqs = UserFreqs(by = 'amino', freqs = {'AA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
		
	def test_UserFreqs_calcFreqs_bycodon_badkey_letters(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'ABC':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_bycodon_badkey_numbers(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {123:1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toosmall(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toobigsingle(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':1.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toobigmultiple(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0.7, 'AAT':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
		
	def test_UserFreqs_calcFreqs_bycodon_badvalues_zero(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_negativedecimal(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':-0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_negativeaboveabs1(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':-2})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)


	################## do not provide by or type. assertions shoud be raised. ####################

	def test_UserFreqs_calcFreqs_noby_notype_badtriplet(self):
		self.uFreqs = UserFreqs(freqs = {'WWW':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_noby_notype_badsingleletter(self):
		self.uFreqs = UserFreqs(freqs = {'X':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)

	def test_UserFreqs_calcFreqs_noby_notype_ambignucaa(self):
		self.uFreqs = UserFreqs(freqs = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)


	############ provide type but not by. assertions should be raised. ##################
	def test_UserFreqs_calcFreqs_noby_aminotype_wrongdict_singlecodon(self):
		self.uFreqs = UserFreqs(type = 'amino', freqs = {'GCA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)

	def test_UserFreqs_calcFreqs_noby_aminotype_wrongdict_multiplecodons(self):
		self.uFreqs = UserFreqs(type = 'amino', freqs = {'GCA':0.5, 'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)

	def test_UserFreqs_calcFreqs_noby_nuctype_wrongdict_singlecodon(self):
		self.uFreqs = UserFreqs(type = 'nuc', freqs = {'GCA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_noby_nuctype_wrongdict_multiplecodons(self):
		self.uFreqs = UserFreqs(type = 'nuc', freqs = {'GCA':0.5, 'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
		
	def test_UserFreqs_calcFreqs_noby_codontype_wrongdict_singleaa(self):
		self.uFreqs = UserFreqs(type = 'codon', freqs = {'W':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_noby_codontype_wrongdict_multipleaa(self):
		self.uFreqs = UserFreqs(type = 'codon', freqs = {'W':0.5, 'D':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
		
	def test_UserFreqs_calcFreqs_noby_codontype_wrongdict_ambignucaa(self):
		self.uFreqs = UserFreqs(type = 'codon', freqs = {'A':0.5, 'T':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
			
	def test_UserFreqs_calcFreqs_noby_notype_baddict_singleambig(self):
		self.uFreqs = UserFreqs(freqs = {'A':1.0})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
	
	def test_UserFreqs_calcFreqs_noby_notype_baddict_multipleambig(self):
		self.uFreqs = UserFreqs(freqs = {'A':0.5, 'C':0.5})
		self.assertRaises(AssertionError, self.uFreqs.sanityCheck)
    


class stateFreqs_ReadFreqs_Tests(unittest.TestCase):
	''' Set of "unittests" for the ReadFreqs subclass of StateFreqs.'''
	
	def setUp(self):
		self.dec = 8 # For accuracy


	########### file/column specifications which should raise assertions ##############

	def test_ReadFreqs_sanityCheck_bycodon_nocol_badfile(self):
		self.rFreqs = ReadFreqs(by = 'codon', file = 'NonexistentFilename.txt')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_byamino_nocol_badfile(self):
		self.rFreqs = ReadFreqs(by = 'amino', file = 'NonexistentFilename.txt')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_bycodon_nocol_aminofile(self):
		self.rFreqs = ReadFreqs(by = 'codon', file = 'freqFiles/testFreq_amino_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_bynuc_nocol_aminofile(self):
		self.rFreqs = ReadFreqs(by = 'nuc', file = 'freqFiles/testFreq_amino_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)

	def test_ReadFreqs_sanityCheck_byamino_nocol_codonfile(self):
		self.rFreqs = ReadFreqs(by = 'amino', file = 'freqFiles/testFreq_codon_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_byamino_nocol_bonkersfile(self):
		self.rFreqs = ReadFreqs(by = 'amino', file = 'freqFiles/testFreq_bonkers_letters_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_bynuc_nocol_bonkersfile(self):
		self.rFreqs = ReadFreqs(by = 'nuc', file = 'freqFiles/testFreq_bonkers_letters_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
		
	def test_ReadFreqs_sanityCheck_byamino_goodcol_notaln(self):
		self.rFreqs = ReadFreqs(by = 'amino', columns=[0,1,2], file = 'freqFiles/testFreq_amino_notaln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_bycodon_goodcol_notaln(self):
		self.rFreqs = ReadFreqs(by = 'codon', columns=[0,1,2], file = 'freqFiles/testFreq_codon_notaln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_bycodon_badcol_colstring(self):
		self.rFreqs = ReadFreqs(by = 'amino', columns = "thisisnotacolumnlist", file = 'freqFiles/testFreq_codon_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)
	
	def test_ReadFreqs_sanityCheck_bycodon_badcol_colfloat(self):
		self.rFreqs = ReadFreqs(by = 'amino', columns = 8.345193, file = 'freqFiles/testFreq_codon_aln.fasta')
		self.assertRaises(AssertionError, self.rFreqs.sanityCheck)