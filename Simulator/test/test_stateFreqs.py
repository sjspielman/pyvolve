import unittest
import sys
import numpy as np

SRC_CODE = "../src/"  # Path to source code.
sys.path.append(SRC_CODE)
from misc import Genetics
from stateFreqs import *




class stateFreqs_RandFreqs_Tests(unittest.TestCase):
	''' Set of "unittests" for the EqualFreqs subclass of StateFreqs. Note that since random, cannot test exact values.'''
	
	def setUp(self):
		self.dec = 8

	############### type is not provided - should raise assertion #################
	def test_RandFreqs_calcFreqs_bycodon_notype(self):
		self.rFreqs = RandFreqs(by = 'codon')
		self.assertRaises(AssertionError, self.rFreqs.calcFreqs)
	def test_RandFreqs_calcFreqs_byamino_notype(self):
		self.rFreqs = RandFreqs(by = 'amino')
		self.assertRaises(AssertionError, self.rFreqs.calcFreqs)
	def test_RandFreqs_calcFreqs_bynuc_notype(self):
		self.rFreqs = RandFreqs(by = 'nuc')
		self.assertRaises(AssertionError, self.rFreqs.calcFreqs)
	
	############################# by=codon tests ##################################
	def test_RandFreqs_calcFreqs_bycodon_typecodon(self):
		correct_len = 61
		self.rFreqs = RandFreqs( by = 'codon', type = 'codon' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=codon.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")
	
	def test_RandFreqs_calcFreqs_bycodon_typeamino(self):
		correct_len = 20
		self.rFreqs = RandFreqs( by = 'codon', type = 'amino' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=amino.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")

	def test_RandFreqs_calcFreqs_bycodon_typenuc(self):
		correct_len = 4
		self.rFreqs = RandFreqs( by = 'codon', type = 'nuc' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=amino.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")
	
	############################# by=amino tests ##################################
	def test_RandFreqs_calcFreqs_byamino_typecodon(self):
		correct_len = 61
		self.rFreqs = RandFreqs( by = 'amino', type = 'codon' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=codon.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=codon.")
	
	def test_RandFreqs_calcFreqs_byamino_typeamino(self):
		correct_len = 20
		self.rFreqs = RandFreqs( by = 'amino', type = 'amino' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=amino.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=amino.")

	def test_RandFreqs_calcFreqs_byamino_typenuc(self):
		correct_len = 4
		self.rFreqs = RandFreqs( by = 'amino', type = 'nuc' )
		freqs = self.rFreqs.calcFreqs()
		print "here",freqs
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=nuc.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=nuc.")

	############################# by=nuc tests ##################################
	def test_RandFreqs_calcFreqs_bynuc_typecodon(self):
		correct_len = 4
		self.rFreqs = RandFreqs( by = 'nuc', type = 'codon' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=codon.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=nuc, type=codon.")

	def test_RandFreqs_calcFreqs_bynuc_typeamino(self):
		correct_len = 4
		self.rFreqs = RandFreqs( by = 'nuc', type = 'amino' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=amino.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=nuc, type=amino.")
	
	def test_RandFreqs_calcFreqs_bynuc_typenuc(self):
		correct_len = 4
		self.rFreqs = RandFreqs( by = 'nuc', type = 'nuc' )
		freqs = self.rFreqs.calcFreqs()
		np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=nuc.")
		self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=nuc, type=nuc.")









class stateFreqs_EqualFreqs_Tests(unittest.TestCase):
	''' Set of "unittests" for the EqualFreqs subclass of StateFreqs.'''
	
	def setUp(self):
		self.dec = 8 # For accuracy
	
	############### type is not provided - should raise assertion #################
	def test_EqualFreqs_calcFreqs_bycodon_notype(self):
		self.eqFreqs = EqualFreqs(by = 'codon')
		self.assertRaises(AssertionError, self.eqFreqs.calcFreqs)
	def test_EqualFreqs_calcFreqs_byamino_notype(self):
		self.eqFreqs = EqualFreqs(by = 'amino')
		self.assertRaises(AssertionError, self.eqFreqs.calcFreqs)
	def test_EqualFreqs_calcFreqs_bynuc_notype(self):
		self.eqFreqs = EqualFreqs(by = 'nuc')
		self.assertRaises(AssertionError, self.eqFreqs.calcFreqs)
	
	
	############################# by=codon tests ##################################
	
	def test_EqualFreqs_calcFreqs_bycodon_typecodon(self):
		correct = np.array(np.repeat(1./61., 61))
		self.eqFreqs = EqualFreqs( by = 'codon', type = 'codon' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=codon.")
	
	def test_EqualFreqs_calcFreqs_bycodon_typeamino(self):
		correct = np.array([0.0655737704918, 0.0327868852459, 0.0327868852459, 0.0327868852459, 0.0327868852459, 0.0655737704918, 0.0327868852459, 0.0491803278689, 0.0327868852459, 0.0983606557377, 0.016393442623, 0.0327868852459, 0.0655737704918, 0.0327868852459, 0.0983606557377, 0.0983606557377, 0.0655737704918, 0.0655737704918, 0.016393442623, 0.0327868852459])
		self.eqFreqs = EqualFreqs( by = 'codon', type = 'amino' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=amino.")
	
	def test_EqualFreqs_calcFreqs_bycodon_typenuc(self):
		correct = np.array([0.24043715846994534, 0.26229508196721313, 0.25136612021857924, 0.2459016393442623])
		self.eqFreqs = EqualFreqs( by = 'codon', type = 'nuc' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=nuc.")

	############################# by=amino tests ##################################
	
	def test_EqualFreqs_calcFreqs_byamino_typecodon(self):
		correct = np.array([0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.01666667, 0.01666667, 0.05, 0.01666667, 0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.025, 0.025, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.025, 0.05, 0.025, 0.00833333, 0.025, 0.00833333, 0.025])
		self.eqFreqs = EqualFreqs( by = 'amino', type = 'codon' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=codon.")
	
	def test_EqualFreqs_calcFreqs_byamino_typeamino(self):
		correct = np.array(np.repeat(1./20., 20))
		self.eqFreqs = EqualFreqs( by = 'amino', type = 'amino' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=codon.")
	
	def test_EqualFreqs_calcFreqs_byamino_typenuc(self):
		correct = np.array(np.repeat(0.25, 4))
		self.eqFreqs = EqualFreqs( by = 'amino', type = 'nuc' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=nuc.")
	
	
	############################# by=nuc tests ##################################
	def test_EqualFreqs_calcFreqs_bynuc_typecodon(self):
		correct = np.array(np.repeat(0.25, 4))
		self.eqFreqs = EqualFreqs( by = 'nuc', type = 'codon' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=nuc.")

	def test_EqualFreqs_calcFreqs_bynuc_typeamino(self):
		correct = np.array(np.repeat(0.25, 4))
		self.eqFreqs = EqualFreqs( by = 'nuc', type = 'amino' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=nuc.")

	def test_EqualFreqs_calcFreqs_bynuc_typenuc(self):
		correct = np.array(np.repeat(0.25, 4))
		self.eqFreqs = EqualFreqs( by = 'nuc', type = 'nuc' )
		freqs = self.eqFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=nuc.")


























class stateFreqs_UserFreqs_Tests(unittest.TestCase):
	''' Set of "unittests" for the UserFreqs subclass of StateFreqs.'''
	
	def setUp(self):
		self.dec = 8 # For accuracy


	########### incorrect dictionaries provided. should raise assertions. ##############
	def test_UserFreqs_calcFreqs_bycodon_badkey_length(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'A':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_byamino_badkey_length(self):
		self.uFreqs = UserFreqs(by = 'amino', freqs = {'AAA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
		
	def test_UserFreqs_calcFreqs_bycodon_keytoolong(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAAAAA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
		
	def test_UserFreqs_calcFreqs_bycodon_keytooshort(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_byamino_keytoolong(self):
		self.uFreqs = UserFreqs(by = 'amino', freqs = {'AA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
		
	def test_UserFreqs_calcFreqs_bycodon_badkey_letters(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'ABC':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_bycodon_badkey_numbers(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {123:1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toosmall(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toobigsingle(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':1.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toobigmultiple(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0.7, 'AAT':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
		
	def test_UserFreqs_calcFreqs_bycodon_badvalues_zero(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_negativedecimal(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':-0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_bycodon_badvalues_negativeaboveabs1(self):
		self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':-2})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)


	################## do not provide by or type. assertions shoud be raised. ####################

	def test_UserFreqs_calcFreqs_noby_notype_badtriplet(self):
		self.uFreqs = UserFreqs(freqs = {'WWW':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_noby_notype_badsingleletter(self):
		self.uFreqs = UserFreqs(freqs = {'X':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)

	def test_UserFreqs_calcFreqs_noby_notype_ambignucaa(self):
		self.uFreqs = UserFreqs(freqs = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)


	############ provide type but not by. assertions should be raised. ##################
	def test_UserFreqs_calcFreqs_noby_aminotype_wrongdict_singlecodon(self):
		self.uFreqs = UserFreqs(type = 'amino', freqs = {'GCA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)

	def test_UserFreqs_calcFreqs_noby_aminotype_wrongdict_multiplecodons(self):
		self.uFreqs = UserFreqs(type = 'amino', freqs = {'GCA':0.5, 'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)

	def test_UserFreqs_calcFreqs_noby_nuctype_wrongdict_singlecodon(self):
		self.uFreqs = UserFreqs(type = 'nuc', freqs = {'GCA':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_noby_nuctype_wrongdict_multiplecodons(self):
		self.uFreqs = UserFreqs(type = 'nuc', freqs = {'GCA':0.5, 'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
		
	def test_UserFreqs_calcFreqs_noby_codontype_wrongdict_singleaa(self):
		self.uFreqs = UserFreqs(type = 'codon', freqs = {'W':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_noby_codontype_wrongdict_multipleaa(self):
		self.uFreqs = UserFreqs(type = 'codon', freqs = {'W':0.5, 'D':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
		
	def test_UserFreqs_calcFreqs_noby_codontype_wrongdict_ambignucaa(self):
		self.uFreqs = UserFreqs(type = 'codon', freqs = {'A':0.5, 'T':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
			
	def test_UserFreqs_calcFreqs_noby_notype_baddict_singleambig(self):
		self.uFreqs = UserFreqs(freqs = {'A':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	def test_UserFreqs_calcFreqs_noby_notype_baddict_multipleambig(self):
		self.uFreqs = UserFreqs(freqs = {'A':0.5, 'C':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
    
    
	
	############# do not provide by, but everything should work out ok ###################
	def test_UserFreqs_calcFreqs_noby_notype_gooddict_singlecodon(self):
		correct = np.zeros(61)
		correct[0] = 1.0
		self.uFreqs = UserFreqs(freqs = {'AAA':1.0})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_noby_notype_gooddict_mulcodons(self):
		correct = np.zeros(61)
		correct[0] = 0.5
		correct[60] =0.5
		self.uFreqs = UserFreqs(freqs = {'AAA':0.5, 'TTT':0.5})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_noby_notype_gooddict_singleaa(self):
		correct = np.zeros(20)
		correct[2] = 1.0
		self.uFreqs = UserFreqs(freqs = {'D':1.0})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
		
	def test_UserFreqs_calcFreqs_noby_notype_gooddict_multipleaa(self):
		correct = np.zeros(20)
		correct[2] = 0.5
		correct[3] = 0.5
		self.uFreqs = UserFreqs(freqs = {'D':0.5, 'E':0.5})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_noby_aminotype_gooddict_ambignucaa(self):
		correct = np.zeros(20)
		correct[0] = 0.5
		correct[1] = 0.5
		self.uFreqs = UserFreqs(type = 'amino', freqs = {'A':0.5, 'C':0.5})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_noby_nuctype_gooddict_ambignucaa(self):
		correct = np.zeros(4)
		correct[0] = 0.5
		correct[1] = 0.5
		self.uFreqs = UserFreqs(type = 'nuc', freqs = {'A':0.5, 'C':0.5})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")

    ############## provide everything correctly, and calculate everything correctly ############
	def test_UserFreqs_calcFreqs_byamino_typeamino_gooddict_singleaa(self):
		correct = np.zeros(20)
		correct[1] = 1.0
		self.uFreqs = UserFreqs(by='amino', type='amino', freqs = {'C':1.0})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_byamino_typeamino_gooddict_multipleaa(self):
		correct = np.zeros(20)
		correct[1] = 0.25
		correct[2] = 0.25
		correct[3] = 0.25
		correct[4] = 0.25
		self.uFreqs = UserFreqs(by='amino', type='amino', freqs = {'C':0.25, 'D':0.25, 'E':0.25, 'F':0.25})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_bynuc_typenuc_gooddict_singlenuc(self):
		correct = np.zeros(4)
		correct[1] = 1.0
		self.uFreqs = UserFreqs(by='nuc', type='nuc', freqs = {'C':1.0})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_bynuc_typenuc_gooddict_multiplenuc(self):
		correct = np.zeros(4)
		correct[0] = 0.25
		correct[1] = 0.25
		correct[2] = 0.25
		correct[3] = 0.25
		self.uFreqs = UserFreqs(by='nuc', type='nuc', freqs = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	
	def test_UserFreqs_calcFreqs_bycodon_typecodon_gooddict_singlecodon(self):
		correct = np.zeros(61)
		correct[1] = 1.0
		self.uFreqs = UserFreqs(by='codon', type='codon', freqs = {'AAC':1.0})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
	def test_UserFreqs_calcFreqs_bycodon_typecodon_gooddict_multiplecodons(self):
		correct = np.zeros(61)
		correct[1] = 0.5
		correct[2] = 0.5
		self.uFreqs = UserFreqs(by='codon', type='codon', freqs = {'AAC':0.5, 'AAG':0.5})
		freqs = self.uFreqs.calcFreqs()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for no by provided with correct codon data.")
	
    
    
  

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
		self.rFreqs = ReadFreqs(by = 'amino', columns[0,1,2], file = 'freqFiles/testFreq_amino_notaln.fasta')
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

	### above here are all ok.
	# I still need to confirm that the columns selected are actually within range in the alignment (are there that many columns...)
	
   
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
	run_tests = unittest.TextTestRunner()
	#test_suite_Equal = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_EqualFreqs_Tests)
	#run_tests.run(test_suite_Equal)
	#test_suite_Rand = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_RandFreqs_Tests)
	#run_tests.run(test_suite_Rand)
	#test_suite_User = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_UserFreqs_Tests)
	#run_tests.run(test_suite_User)
	test_suite_Read = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_ReadFreqs_Tests)
	run_tests.run(test_suite_Read)
	
	
	
	
	
	
	
	
	
	
	
	