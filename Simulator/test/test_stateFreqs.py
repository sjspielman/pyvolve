import unittest
import sys
import numpy as np

SRC_CODE = "../src/"  # Path to source code.
sys.path.append(SRC_CODE)
from misc import Genetics
from stateFreqs import *




class stateFreqs_RandFreqs_Tests(unittest.TestCase):
	''' Set of unittests for the EqualFreqs subclass of StateFreqs. Note that since random, cannot test exact values.'''
	
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
	''' Set of unittests for the EqualFreqs subclass of StateFreqs.'''
	
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
	''' Set of unittests for the UserFreqs subclass of StateFreqs.'''
	
	def setUp(self):
		self.dec = 8 # For accuracy
	
	########### incorrect dictionaries provided. should raise assertions. ##############
	def test_UserFreqs_calcFreqs_bycodon_badkey_length(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {'A':0.5, 'T':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	def test_UserFreqs_calcFreqs_bycodon_badkey_letters(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {'ABC':1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	def test_UserFreqs_calcFreqs_bycodon_badkey_numbers(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {123:1.0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toobig(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {'AAA':1.1})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	def test_UserFreqs_calcFreqs_bycodon_badvalues_zero(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {'AAA':0})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	def test_UserFreqs_calcFreqs_bycodon_badvalues_negative(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {'AAA':-1})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	def test_UserFreqs_calcFreqs_bycodon_badvalues_toosmall(self):
		self.uFreqs = UserFreqs(by = 'codon', type = 'codon', freqs = {'AAA':0.5})
		self.assertRaises(AssertionError, self.uFreqs.calcFreqs)
	
	########################### guess the by= from the keys ############################
	#def test_UserFreqs_calcFreqs_bycodon_noby_badtriplet(self):
	#	self.uFreqs = UserFreqs(type = 'codon', freqs = {'WWW':1.0})
	#	self.assertRaises(AssertionError, self.uFreqs.calcFreqs)





    
if __name__ == '__main__':
	run_tests = unittest.TextTestRunner()
	#test_suite_Equal = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_EqualFreqs_Tests)
	#run_tests.run(test_suite_Equal)
	#test_suite_Rand = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_RandFreqs_Tests)
	#run_tests.run(test_suite_Rand)
	test_suite_User = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_UserFreqs_Tests)
	run_tests.run(test_suite_User)
	
	
	
	
	
	
	
	
	
	
	
	
	