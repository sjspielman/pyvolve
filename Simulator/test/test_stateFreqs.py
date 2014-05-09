import unittest
import sys
import numpy as np

SRC_CODE = "../src/"  # Path to source code.
sys.path.append(SRC_CODE)
from misc import Genetics
from stateFreqs import *

'''
class stateFreqs_basic_functions_bycodon(unittest.TestCase):
	''''''Unit testing for the base class functions in stateFreqs when by=codon '''''' 
	
	def setUp(self):
		self.codonFreqs = 
'''

class stateFreqs_EqualFreqs_Tests(unittest.TestCase):
	''' Set of unittests for the EqualFreqs subclass of StateFreqs.'''
	#correct_nuc = [0.24043715846994534, 0.26229508196721313, 0.25136612021857924, 0.2459016393442623]
		
	def test_EqualFreqs_generate_bycodon(self):
		print "Testing EqualFreqs generate function, for by=codon."
		correct = np.array(np.repeat(1./61., 61))
		self.eqFreqs = EqualFreqs( by = 'codon' )
		freqs = self.eqFreqs.generate()
		np.testing.assert_array_almost_equal(correct, freqs, decimal = 10, err_msg = "Codon EqualFreqs not generated properly for by=codon.")
	
	def test_EqualFreqs_generate_byamino(self):
		print "Testing EqualFreqs generate function, for by=amino."
		correct = np.array(np.repeat(1./20., 20))
		self.eqFreqs = EqualFreqs( by = 'amino' )
		freqs = self.eqFreqs.generate()	
		np.testing.assert_array_almost_equal(correct, freqs, decimal = 10, err_msg = "Codon EqualFreqs not generated properly for by=amino.")

	
	
	
	
	
	
	
if __name__ == '__main__':
	run_tests = unittest.TextTestRunner()
	test_suite = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_EqualFreqs_Tests)
	#test_suite = suite()
	run_tests.run(test_suite)
	