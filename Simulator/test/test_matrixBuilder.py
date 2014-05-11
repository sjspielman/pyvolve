####### ALL TESTS HERE PROVIDE INPUTS PROPERLY ###########
import unittest
import sys
import numpy as np

SRC_CODE = "../src/"  # Path to source code.
sys.path.append(SRC_CODE)
from matrixBuilder import *
from misc import Model


class matrixBuilder_baseClass_Tests(unittest.TestCase):
	''' 
		Set of unittests for simple base class functions in matrixBuilder.
		Functions tested here include isTI, isSyn, getCodonFreq. 
		All other functions require full model specification, so they are tested elsewhere.
		Note: stateFreqs specification is required - we will simply give it equal frequencies.
	'''
	
	def setUp(self):
		basicModel = Model()
		freqs = np.array(np.repeat(1./61., 61))
		basicModel.params = {'stateFreqs':freqs}
		self.baseObject = MatrixBuilder(basicModel)
		self.dec = 8
		
	def test_matrixBuilder_baseClass_isTI(self):	
		''' Test that transitions can be properly identified. '''
		
		self.assertTrue( self.baseObject.isTI('A', 'G'), msg = "matrixBuilder.isTI() does not think A -> G is a transition.")
		self.assertTrue( self.baseObject.isTI('G', 'A'), msg = "matrixBuilder.self.baseObject.isTI() does not think G -> A is a transition.")
		self.assertTrue( self.baseObject.isTI('C', 'T'), msg = "matrixBuilder.self.baseObject.isTI() does not think C -> T is a transition.")
		self.assertTrue( self.baseObject.isTI('T', 'C'), msg = "matrixBuilder.self.baseObject.isTI() does not think C -> T is a transition.")
		self.assertFalse( self.baseObject.isTI('A', 'C'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks A -> C is a transition.")
		self.assertFalse( self.baseObject.isTI('C', 'A'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks C -> A is a transition.")
		self.assertFalse( self.baseObject.isTI('A', 'T'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks A -> T is a transition.")
		self.assertFalse( self.baseObject.isTI('T', 'A'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks T -> A is a transition.")
		self.assertFalse( self.baseObject.isTI('G', 'C'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks G -> C is a transition.")
		self.assertFalse( self.baseObject.isTI('C', 'G'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks C -> G is a transition.")
		self.assertFalse( self.baseObject.isTI('G', 'T'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks G -> T is a transition.")
		self.assertFalse( self.baseObject.isTI('T', 'G'), msg = "matrixBuilder.self.baseObject.isTI() mistakenly thinks T -> G is a transition.")

if __name__ == '__main__':
	run_tests = unittest.TextTestRunner()
	

	print "Testing the simple functions in the base class matrixBuilder"
	test_suite_baseMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_baseClass_Tests)
	run_tests.run(test_suite_baseMatrix)