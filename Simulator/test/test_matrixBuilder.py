####### ALL TESTS HERE PROVIDE INPUTS PROPERLY ###########
import unittest
import sys
import numpy as np
from Bio import Seq
from Bio.Alphabet import generic_dna

SRC_CODE = "../src/"  # Path to source code.
sys.path.append(SRC_CODE)
from matrixBuilder import *
from misc import Model


class matrixBuilder_baseClass_Tests(unittest.TestCase):
	''' 
		Set of unittests for simple base class functions in matrixBuilder.
		Functions tested here include isTI, isSyn, getCodonFreq. 
		All other functions require full model specification, so they are tested elsewhere.
		Note: stateFreqs specification is required. To test this, all frequencies must be unique!.
	'''
	
	def setUp(self):
		# Do not rely on misc for codons in case something happens to it!
		self.codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
		basicModel = Model()
		self.codonfreqs = [0.01617666, 0.00291771, 0.02664918, 0.02999061, 0.00717921, 0.00700012, 0.01435559, 0.0231568, 0.02403056, 0.00737008, 0.03185765, 0.0193576, 0.03277142, 0.02141258, 0.0127537, 0.00298803, 0.0256333, 0.02312437, 0.01861465, 0.01586447, 0.00373147, 0.02662654, 0.00082524, 0.00048916, 0.01191673, 0.00512658, 0.00050502, 0.01688169, 0.01843001, 0.00215437, 0.02659356, 0.02377742, 0.01169375, 0.00097256, 0.02937344, 0.00268204, 0.01414414, 0.02781933, 0.00070877, 0.02370841, 0.02984617, 0.01828081, 0.01002825, 0.00870788, 0.00728006, 0.02179328, 0.00379049, 0.01978996, 0.00443774, 0.01201798, 0.02030269, 0.01238501, 0.01279963, 0.02094385, 0.02810987, 0.00918507, 0.02880549, 0.0029311, 0.0237658, 0.03194712, 0.06148723]
		basicModel.params = {'stateFreqs':self.codonfreqs}
		self.baseObject = MatrixBuilder(basicModel)
		self.zero = 1e-8
		
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

	def test_matrixBuilder_baseClass_isSyn(self):	
		''' Test that synonymous vs nonsynymous changes can be properly identified. 
			Assumes that biopython is not broken. This is (theoretically...) a very safe assumption.
		'''
		for source in self.codons:
			for target in self.codons:
				source_aa = str( Seq.Seq(source, generic_dna).translate() )
				target_aa = str( Seq.Seq(target, generic_dna).translate() )
				if source_aa == target_aa:
					self.assertTrue( self.baseObject.isSyn(source, target), msg = ("matrixBuilder.isSyn() does not think", source, " -> ", target, " is synonymous.") )
				else:
					self.assertFalse( self.baseObject.isSyn(source, target), msg = ("matrixBuilder.isSyn() mistakenly thinks", source, " -> ", target, " is synonymous.") )
	
	def test_matrixBuilder_baseClass_getCodonFreq(self):
		''' Test that, given a codon, base frequency is properly identified. '''
		for i in range(61):
			codon = self.codons[i]
			correct_freq = self.codonfreqs[i]
			self.assertTrue( (abs(self.baseObject.getCodonFreq(codon) - correct_freq) < self.zero), msg = "matrixBuilder.getCodonFreq doesn't work properly.")
			
if __name__ == '__main__':
	run_tests = unittest.TextTestRunner()
	

	print "Testing the simple functions in the base class matrixBuilder"
	test_suite_baseMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_baseClass_Tests)
	run_tests.run(test_suite_baseMatrix)