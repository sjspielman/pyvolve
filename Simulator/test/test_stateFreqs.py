####### ALL TESTS HERE PROVIDE INPUTS PROPERLY ###########
import unittest
import sys
import numpy as np

SRC_CODE = "../src/"  # Path to source code.
sys.path.append(SRC_CODE)
from stateFreqs import *


class stateFreqs_RandFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the EqualFreqs subclass of StateFreqs. Note that since random, cannot test exact values.'''
    
    def setUp(self):
        self.dec = 8
    
    ############################# by=codon tests ##################################
    def test_RandFreqs_calcFreqs_bycodon_typecodon_restriction(self):
        correct_len = 61
        self.rFreqs = RandFreqs( by = 'codon', restrict = ['AGG', 'CCT', 'GCG', 'TGC'] )
        freqs = self.rFreqs.calcFreqs(type = 'codon')
        sum = freqs[10] + freqs[23] + freqs[38] + freqs[54]
        np.testing.assert_almost_equal(sum, 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=codon, with restriction.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")
           
    
    def test_RandFreqs_calcFreqs_bycodon_typecodon(self):
        correct_len = 61
        self.rFreqs = RandFreqs( by = 'codon' )
        freqs = self.rFreqs.calcFreqs(type = 'codon')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=codon.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")
    
    def test_RandFreqs_calcFreqs_bycodon_typeamino(self):
        correct_len = 20
        self.rFreqs = RandFreqs( by = 'codon' )
        freqs = self.rFreqs.calcFreqs(type = 'amino')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=amino.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")

    def test_RandFreqs_calcFreqs_bycodon_typenuc(self):
        correct_len = 4
        self.rFreqs = RandFreqs( by = 'codon' )
        freqs = self.rFreqs.calcFreqs(type = 'nuc')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=nuc.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")
    

    ############################# by=amino tests ##################################
    def test_RandFreqs_calcFreqs_byamino_typecodon(self):
        correct_len = 61
        self.rFreqs = RandFreqs( by = 'amino' )
        freqs = self.rFreqs.calcFreqs(type = 'codon')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=codon.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=codon.")
    
    def test_RandFreqs_calcFreqs_byamino_typeamino(self):
        correct_len = 20
        self.rFreqs = RandFreqs( by = 'amino' )
        freqs = self.rFreqs.calcFreqs(type = 'amino')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=amino.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=amino.")

    def test_RandFreqs_calcFreqs_byamino_typenuc(self):
        correct_len = 4
        self.rFreqs = RandFreqs( by = 'amino' )
        freqs = self.rFreqs.calcFreqs(type = 'nuc')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=nuc.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=nuc.")

 

    ############################# by=nuc tests ##################################
    def test_RandFreqs_calcFreqs_bynuc_typenuc(self):
        correct_len = 4
        self.rFreqs = RandFreqs( by = 'nuc'  )
        freqs = self.rFreqs.calcFreqs(type = 'nuc')
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=nuc.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=nuc, type=nuc.")

 
 
 
 
 
 
 
 
 
 

class stateFreqs_EqualFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the EqualFreqs subclass of StateFreqs.'''
    
    def setUp(self):
        self.dec = 8 # For accuracy
    
    
    ############################# by=codon tests ##################################  
    def test_EqualFreqs_calcFreqs_bycodon_typecodon_restrict(self):
        correct = np.zeros(61)
        correct[10] = 0.25
        correct[23] = 0.25
        correct[38] = 0.25
        correct[54] = 0.25
        self.eqFreqs = EqualFreqs( by = 'codon', restrict = ['AGG', 'CCT', 'GCG', 'TGC'])
        freqs = self.eqFreqs.calcFreqs(type = 'codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=codon with restriction.")
    
    def test_EqualFreqs_calcFreqs_bycodon_typecodon(self):
        correct = np.array(np.repeat(1./61., 61))
        self.eqFreqs = EqualFreqs( by = 'codon' )
        freqs = self.eqFreqs.calcFreqs(type = 'codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=codon.")
    
    def test_EqualFreqs_calcFreqs_bycodon_typeamino(self):
        correct = np.array([0.0655737704918, 0.0327868852459, 0.0327868852459, 0.0327868852459, 0.0327868852459, 0.0655737704918, 0.0327868852459, 0.0491803278689, 0.0327868852459, 0.0983606557377, 0.016393442623, 0.0327868852459, 0.0655737704918, 0.0327868852459, 0.0983606557377, 0.0983606557377, 0.0655737704918, 0.0655737704918, 0.016393442623, 0.0327868852459])
        self.eqFreqs = EqualFreqs( by = 'codon' )
        freqs = self.eqFreqs.calcFreqs(type = 'amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=amino.")
    
    def test_EqualFreqs_calcFreqs_bycodon_typenuc(self):
        correct = np.array([0.24043715846994534, 0.26229508196721313, 0.25136612021857924, 0.2459016393442623])
        self.eqFreqs = EqualFreqs( by = 'codon' )
        freqs = self.eqFreqs.calcFreqs(type = 'nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=nuc.")

    ############################# by=amino tests ##################################
    def test_EqualFreqs_calcFreqs_byamino_typeamino_restrict(self):
        correct = [1./3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1./3., 0., 1./3., 0., 0., 0., 0., 0., 0., 0.]
        self.eqFreqs = EqualFreqs( by = 'amino', restrict = ['A', 'M', 'P'])
        freqs = self.eqFreqs.calcFreqs(type = 'amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=amino with restriction.")
    
    
    def test_EqualFreqs_calcFreqs_byamino_typecodon(self):
        correct = np.array([0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.01666667, 0.01666667, 0.05, 0.01666667, 0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.025, 0.025, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.025, 0.05, 0.025, 0.00833333, 0.025, 0.00833333, 0.025])
        self.eqFreqs = EqualFreqs( by = 'amino' )
        freqs = self.eqFreqs.calcFreqs( type = 'codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=codon.")
    
    def test_EqualFreqs_calcFreqs_byamino_typeamino(self):
        correct = np.array(np.repeat(1./20., 20))
        self.eqFreqs = EqualFreqs( by = 'amino'  )
        freqs = self.eqFreqs.calcFreqs(type = 'amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=amino.")
    
    def test_EqualFreqs_calcFreqs_byamino_typenuc(self):
        correct = [ 0.27638889,  0.22083333,  0.24861111,  0.25416667]
        self.eqFreqs = EqualFreqs( by = 'amino' )
        freqs = self.eqFreqs.calcFreqs(type = 'nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=nuc.")
    
    
    ############################# by=nuc tests ##################################
    def test_EqualFreqs_calcFreqs_bynuc_typenuc(self):
        correct = np.array(np.repeat(0.25, 4))
        self.eqFreqs = EqualFreqs( by = 'nuc' )
        freqs = self.eqFreqs.calcFreqs(type = 'nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=nuc, type=nuc.")




class stateFreqs_UserFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the UserFreqs subclass of StateFreqs.'''
    
    def setUp(self):
        self.dec = 8 # For accuracy

  
    ############################## by = nuc #############################################
    def test_UserFreqs_calcFreqs_bynuc_typenuc_singlenuc(self):
        correct = np.zeros(4)
        correct[1] = 1.0
        self.uFreqs = UserFreqs(by='nuc', freqs = {'C':1.0})
        freqs = self.uFreqs.calcFreqs(type='nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly by=nuc, type=nuc with single nucleotide freq specified.")
    
    def test_UserFreqs_calcFreqs_bynuc_typenuc_multiplenuc(self):
        correct = np.zeros(4)
        correct[0] = 0.25
        correct[1] = 0.25
        correct[2] = 0.25
        correct[3] = 0.25
        self.uFreqs = UserFreqs(by='nuc', freqs = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25})
        freqs = self.uFreqs.calcFreqs(type='nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly by=nuc, type=nuc with multiple nucleotide freqs specified.")
    
    ################################ by = codon ##########################################   
    def test_UserFreqs_calcFreqs_bycodon_typecodon_singlecodon(self):
        correct = np.zeros(61)
        correct[1] = 1.0
        self.uFreqs = UserFreqs(by='codon', freqs = {'AAC':1.0})
        freqs = self.uFreqs.calcFreqs( type='codon' )
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=codon, single codon freq specified.")
    
    def test_UserFreqs_calcFreqs_bycodon_typecodon_multiplecodons(self):
        correct = np.zeros(61)
        correct[1] = 0.5
        correct[2] = 0.5
        self.uFreqs = UserFreqs(by='codon', freqs = {'AAC':0.5, 'AAG':0.5})
        freqs = self.uFreqs.calcFreqs(type='codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=codon, multiple codon freqs specified.")
    
    def test_UserFreqs_calcFreqs_bycodon_typeamino(self):
        correct = np.zeros(20)
        correct[8] = 0.5
        correct[9] = 0.25
        correct[10] = 0.25
        self.uFreqs = UserFreqs(by = 'codon', freqs = {'AAA':0.5, 'ATG':0.25, 'CTT':0.25})
        freqs = self.uFreqs.calcFreqs(type = 'amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=amino, multiple codon freqs specified.")
       
       
       
    ################################ by = amino ##########################################
    def test_UserFreqs_calcFreqs_byamino_typeamino_singleaa(self):
        correct = np.zeros(20)
        correct[1] = 1.0
        self.uFreqs = UserFreqs(by='amino', freqs = {'C':1.0})
        freqs = self.uFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=amino with single aa specified.")
    
    def test_UserFreqs_calcFreqs_byamino_typeamino_multipleaa(self):
        correct = np.zeros(20)
        correct[1] = 0.25
        correct[2] = 0.25
        correct[3] = 0.25
        correct[4] = 0.25
        self.uFreqs = UserFreqs(by='amino', freqs = {'C':0.25, 'D':0.25, 'E':0.25, 'F':0.25})
        freqs = self.uFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=amino with multiple aas specified.")
    
    def test_UserFreqs_calcFreqs_byamino_typeamino_constraint(self):
        myFreqs = {'I': 0.33, 'L':0.33, 'V':0.34}
        correct = [0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.2805, 0.00882353, 0.2805, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.289, 0.00882353, 0.00882353]
        self.uFreqs = UserFreqs(by='amino', freqs = myFreqs, constraint = 0.85)
        freqs = self.uFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type = amino, with constraint.")
    
    def test_UserFreqs_calcFreqs_byamino_typecodon_multipleaa(self):
        correct = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.08333333, 0.08333333, 0.08333333, 0.08333333, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.125, 0.125, 0.125, 0.125, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.08333333, 0., 0.08333333, 0.]
        self.uFreqs = UserFreqs(by='amino', freqs = {'L':0.5, 'V':0.5})
        freqs = self.uFreqs.calcFreqs(type='codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=amino with multiple aas specified.")
    
    
  
class stateFreqs_ReadFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the ReadFreqs subclass of StateFreqs.'''
    
    def setUp(self):
        self.dec = 8 # For accuracy

    ################################# by = nuc ###########################################  
    def test_ReadFreqs_calcFreqs_bynuc_typenuc_nocol(self):
        correct = np.array([5./54., 12./54., 18./54., 19./54.])
        self.rFreqs = ReadFreqs(by='nuc', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=nuc, type=nuc, no columns.")

    def test_ReadFreqs_calcFreqs_bynuc_typenuc_col(self):
        correct = np.array([2./9., 7./9., 0, 0])
        self.rFreqs = ReadFreqs(by='nuc', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=nuc, type=nuc, with columns.")
    
    
    ################################## by = codon ########################################   
    def test_ReadFreqs_calcFreqs_bycodon_typecodon_nocol(self):
        correct = np.array([0, 0, 0, 0, 0, 1./18., 0, 0, 0, 0, 0, 0, 0, 0, 1./6., 0, 0, 0, 0, 0, 1./18., 1./18., 0, 0, 0, 0, 1./9., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1./6., 2./9., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1./6.])
        self.rFreqs = ReadFreqs(by='codon', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=codon, no columns.")


    def test_ReadFreqs_calcFreqs_bycodon_typecodon_col(self):
        correct = np.array([0, 0, 0, 0, 0, 1./9., 0, 0, 0, 0, 0, 0, 0, 0, 1./3., 0, 0, 0, 0, 0, 1./9., 1./9., 0, 0, 0, 0, 1./9., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2./9.])
        self.rFreqs = ReadFreqs(by='codon', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=codon, with columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typeamino_nocol(self):
        correct = np.array([0, 0, 0, 0, 1./6., 0, 0, 0, 0, 0, 1./6., 0, 1./9., 0, 1./9., 0, 1./18., 7./18., 0, 0])
        self.rFreqs = ReadFreqs(by='codon', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=amino, no columns.")
        
    def test_ReadFreqs_calcFreqs_bycodon_typeamino_col(self):
        correct = np.array([0, 0, 0, 0, 2./9., 0, 0, 0, 0, 0, 1./3., 0, 2./9., 0, 1./9., 0, 1./9., 0, 0, 0])
        self.rFreqs = ReadFreqs(by='codon', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=amino, with columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typenuc_nocol(self):
        correct = np.array([5./54., 12./54., 18./54., 19./54.])
        self.rFreqs = ReadFreqs(by='codon', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=nuc, no columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typenuc_col(self):
        correct = np.array([5./27., 8./27., 5./27., 9./27.])
        self.rFreqs = ReadFreqs(by='codon', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='nuc')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=nuc, with columns.")


    ################################## by = amino ########################################   
    def test_ReadFreqs_calcFreqs_byamino_typeamino_nocol(self):
        correct = np.array([1./6., 1./6., 1./6., 1./6., 1./6., 1./6., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.rFreqs = ReadFreqs(by='amino', file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=amino, no columns.")
    
    def test_ReadFreqs_calcFreqs_byamino_typeamino_col(self):
        correct = np.array([2./3., 0, 0, 1./3., 0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.rFreqs = ReadFreqs(by='amino', columns = [0,1,2], file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='amino')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=amino, with columns.")

    def test_ReadFreqs_calcFreqs_byamino_typecodon_nocol(self):
        correct = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08333333, 0.08333333, 0.08333333, 0.08333333, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08333333, 0, 0.08333333, 0, 0.08333333, 0, 0.08333333])
        self.rFreqs = ReadFreqs(by='amino', file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=codon, no columns.")

    def test_ReadFreqs_calcFreqs_byamino_typecodon_col(self):
        correct = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1./6., 0,  1./6., 0,  1./6.,  1./6.,  1./6., 1./6., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.rFreqs = ReadFreqs(by='amino', columns = [0,1,2], file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs(type='codon')
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=codon, with columns.")

    
if __name__ == '__main__':
    run_tests = unittest.TextTestRunner()
    
    print "Testing the RandFreqs subclass of StateFreqs"
    test_suite_Rand = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_RandFreqs_Tests)
    run_tests.run(test_suite_Rand)

    print "Testing the EqualFreqs subclass of StateFreqs"
    test_suite_Equal = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_EqualFreqs_Tests)
    run_tests.run(test_suite_Equal)

    print "Testing the UserFreqs subclass of StateFreqs"
    test_suite_User = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_UserFreqs_Tests)
    run_tests.run(test_suite_User)

    print "Testing the ReadFreqs subclass of StateFreqs"
    test_suite_Read = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_ReadFreqs_Tests)
    run_tests.run(test_suite_Read)

    
    
    
    
    
    
    
    