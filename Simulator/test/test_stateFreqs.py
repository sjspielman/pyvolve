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
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=nuc.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=codon, type=codon.")
    
    def test_RandFreqs_calcFreqs_bycodon_typeposNuc(self):
        correct_shape = (3,4)
        self.rFreqs = RandFreqs( by = 'codon', type = 'posNuc' )
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_almost_equal(np.sum(freqs, axis = 1), ([1., 1., 1.]), decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=codon, type=posNuc.")
        self.assertEqual(freqs.shape, correct_shape, msg= "RandFreqs has incorrect size for by=codon, type=posNuc.")

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
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=nuc.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=amino, type=nuc.")

    def test_RandFreqs_calcFreqs_byamino_typeposNuc(self):
        correct_shape = (3,4)
        self.rFreqs = RandFreqs( by = 'amino', type = 'posNuc' )
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_almost_equal(np.sum(freqs, axis = 1), ([1., 1., 1.]), decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=amino, type=posNuc.")
        self.assertEqual(freqs.shape, correct_shape, msg= "RandFreqs has incorrect size for by=amino, type=posNuc.")


    ############################# by=nuc tests ##################################
    def test_RandFreqs_calcFreqs_bynuc_typenuc(self):
        correct_len = 4
        self.rFreqs = RandFreqs( by = 'nuc', type = 'nuc' )
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_almost_equal(np.sum(freqs), 1., decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=nuc.")
        self.assertEqual(len(freqs), correct_len, msg= "RandFreqs has incorrect size for by=nuc, type=nuc.")

    def test_RandFreqs_calcFreqs_bynuc_typeposNuc(self):
        correct_shape = (3,4)
        self.rFreqs = RandFreqs( by = 'nuc', type = 'posNuc' )
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_almost_equal(np.sum(freqs, axis = 1), ([1., 1., 1.]), decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=nuc.")
        self.assertEqual(freqs.shape, correct_shape, msg= "RandFreqs has incorrect size for by=nuc, type=posNuc.")

    ############################# by=posNuc tests ##################################
    def test_RandFreqs_calcFreqs_byposNuc_typeposNuc(self):
        correct_shape = (3,4)
        self.rFreqs = RandFreqs( by = 'posNuc', type = 'posNuc' )
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_almost_equal(np.sum(freqs, axis=1), ([1., 1., 1.]), decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=posNuc, type=posNuc.")
        self.assertEqual(freqs.shape, correct_shape, msg= "RandFreqs has incorrect size for by=posNuc, type=posNuc.")

    def test_RandFreqs_calcFreqs_bynuc_typeposNuc(self):
        correct_shape = (3,4)
        self.rFreqs = RandFreqs( by = 'nuc', type = 'posNuc' )
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_almost_equal(np.sum(freqs, axis = 1), ([1., 1., 1.]), decimal = self.dec, err_msg = "RandFreqs do not sum to 1 for by=nuc, type=nuc.")
        self.assertEqual(freqs.shape, correct_shape, msg= "RandFreqs has incorrect size for by=nuc, type=posNuc.")


    


class stateFreqs_EqualFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the EqualFreqs subclass of StateFreqs.'''
    
    def setUp(self):
        self.dec = 8 # For accuracy
    
    ############################# by=posNuc tests ##################################
    def test_EqualFreqs_calcFreqs_byposNuc_typeposNuc(self):
        correct = [[ 0.25,   0.25,   0.25,  0.25 ],
                   [ 0.25,   0.25,   0.25,  0.25 ],
                   [ 0.25,   0.25,   0.25,  0.25 ]]
        self.eqFreqs = EqualFreqs( by = 'posNuc', type = 'posNuc' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=posNuc, type=posNuc.")

    
    ############################# by=codon tests ##################################
    def test_EqualFreqs_calcFreqs_bycodon_typeposNuc(self):
        correct = [[ 0.26229508,  0.26229508,  0.26229508,  0.21311475],
                   [ 0.2295082,   0.26229508,  0.24590164,  0.26229508],
                   [ 0.2295082,   0.26229508,  0.24590164,  0.26229508]]
        self.eqFreqs = EqualFreqs( by = 'codon', type = 'posNuc' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=codon, type=posNuc.")
        
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
    
    def test_EqualFreqs_calcFreqs_byamino_typeposNuc(self):
        correct = [[ 0.28333333,  0.21666667,  0.25,        0.25      ],
                   [ 0.35,        0.18333333,  0.21666667,  0.25      ],
                   [ 0.19583333,  0.2625,      0.27916667,  0.2625    ]]
        self.eqFreqs = EqualFreqs( by = 'amino', type = 'posNuc' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=posNuc.")
    
    def test_EqualFreqs_calcFreqs_byamino_typecodon(self):
        correct = np.array([0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.01666667, 0.01666667, 0.05, 0.01666667, 0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.025, 0.025, 0.025, 0.025, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.025, 0.025, 0.00833333, 0.00833333, 0.00833333, 0.00833333, 0.025, 0.05, 0.025, 0.00833333, 0.025, 0.00833333, 0.025])
        self.eqFreqs = EqualFreqs( by = 'amino', type = 'codon' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=codon.")
    
    def test_EqualFreqs_calcFreqs_byamino_typeamino(self):
        correct = np.array(np.repeat(1./20., 20))
        self.eqFreqs = EqualFreqs( by = 'amino', type = 'amino' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=amino.")
    
    def test_EqualFreqs_calcFreqs_byamino_typenuc(self):
        correct = [ 0.27638889,  0.22083333,  0.24861111,  0.25416667]
        self.eqFreqs = EqualFreqs( by = 'amino', type = 'nuc' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=amino, type=nuc.")
    
    
    ############################# by=nuc tests ##################################
    def test_EqualFreqs_calcFreqs_bynuc_typeposNuc(self):
        correct = [[ 0.25,   0.25,   0.25,  0.25 ],
                   [ 0.25,   0.25,   0.25,  0.25 ],
                   [ 0.25,   0.25,   0.25,  0.25 ]]
        self.eqFreqs = EqualFreqs( by = 'nuc', type = 'posNuc' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=nuc, type=posNuc.")
        
   
    def test_EqualFreqs_calcFreqs_bynuc_typenuc(self):
        correct = np.array(np.repeat(0.25, 4))
        self.eqFreqs = EqualFreqs( by = 'nuc', type = 'nuc' )
        freqs = self.eqFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "EqualFreqs not calculated properly for by=nuc, type=nuc.")




class stateFreqs_UserFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the UserFreqs subclass of StateFreqs.'''
    
    def setUp(self):
        self.dec = 8 # For accuracy

    ############################### by = posNuc ##########################################
    def test_UserFreqs_calcFreqs_byposNuc_typeposNuc(self):
        correct = [[ 0.5,   0.1,   0.25,  0.15],
                   [ 0.1,   0.4,   0.3,   0.2 ],
                   [ 0.2,   0.5,   0.05,  0.25]]
        myFreqs = {'A':[0.5, 0.1, 0.2], 'C':[0.1, 0.4, 0.5], 'G':[0.25, 0.3, 0.05], 'T':[0.15, 0.2, 0.25] }
        self.uFreqs = UserFreqs(by='posNuc', type='posNuc', freqs = myFreqs)
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=posNuc, type=posNuc.")
        
    ############################## by = nuc #############################################
    def test_UserFreqs_calcFreqs_bynuc_typeposNuc(self):
        correct = [[ 0.2,   0.3,   0.1,  0.4],
                   [ 0.2,   0.3,   0.1,  0.4],
                   [ 0.2,   0.3,   0.1,  0.4]]
        myFreqs = {'A':0.2, 'C':0.3, 'G':0.1, 'T':0.4}
        self.uFreqs = UserFreqs(by='nuc', type='posNuc', freqs = myFreqs)
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=nuc, type=posNuc.")
   
    def test_UserFreqs_calcFreqs_bynuc_typenuc_singlenuc(self):
        correct = np.zeros(4)
        correct[1] = 1.0
        self.uFreqs = UserFreqs(by='nuc', type='nuc', freqs = {'C':1.0})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly by=nuc, type=nuc with single nucleotide freq specified.")
    
    def test_UserFreqs_calcFreqs_bynuc_typenuc_multiplenuc(self):
        correct = np.zeros(4)
        correct[0] = 0.25
        correct[1] = 0.25
        correct[2] = 0.25
        correct[3] = 0.25
        self.uFreqs = UserFreqs(by='nuc', type='nuc', freqs = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly by=nuc, type=nuc with multiple nucleotide freqs specified.")
    
    ################################ by = codon ##########################################
    def test_UserFreqs_calcFreqs_bycodon_typeposNuc(self):
        correct = [[ 0.5,  0.25,   0.25, 0. ],
                   [ 0.,   0.25,   0.25, 0.5],
                   [ 0.25, 0.25,   0.5,  0. ]]
        myFreqs = {'ATG':0.25, 'ACC':0.25, 'CTA':0.25, 'GGG':0.25}
        self.uFreqs = UserFreqs(by='codon', type='posNuc', freqs = myFreqs)
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=posNuc.")
    
    def test_UserFreqs_calcFreqs_bycodon_typecodon_singlecodon(self):
        correct = np.zeros(61)
        correct[1] = 1.0
        self.uFreqs = UserFreqs(by='codon', type='codon', freqs = {'AAC':1.0})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=codon, single codon freq specified.")
    
    def test_UserFreqs_calcFreqs_bycodon_typecodon_multiplecodons(self):
        correct = np.zeros(61)
        correct[1] = 0.5
        correct[2] = 0.5
        self.uFreqs = UserFreqs(by='codon', type='codon', freqs = {'AAC':0.5, 'AAG':0.5})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=codon, multiple codon freqs specified.")
    
    def test_UserFreqs_calcFreqs_bycodon_typeamino(self):
        correct = np.zeros(20)
        correct[8] = 0.5
        correct[9] = 0.25
        correct[10] = 0.25
        self.uFreqs = UserFreqs(by = 'codon', type = 'amino', freqs = {'AAA':0.5, 'ATG':0.25, 'CTT':0.25})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=codon, type=amino, multiple codon freqs specified.")
       
       
       
    ################################ by = amino ##########################################
    def test_UserFreqs_calcFreqs_byamino_typeposNuc(self):
        correct = [[ 0.5,   0.5,   0.,    0. ],
                   [ 0.,    0.5,   0.,    0.5],
                   [ 0.125, 0.125, 0.625, 0.125 ]]
        self.uFreqs = UserFreqs(by='amino', type='posNuc', freqs = {'M':0.5, 'P':0.5})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=posNuc.")
    
    
    def test_UserFreqs_calcFreqs_byamino_typeamino_singleaa(self):
        correct = np.zeros(20)
        correct[1] = 1.0
        self.uFreqs = UserFreqs(by='amino', type='amino', freqs = {'C':1.0})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=amino with single aa specified.")
    
    def test_UserFreqs_calcFreqs_byamino_typeamino_multipleaa(self):
        correct = np.zeros(20)
        correct[1] = 0.25
        correct[2] = 0.25
        correct[3] = 0.25
        correct[4] = 0.25
        self.uFreqs = UserFreqs(by='amino', type='amino', freqs = {'C':0.25, 'D':0.25, 'E':0.25, 'F':0.25})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=amino with multiple aas specified.")
    
    def test_UserFreqs_calcFreqs_byamino_typeamino_constraint(self):
        myFreqs = {'I': 0.33, 'L':0.33, 'V':0.34}
        correct = [0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.2805, 0.00882353, 0.2805, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.00882353, 0.289, 0.00882353, 0.00882353]
        self.uFreqs = UserFreqs(by='amino', type='amino', freqs = myFreqs, constraint = 0.85)
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type = amino, with constraint.")
    
    def test_UserFreqs_calcFreqs_byamino_typecodon_multipleaa(self):
        correct = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.08333333, 0.08333333, 0.08333333, 0.08333333, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.125, 0.125, 0.125, 0.125, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.08333333, 0., 0.08333333, 0.]
        self.uFreqs = UserFreqs(by='amino', type='codon', freqs = {'L':0.5, 'V':0.5})
        freqs = self.uFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "UserFreqs not calculated properly for by=amino, type=amino with multiple aas specified.")
    
    
### NEED ###
# 2. type posNuc, by amino.   
class stateFreqs_ReadFreqs_Tests(unittest.TestCase):
    ''' Set of "unittests" for the ReadFreqs subclass of StateFreqs.'''
    
    def setUp(self):
        self.dec = 8 # For accuracy

    ############################# by = posNuc ############################################
    def test_ReadFreqs_calcFreqs_byposNuc_typeposNuc_nocol(self):
        correct = [ [ 0.22222222,  0.22222222,  0.38888889,  0.16666667],
                    [ 0.,          0.16666667,  0.11111111,  0.72222222],
                    [ 0.05555556,  0.27777778,  0.5,         0.16666667]]
        self.rFreqs = ReadFreqs(by = 'posNuc', type = 'posNuc', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=posNuc, type=posNuc, no columns.")
    
    def test_ReadFreqs_calcFreqs_byposNuc_typeposNuc_col(self):
        correct = [ [ 4./9.,   3./9.,   0.,    2./9.],
                    [ 0.,      3./9.,   1./9., 5./9.],
                    [ 1./9.,   2./9.,   4./9., 2./9.]]
        self.rFreqs = ReadFreqs(by = 'posNuc', type = 'posNuc', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=posNuc, type=posNuc, with columns.")
  
  
    ################################# by = nuc ###########################################  
    def test_ReadFreqs_calcFreqs_bynuc_typeposNuc_nocol(self):
        correct = [[ 0.09259259, 0.22222222, 0.33333333, 0.35185185],
                    [ 0.09259259, 0.22222222, 0.33333333, 0.35185185],
                    [ 0.09259259, 0.22222222, 0.33333333, 0.35185185]]
        self.rFreqs = ReadFreqs(by = 'nuc', type = 'posNuc', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=nuc, type=posNuc, no columns.")
    
    def test_ReadFreqs_calcFreqs_bynuc_typeposNuc_col(self):
        correct = [[ 2./9., 7./9., 0., 0. ],
                   [ 2./9., 7./9., 0., 0. ],
                   [ 2./9., 7./9., 0., 0. ]]
        self.rFreqs = ReadFreqs(by = 'nuc', type = 'posNuc', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=nuc, type=posNuc, no columns.")
    
    def test_ReadFreqs_calcFreqs_bynuc_typenuc_nocol(self):
        correct = np.array([5./54., 12./54., 18./54., 19./54.])
        self.rFreqs = ReadFreqs(by='nuc', type='nuc', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=nuc, type=nuc, no columns.")

    def test_ReadFreqs_calcFreqs_bynuc_typenuc_col(self):
        correct = np.array([2./9., 7./9., 0, 0])
        self.rFreqs = ReadFreqs(by='nuc', type='nuc', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=nuc, type=nuc, with columns.")
    ################################## by = codon ########################################   
    def test_ReadFreqs_calcFreqs_bycodon_typecodon_nocol(self):
        correct = np.array([0, 0, 0, 0, 0, 1./18., 0, 0, 0, 0, 0, 0, 0, 0, 1./6., 0, 0, 0, 0, 0, 1./18., 1./18., 0, 0, 0, 0, 1./9., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1./6., 2./9., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1./6.])
        self.rFreqs = ReadFreqs(by='codon', type='codon', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=codon, no columns.")


    def test_ReadFreqs_calcFreqs_bycodon_typecodon_col(self):
        correct = np.array([0, 0, 0, 0, 0, 1./9., 0, 0, 0, 0, 0, 0, 0, 0, 1./3., 0, 0, 0, 0, 0, 1./9., 1./9., 0, 0, 0, 0, 1./9., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2./9.])
        self.rFreqs = ReadFreqs(by='codon', type='codon', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=codon, with columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typeamino_nocol(self):
        correct = np.array([0, 0, 0, 0, 1./6., 0, 0, 0, 0, 0, 1./6., 0, 1./9., 0, 1./9., 0, 1./18., 7./18., 0, 0])
        self.rFreqs = ReadFreqs(by='codon', type='amino', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=amino, no columns.")
        
    def test_ReadFreqs_calcFreqs_bycodon_typeamino_col(self):
        correct = np.array([0, 0, 0, 0, 2./9., 0, 0, 0, 0, 0, 1./3., 0, 2./9., 0, 1./9., 0, 1./9., 0, 0, 0])
        self.rFreqs = ReadFreqs(by='codon', type='amino', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=amino, with columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typenuc_nocol(self):
        correct = np.array([5./54., 12./54., 18./54., 19./54.])
        self.rFreqs = ReadFreqs(by='codon', type='nuc', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=nuc, no columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typenuc_col(self):
        correct = np.array([5./27., 8./27., 5./27., 9./27.])
        self.rFreqs = ReadFreqs(by='codon', type='nuc', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=nuc, with columns.")

    def test_ReadFreqs_calcFreqs_bycodon_typeposNuc_nocol(self):
        correct = [ [ 0.22222222,  0.22222222,  0.38888889,  0.16666667],
                    [ 0.,          0.16666667,  0.11111111,  0.72222222],
                    [ 0.05555556,  0.27777778,  0.5,         0.16666667]]
        self.rFreqs = ReadFreqs(by = 'codon', type = 'posNuc', file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=posNuc, no columns.")
    
    def test_ReadFreqs_calcFreqs_bycodon_typeposNuc_col(self):
        correct = [ [ 4./9.,   3./9.,   0.,    2./9.],
                    [ 0.,      3./9.,   1./9., 5./9.],
                    [ 1./9.,   2./9.,   4./9., 2./9.]]
        self.rFreqs = ReadFreqs(by = 'codon', type = 'posNuc', columns = [0,1,2], file = 'freqFiles/testFreq_codon_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=codon, type=posNuc, with columns.")
     

    ################################## by = amino ########################################   
    def test_ReadFreqs_calcFreqs_byamino_typeamino_nocol(self):
        correct = np.array([1./6., 1./6., 1./6., 1./6., 1./6., 1./6., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.rFreqs = ReadFreqs(by='amino', type='amino', file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=amino, no columns.")
    
    def test_ReadFreqs_calcFreqs_byamino_typeamino_col(self):
        correct = np.array([2./3., 0, 0, 1./3., 0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.rFreqs = ReadFreqs(by='amino', type='amino', columns = [0,1,2], file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=amino, with columns.")

    def test_ReadFreqs_calcFreqs_byamino_typecodon_nocol(self):
        correct = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08333333, 0.08333333, 0.08333333, 0.08333333, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08333333, 0, 0.08333333, 0, 0.08333333, 0, 0.08333333])
        self.rFreqs = ReadFreqs(by='amino', type='codon', file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=codon, no columns.")

    def test_ReadFreqs_calcFreqs_byamino_typecodon_col(self):
        correct = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1./6., 0,  1./6., 0,  1./6.,  1./6.,  1./6., 1./6., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.rFreqs = ReadFreqs(by='amino', type='codon', columns = [0,1,2], file = 'freqFiles/testFreq_amino_aln.fasta')
        freqs = self.rFreqs.calcFreqs()
        np.testing.assert_array_almost_equal(correct, freqs, decimal = self.dec, err_msg = "ReadFreqs not calculated properly for by=amino, type=codon, with columns.")

    ########### CAN ALWAYS ADD IN SOME TESTING FOR BYAMINO TYPEPOSNUC, BUT THIS SEEMS LIKE IT WILL BE THE MOST UNUSED FUNCTIONALITY OF EVER. ##########################
    
if __name__ == '__main__':
    run_tests = unittest.TextTestRunner()
    

    print "Testing the EqualFreqs subclass of StateFreqs"
    test_suite_Equal = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_EqualFreqs_Tests)
    run_tests.run(test_suite_Equal)
    
    print "Testing the RandFreqs subclass of StateFreqs"
    test_suite_Rand = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_RandFreqs_Tests)
    run_tests.run(test_suite_Rand)
    
    print "Testing the UserFreqs subclass of StateFreqs"
    test_suite_User = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_UserFreqs_Tests)
    run_tests.run(test_suite_User)

    print "Testing the ReadFreqs subclass of StateFreqs"
    test_suite_Read = unittest.TestLoader().loadTestsFromTestCase(stateFreqs_ReadFreqs_Tests)
    run_tests.run(test_suite_Read)
    
    
    
    
    
    
    
    
    
    
    
    