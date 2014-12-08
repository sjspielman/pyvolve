#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    Test model module.
'''


import unittest
#from pyvolve import *
import numpy as np
ZERO    = 1e-8
DECIMAL = 8


import sys
sys.path.append("/Users/sjspielman/Research/pyvolve/src/")
from genetics import *
from models import *
from partition import *
from newick import *
from state_freqs import *
from matrix_builder import *
from evolver import *

class model_nohet_tests(unittest.TestCase):
    ''' 
        Suite of tests for Model without site heterogeneity.
    '''
    
    def setUp(self):
        
        mu_dict     = {'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}
        nuc_freqs   = np.repeat(0.25, 4)
        codon_freqs = np.repeat(1./61., 61)
        amino_freqs = np.repeat(0.05, 20)
        
        self.nuc_model    = Model( {'state_freqs':nuc_freqs, 'mu':mu_dict}, "nucleotide")
        self.aa_model     = Model( {'state_freqs':amino_freqs, 'aa_model':'wag'}, "amino_acid")
        self.gy_model     = Model( {'state_freqs':codon_freqs, 'mu':mu_dict, 'beta':2.5, 'alpha':1.0}, "GY94")
        self.codon_model  = Model( {'state_freqs':codon_freqs, 'mu':mu_dict, 'beta':2.5, 'alpha':1.0}, "codon")
        self.mg_model     = Model( {'state_freqs':codon_freqs, 'mu':mu_dict, 'beta':2.5, 'alpha':1.0}, "MG94")
        self.mutsel_model = Model( {'state_freqs':codon_freqs, 'mu':mu_dict}, "mutsel")
   
        self.nuc_model.construct_model()
        self.aa_model.construct_model()
        self.gy_model.construct_model()
        self.codon_model.construct_model()
        self.mg_model.construct_model()
        self.mutsel_model.construct_model()
    
    def test_model_nohet_construct_model_matrix_type(self):
        '''
            Are matrices of correct dimension created?"
        '''
        self.assertTrue( self.nuc_model.matrix.shape == (4,4), msg = "nucleotide model without heterogeneity built wrong matrix type.")
        self.assertTrue( self.aa_model.matrix.shape == (20,20), msg = "amino acid model without heterogeneity built wrong matrix type.")
        self.assertTrue( self.codon_model.matrix.shape == (61,61), msg = "codon model without heterogeneity built wrong matrix type.")
        self.assertTrue( self.gy_model.matrix.shape == (61,61), msg = "gy model without heterogeneity built wrong matrix type.")
        self.assertTrue( self.mg_model.matrix.shape == (61,61), msg = "mg model without heterogeneity built wrong matrix type.")
        self.assertTrue( self.mutsel_model.matrix.shape == (61,61), msg = "mutset model without heterogeneity built wrong matrix type.")


    def test_model_nohet_rates(self):
        '''
            Are matrices of correct dimension created?"
        '''
        np.testing.assert_array_almost_equal(self.nuc_model.rate_factors, np.array([1.]), decimal = DECIMAL, err_msg = "rate_factor doesn't==1. for no heterogeneity.")
        np.testing.assert_array_almost_equal(self.nuc_model.rate_probs, np.array([1.]), decimal = DECIMAL, err_msg = "rate_probs doesn't==1. for no heterogeneity.")
        
        
 
 
 
 
class model_gammahet_tests(unittest.TestCase):
    ''' 
        Suite of tests for Model without gamma heterogeneity.
        Tests conducted with nucleotides.
    ''' 
    
    def setUp(self):
        
        mu_dict        = {'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}
        nuc_freqs      = np.repeat(0.25, 4)
        self.nuc_model = Model( {'state_freqs':nuc_freqs, 'mu':mu_dict}, "nucleotide")
    
    
    def test_model_het_gamma_rates_simprobs(self):
        '''
            Are gamma rates and probabilities assigned correctly?"
        '''    
        self.nuc_model.construct_model(alpha = 0.5, num_categories = 4)

        self.assertTrue( len(self.nuc_model.rate_probs) == 4, msg = "incorrect number of rate probabilities for gamma heterogeneity.")
        self.assertTrue( len(self.nuc_model.rate_factors) == 4, msg = "incorrect number of rate factors for gamma heterogeneity.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs)) < ZERO, msg = "rate probabilities don't sum to 1 for gamma hetereogenity.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs * self.nuc_model.rate_factors)) < ZERO, msg = "rate probabilities and factors improperly normalized for gamma hetereogenity.")
        
        
        
    def test_model_het_gamma_rates_userprobs(self):
        '''
            Are gamma rates and probabilities assigned correctly when users provide rate probabilities?"
        ''' 
        self.nuc_model.construct_model(alpha = 0.5, rate_probs = [0.25, 0.25, 0.3, 0.2])

        self.assertTrue( len(self.nuc_model.rate_probs) == 4, msg = "incorrect number of rate probabilities for gamma heterogeneity with user-provided probabilties.")
        self.assertTrue( len(self.nuc_model.rate_factors) == 4, msg = "incorrect number of rate factors for gamma heterogeneity with user-provided probabilties.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs)) < ZERO, msg = "rate probabilities don't sum to 1 for gamma hetereogenity with user-provided probabilties.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs * self.nuc_model.rate_factors)) < ZERO, msg = "rate probabilities and factors improperly normalized for gamma hetereogenity with user-provided probabilties.")

 




class model_userhet_tests(unittest.TestCase):
    ''' 
        Suite of tests for Model with user-specified heterogeneity.
        Tests conducted with nucleotides.
    ''' 
    
    def setUp(self):
        
        mu_dict        = {'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}
        nuc_freqs      = np.repeat(0.25, 4)
        self.nuc_model = Model( {'state_freqs':nuc_freqs, 'mu':mu_dict}, "nucleotide")
        self.rate_factors = [2.5, 4.5, 3.5, 0.8]
        self.rate_probs = [0.1, 0.1, 0.1, 0.7]
    
    
    
    def test_model_het_userrates_simprobs(self):
        '''
            Are rates and probabilities assigned correctly when users give rates but not probs?"
        '''    
        self.nuc_model.construct_model(rate_factors = self.rate_factors)
        
        self.assertTrue( len(self.nuc_model.rate_probs) == 4, msg = "incorrect number of rate probabilities for user-specified factors, default rates.")
        self.assertTrue( len(self.nuc_model.rate_factors) == 4, msg = "incorrect number of rate factors for user-specified factors, default rates.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs)) < ZERO, msg = "rate probabilities don't sum to 1 for user-specified factors, default rates.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs * self.nuc_model.rate_factors)) < ZERO, msg = "rate probabilities and factors improperly normalized for user-specified factors, default rates.")

    
    
    def test_model_het_userrates_userprobs(self):
        '''
            Are rates and probabilities assigned correctly when users give compatible rates and probs?"
        '''    
        
        self.nuc_model.construct_model(rate_probs = self.rate_probs, rate_factors = self.rate_factors)
        
        np.testing.assert_array_almost_equal( self.nuc_model.rate_probs, np.array(self.rate_probs), decimal = DECIMAL, err_msg = "incorrect rate probs from user-specified list.")               
        self.assertTrue( len(self.nuc_model.rate_probs) == 4, msg = "incorrect number of rate probabilities for user-specified factors, default rates.")
        self.assertTrue( len(self.nuc_model.rate_factors) == 4, msg = "incorrect number of rate factors for user-specified factors, default rates.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs)) < ZERO, msg = "rate probabilities don't sum to 1 for user-specified factors, default rates.")
        self.assertTrue( abs(1. - np.sum(self.nuc_model.rate_probs * self.nuc_model.rate_factors)) < ZERO, msg = "rate probabilities and factors improperly normalized for user-specified factors, default rates.")


    def test_model_het_bad_ratesprobs(self):
        '''
            Errors thrown when rates and probabilities specified incorrectly?"
        '''    
        self.assertRaises(AssertionError, self.nuc_model.construct_model(), rate_probs = [0.6, 0.3, 0.2], rate_factors = [1., 2, 3.], msg = "Assertion not raised when user-specified Model rate_probs sum > 1.")
        self.assertRaises(AssertionError, self.nuc_model.construct_model(), rate_probs = [0.6, 0.3, 0.1], rate_factors = [1., 2.], msg = "Assertion not raised when user-specified Model rate_probs diff size from rate_factors.")
 
 









class model_codonmodel_tests(unittest.TestCase):
    ''' 
        Suite of tests for CodonModel with user-specified heterogeneity.
        Tests conducted with GY.
    ''' 
    
    
    def setUp(self):
        
        mu_dict       = {'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}
        codon_freqs   = np.repeat(1./61., 61)
        self.gy_model = CodonModel( {'state_freqs':codon_freqs, 'mu':mu_dict, 'beta':[2.5, 1.5], 'alpha':[1.0, 0.75]}, "GY94")


    def test_codonmodel_simprobs(self):
        '''
            CodonModel created properly when probs are default?"
        '''    
        self.gy_model.construct_model()
        
        np.testing.assert_array_almost_equal(self.gy_model.rate_probs, np.array([0.5, 0.5]), decimal=DECIMAL, err_msg = "incorrect default rate_probs for CodonModel.")
        self.assertTrue( len(self.gy_model.matrices) == 2, msg = "incorrect number of matrices created for CodonModel.")
        self.assertTrue( self.gy_model.matrices[0].shape == (61,61) and self.gy_model.matrices[1].shape == (61,61), msg = "incorrect matrix dimensions for CodonModel.")


    def test_codonmodel_userprobs(self):
        '''
            CodonModel created properly when probs specified correctly?"
        '''    
        self.gy_model.construct_model(rate_probs = [0.9, 0.1])
        
        np.testing.assert_array_almost_equal(self.gy_model.rate_probs, np.array([0.9, 0.1]), decimal=DECIMAL, err_msg = "incorrect default rate_probs for CodonModel.")
        self.assertTrue( len(self.gy_model.matrices) == 2, msg = "incorrect number of matrices created for CodonModel.")
        self.assertTrue( self.gy_model.matrices[0].shape == (61,61) and self.gy_model.matrices[1].shape == (61,61), msg = "incorrect matrix dimensions for CodonModel.")
        

    def test_codonmodel_badprobs(self):
        '''
            CodonModel created properly when probs specified incorrectly?"
        '''    
        self.assertRaises(AssertionError, self.gy_model.construct_model(), rate_probs = [0.6, 0.9], msg = "Assertion not raised when user-specified CodonModel rate_probs sum > 1.")
        self.assertRaises(AssertionError, self.gy_model.construct_model(), rate_probs = [0.5, 0.25, 0.5], msg = "Assertion not raised when user-specified CodonModel rate_probs size diff from number of dN/dS values.")
    
        
def run_models_test():
       
    run_tests = unittest.TextTestRunner()
    
    print "Testing Model construction without site heterogeneity."
    test_suite_call = unittest.TestLoader().loadTestsFromTestCase(model_nohet_tests)
    run_tests.run(test_suite_call)
    
    print "Testing Model construction with gamma site heterogeneity."
    test_suite_call = unittest.TestLoader().loadTestsFromTestCase(model_gammahet_tests)
    run_tests.run(test_suite_call)

    print "Testing Model construction with user-specified site heterogeneity."
    test_suite_call = unittest.TestLoader().loadTestsFromTestCase(model_userhet_tests)
    run_tests.run(test_suite_call)
         
    print "Testing CodonModel construction."
    test_suite_call = unittest.TestLoader().loadTestsFromTestCase(model_codonmodel_tests)
    run_tests.run(test_suite_call)
        
         
        
        
        
        
        
        
        
        
