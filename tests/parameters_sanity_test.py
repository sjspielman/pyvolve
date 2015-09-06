#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

''' Suite of unit tests for parameters_sanity module.'''

import unittest
import numpy as np
from pyvolve import *
MOLECULES = Genetics()

ZERO=1e-8
DECIMAL=8


class ParametersSanity_sanity_freqs(unittest.TestCase):
    ''' 
        Set of unittests for state frequency sanity testing.
    '''
    
    
    def test_sanity_state_freqs_incorrect_size(self):    
        params = {'state_freqs':np.repeat(0.05, 20)}
        pm = Nucleotide_Sanity("nucleotide", params, size = 4)
        self.assertRaises(AssertionError, lambda: pm._sanity_state_freqs())
        
    def test_sanity_state_freqs_missing(self):    
        params = {}
        pm = Nucleotide_Sanity("nucleotide", params, size = 4)
        pm._sanity_state_freqs()
        np.testing.assert_array_almost_equal(pm.params['state_freqs'], np.array([0.25, 0.25, 0.25, 0.25]), decimal = DECIMAL, err_msg = "state_freqs not set to equal when not provided in params dict.")

    def test_sanity_state_freqs_incorrect_sum(self):    
        params = {"state_freqs": [0.1, 0.2, 0.1, 0.3]}
        pm = Nucleotide_Sanity("nucleotide", params, size = 4)
        self.assertRaises(AssertionError, lambda: pm._sanity_state_freqs())
        
    def test_sanity_state_freqs_ecm(self):
        ecmrest_freqs = [0.03000113, 0.02017108, 0.02634411, 0.02300609, 0.01552806, 0.02020108, 0.01214205, 0.01342405, 0.01037204, 0.01167905, 0.00819503, 0.01014204, 0.01355105, 0.02344109, 0.02010208, 0.0255761 , 0.01607606, 0.01170305, 0.02021108, 0.01109704, 0.01064204, 0.01010004, 0.01184305, 0.01000704, 0.00480002, 0.01414806, 0.00783703, 0.00831103, 0.00760003, 0.01738607, 0.02883912, 0.01446706, 0.03322313, 0.0245321 , 0.03187813, 0.02819811, 0.01590806, 0.02830711, 0.01885308, 0.01900508, 0.01579606, 0.02298209, 0.01019104, 0.01685207, 0.01090104, 0.01893808, 0.02274709, 0.01904708, 0.01578206, 0.01596506, 0.00975004, 0.01113104, 0.00895604, 0.01188005, 0.00702903, 0.01188005, 0.00602502, 0.01638707, 0.02138309, 0.01542506, 0.02210309] 
        params = {}    
        pm = ECM_Sanity("ecmrest", params, size = 61)
        pm._sanity_state_freqs(empirical = True)
        np.testing.assert_array_almost_equal(pm.params["state_freqs"], ecmrest_freqs, decimal = DECIMAL, err_msg = "ECMrest default state_freqs not set up properly.")

    def test_sanity_state_freqs_aminoacid(self):
        wag_freqs = [0.08662799999999994, 0.0193078, 0.0570451, 0.0580589, 0.0384319, 0.0832518, 0.0244313, 0.048466, 0.0620286, 0.086209, 0.0195027, 0.0390894, 0.0457631, 0.0367281, 0.043972, 0.0695179, 0.0610127, 0.0708956, 0.0143859, 0.0352742]
        params = {}    
        pm = ECM_Sanity("wag", params, size = 20)
        pm._sanity_state_freqs(empirical = True)
        np.testing.assert_array_almost_equal(pm.params["state_freqs"], wag_freqs, decimal = DECIMAL, err_msg = "WAG (aa model) default state_freqs not set up properly.")

    def test_sanity_state_freqs_MG_provided_codon(self):
        mg_nuc = [ 0.23577944,  0.2426595,   0.28132983,  0.24023122]
        mg_codon = [6.48685834e-03, 2.56846182e-02, 3.62821504e-02, 1.78992022e-02, 3.78995203e-03, 1.27376858e-02, 1.11733512e-02, 2.77350952e-03, 3.21357920e-02, 1.15713167e-02, 1.35629025e-03, 2.18752127e-02, 2.50168300e-05, 3.10876439e-02, 3.37328837e-02, 1.65327942e-02, 8.21428969e-03, 1.55572812e-02, 1.76293820e-02, 2.29747386e-02, 1.73289159e-03, 1.13379531e-02, 1.59673812e-02, 2.97921605e-02, 3.52280293e-02, 9.42724910e-03, 2.48052302e-02, 1.47583962e-02, 1.96410901e-02, 9.33795423e-03, 6.06200346e-05, 1.11892264e-02, 1.18844724e-02, 1.75173640e-02, 2.77504624e-02, 1.90563757e-02, 3.16335607e-02, 1.45931598e-02, 1.42666011e-02, 2.96971714e-02, 1.90421569e-02, 3.33563400e-03, 2.94568624e-03, 3.40917454e-02, 1.32157035e-02, 1.94067362e-02, 1.83063743e-02, 1.03705727e-02, 1.09380931e-02, 6.07488829e-03, 1.17862027e-02, 4.53122192e-03, 2.49775815e-02, 8.96308703e-03, 3.10326855e-02, 2.94072089e-02, 1.93306945e-02, 3.42785939e-03, 2.24745563e-02, 7.87119137e-03, 1.92728016e-02]
        params = {"state_freqs": mg_codon}
        pm = MechCodon_Sanity("mg", params, size=61)
        pm._sanity_state_freqs_MG()
        np.testing.assert_array_almost_equal(pm.params["nuc_freqs"], mg_nuc, decimal = DECIMAL, err_msg = "Nucleotide frequencies not properly calc'd from codon frequencies for an MG model.")
            



class ParametersSanity_sanity_mutation(unittest.TestCase):
    ''' 
        Set of unittests for mutation rate sanity testing.
    '''
    def test_sanity_mutation_rates_missing_reverse(self):    
        params = {'state_freqs':np.repeat(0.25, 4), 'mu': {'AC':1.,  'AG':2.5,  'AT':1., 'CG':1.,  'CT':2.5, 'GT':1.}}
        pm = Nucleotide_Sanity("nucleotide",params)
        pm._sanity_mutation_rates()
        self.assertTrue( pm.params['mu'] == {'AC':1., 'CA':1., 'AG':2.5, 'GA': 2.5, 'AT':1., 'TA':1., 'CG':1., 'GC':1., 'CT':2.5, 'TC':2.5, 'GT':1., 'TG':1.},  msg = "Backward mutation rates not properly assigned when all forward rates given.")

    def test_sanity_mutation_rates_missing_some(self):    
        params = {'state_freqs':np.repeat(0.25, 4), 'mu': {'AG':2.5,  'AT':1., 'CG':1.,  'CT':2.5, 'GT':1.}}
        pm = Nucleotide_Sanity("nucleotide",params)
        pm._sanity_mutation_rates()
        self.assertTrue( pm.params['mu'] == {'AC':1., 'CA':1., 'AG':2.5, 'GA': 2.5, 'AT':1., 'TA':1., 'CG':1., 'GC':1., 'CT':2.5, 'TC':2.5, 'GT':1., 'TG':1.},  msg = "Missing mutation rates not filled in when all but 1 pair provided.")

    def test_sanity_mutation_rates_kappa(self):    
        params = {'state_freqs':np.repeat(0.25, 4), 'kappa': 2.5}
        pm = Nucleotide_Sanity("nucleotide",params)
        pm._sanity_mutation_rates()
        self.assertTrue( pm.params['mu'] == {'AC':1., 'CA':1., 'AG':2.5, 'GA': 2.5, 'AT':1., 'TA':1., 'CG':1., 'GC':1., 'CT':2.5, 'TC':2.5, 'GT':1., 'TG':1.},  msg = "Missing mutation rates wrong when only kappa provided.")

    def test_sanity_mutation_rates_all_missing(self):    
        params = {'state_freqs':np.repeat(0.25, 4)}
        pm = Nucleotide_Sanity("nucleotide",params)
        pm._sanity_mutation_rates()
        self.assertTrue( pm.params['mu'] == {'AC':1., 'CA':1., 'AG':1., 'GA':1., 'AT':1., 'TA':1., 'CG':1., 'GC':1., 'CT':1., 'TC':1., 'GT':1., 'TG':1.},  msg = "Missing mutation rates not added in when none provided.")



class ParametersSanity_sanity_dnds(unittest.TestCase):
    ''' 
        Set of unittests for dN/dS rate sanity testing.
    '''
    
    def test_sanity_MechCodon_homo_dNdSmissing(self):     
        params = {}
        pm = MechCodon_Sanity("gy", params)
        self.assertRaises(KeyError, lambda: pm._sanity_dnds())

    def test_sanity_MechCodon_homo_dNmissing(self):
        params = {'alpha':1.5}
        pm = MechCodon_Sanity("gy", params)
        self.assertRaises(KeyError,lambda:  pm._sanity_dnds())
        
    def test_sanity_MechCodon_homo_dSmissing(self):
        params = {'beta':1.5}
        pm = MechCodon_Sanity("gy", params)
        pm._sanity_dnds()
        self.assertEqual(pm.params["alpha"], 1., msg = "dS (alpha key) not set to 1 when not provided.")

    def test_sanity_MechCodon_homo_omega(self):     
        params = {'omega':0.1}
        pm = MechCodon_Sanity("gy", params)
        pm._sanity_dnds()
        self.assertEqual(pm.params["beta"], 0.1, msg = "'omega' key not converted to 'beta' key in codon model sanity.")

    def test_sanity_MechCodon_hetero_dSmissing(self):     
        params = {'beta':[0.2, 0.3, 0.4]}
        pm = MechCodon_Sanity("gy", params)
        pm.hetmodel = True
        pm._sanity_dnds()
        np.testing.assert_array_almost_equal(pm.params["alpha"], np.ones(3), decimal = DECIMAL, err_msg = "dS not properly filled in when not provided for a heterogeneous codon model.")

    def test_sanity_MechCodon_hetero_bad_dNdS_combo(self):     
        params = {'beta':[0.2, 0.3, 0.4], 'alpha':9.}
        pm = MechCodon_Sanity("gy", params)
        pm.hetmodel = True
        self.assertRaises(TypeError, lambda: pm._sanity_dnds())


class ParametersSanity_ECM(unittest.TestCase):
    ''' 
        Unittest for ECM_Sanity.
    '''

    def test_sanity_ECM(self):
        newp = ECM_Sanity("ecmrest", {}, size = 61)()
        self.assertTrue(newp['beta'] == 1. and newp['alpha'] == 1., msg = "dN, dS not properly filled in as 1 when missing from params dictionary for ECM model.")
        self.assertTrue(newp['k_ti'] == 1. and newp['k_tv'] == 1., msg = "kappa params not properly initialized to 1. when not provided for ECM model.")
        

class ParametersSanity_MutSel(unittest.TestCase):
    ''' 
        Unittest for MutSel_Sanity.
    '''

    def test_sanity_MutSel_aa_fitness(self):
        aa_fitness = np.array([-1.309, -0.881,  0.662, -1.245,  3.909, 3.276,  3.578,  0.785, -0.957,  2.862, 2.919,  2.174,  2.907,  1.438,  3.231,0.733,  3.098,  2.505,  1.470, 0.])
        codon_from_aa_fitness = np.array([-0.957, 2.174, -0.957, 2.174, 3.098, 3.098, 3.098, 3.098, 3.231, 0.733, 3.231, 0.733, 0.785, 0.785, 2.919, 0.785, 1.438, 3.578, 1.438, 3.578, 2.907, 2.907, 2.907, 2.907, 3.231, 3.231, 3.231, 3.231, 2.862, 2.862, 2.862, 2.862, -1.245, 0.662, -1.245, 0.662, -1.309, -1.309, -1.309, -1.309, 3.276, 3.276, 3.276, 3.276, 2.505, 2.505, 2.505, 2.505, 0.0, 0.0, 0.733, 0.733, 0.733, 0.733, -0.881, 1.47, -0.881, 2.862, 3.909, 2.862, 3.909])
        pm = MutSel_Sanity("mutsel", {"fitness": aa_fitness})
        pm._sanity_state_freqs_fitness()
        self.assertTrue(pm.params["codon_model"], msg = "MutSel model not recognized as codon model when amino acid fitnesses provided.")
        self.assertFalse(pm.params["calc_by_freqs"], msg = "MutSel model didnt set calc_by_freqs as False when amino acid fitnesses provided.")
        np.testing.assert_array_almost_equal(pm.params["fitness"], codon_from_aa_fitness, decimal = DECIMAL, err_msg = "Mutsel amino acid fitness not properly converted to codon fitness.")


    def test_sanity_MutSel_codon_fitness(self):
        codon_fitness = np.array([-0.957, 2.174, -0.957, 2.174, 3.098, 3.098, 3.098, 3.098, 3.231, 0.733, 3.231, 0.733, 0.785, 0.785, 2.919, 0.785, 1.438, 3.578, 1.438, 3.578, 2.907, 2.907, 2.907, 2.907, 3.231, 3.231, 3.231, 3.231, 2.862, 2.862, 2.862, 2.862, -1.245, 0.662, -1.245, 0.662, -1.309, -1.309, -1.309, -1.309, 3.276, 3.276, 3.276, 3.276, 2.505, 2.505, 2.505, 2.505, 0.0, 0.0, 0.733, 0.733, 0.733, 0.733, -0.881, 1.47, -0.881, 2.862, 3.909, 2.862, 3.909])
        pm = MutSel_Sanity("mutsel", {"fitness": codon_fitness})
        pm._sanity_state_freqs_fitness()
        self.assertTrue(pm.params["codon_model"], msg = "MutSel model not recognized as codon model when codon fitnesses provided.")
        self.assertFalse(pm.params["calc_by_freqs"], msg = "MutSel model didnt set calc_by_freqs as False when codon fitnesses provided.")


    def test_sanity_MutSel_codon_frequencies(self):
        codon_freqs = [0.01372832, 0.01412889, 0.01638048, 0.0139875, 0.01412889, 0.01454117, 0.01685846, 0.01439566, 0.01638048, 0.01685846, 0.01954503, 0.01668976, 0.0139875, 0.01439566, 0.01668976, 0.0142516, 0.01412889, 0.01454117, 0.01685846, 0.01439566, 0.01454117, 0.01496548, 0.01735039, 0.01481573, 0.01685846, 0.01735039, 0.02011536, 0.01717677, 0.01439566, 0.01481573, 0.01717677, 0.01466747, 0.01638048, 0.01685846, 0.01954503, 0.01668976, 0.01685846, 0.01735039, 0.02011536, 0.01717677, 0.01954503, 0.02011536, 0.02332095, 0.01991406, 0.01668976, 0.01717677, 0.01991406, 0.01700488, 0.01439566, 0.0142516, 0.01439566, 0.01481573, 0.01717677, 0.01466747, 0.01717677, 0.01991406, 0.01700488, 0.0142516, 0.01466747, 0.01700488, 0.01452069]
        pm = MutSel_Sanity("mutsel", {"state_freqs": codon_freqs})
        pm._sanity_state_freqs_fitness()
        self.assertTrue(pm.params["codon_model"], msg = "MutSel model not recognized as codon model when codon frequencies provided.")
        self.assertTrue(pm.params["calc_by_freqs"], msg = "MutSel model didnt set calc_by_freqs as True when codon frequencies provided.")


    def test_sanity_MutSel_nuc_fitness(self):
        nuc_fitness = [0.2, 0.4, 1.6, 1.5]
        pm = MutSel_Sanity("mutsel", {"fitness": nuc_fitness})
        pm._sanity_state_freqs_fitness()
        self.assertFalse(pm.params["codon_model"], msg = "MutSel model recognized as codon model when nucleotide fitness provided.")
        self.assertFalse(pm.params["calc_by_freqs"], msg = "MutSel model didnt set calc_by_freqs as False when nucleotide fitnss provided.")


    def test_sanity_MutSel_nuc_frequencies(self):
        nuc_freqs = [0.25, 0.25, 0.25, 0.25]
        pm = MutSel_Sanity("mutsel", {"state_freqs": nuc_freqs})
        pm._sanity_state_freqs_fitness()
        self.assertFalse(pm.params["codon_model"], msg = "MutSel model not recognized as codon model when nucleotide frequencies provided.")
        self.assertTrue(pm.params["calc_by_freqs"], msg = "MutSel model didnt set calc_by_freqs as True when nucleotide frequencies provided.")
        









