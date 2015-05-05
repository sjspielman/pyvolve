#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    Test newick module.
'''


import unittest
from pyvolve import *
import sys, io, traceback

def capture(func, *args, **kwargs):
    """Capture the output of func when called with the given arguments.

    The function output includes any exception raised. capture returns
    a tuple of (function result, standard output, standard error).
    """
    stdout, stderr = sys.stdout, sys.stderr
    sys.stdout = c1 = io.StringIO()
    sys.stderr = c2 = io.StringIO()
    result = None
    try:
        result = func(*args, **kwargs)
    except:
        traceback.print_exc()
    sys.stdout = stdout
    sys.stderr = stderr
    return (result, c1.getvalue(), c2.getvalue())



class newick_tests(unittest.TestCase):
    ''' 
        Suite of tests for newick tree parsing
    '''
    
    def setUp(self):
        
        
        self.string_noflags = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921):0.207);"
        self.string_flags   = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921_m2_):0.207);"
        self.no_bls         = "(t4,(t3,(t2,(t5,t1))));"
        
        with open('tests/newickFiles/printed_noflags.txt', 'r') as f:
            self.true_noflags = str(f.read())
        with open('tests/newickFiles/printed_flags.txt', 'r') as f:
            self.true_flags = str(f.read())
    
    def test_newick_read_tree_file(self):
        '''
            Test newick reading in tree from file, no flags.
        '''
        t = read_tree(file = 'tests/newickFiles/test_tree.tre')
        printed = capture( print_tree, t )[1].rstrip()
        self.assertMultiLineEqual(printed, self.true_noflags, msg = "Couldn't read and parse tree from file properly.")


    def test_newick_read_tree_string_noflags(self):
        '''
            Test newick reading in tree from string, no flags.
        '''
        t = read_tree(tree = self.string_noflags)
        printed = capture( print_tree, t )[1].rstrip()
        self.assertMultiLineEqual(printed, self.true_noflags, msg = "Couldn't read and parse tree from string properly, no flags.")



    def test_newick_read_tree_string_flags(self):
        '''
            Test newick reading in tree from string, with flags.
        '''
        t = read_tree(tree = self.string_flags)
        printed = capture( print_tree, t )[1].rstrip()
        self.assertMultiLineEqual(printed, self.true_flags, msg = "Couldn't read and parse tree from string properly, no flags.")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
