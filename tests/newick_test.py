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


class newick_tests(unittest.TestCase):
    ''' 
        Suite of tests for newick tree parsing.
    '''
    
    def setUp(self):
        
        
        self.string_noflags = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921):0.207);"
        self.string_flags   = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921_m2_):0.207);"
        self.string_nodenames = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)bobby:0.762)bobbybubby:0.921)robert:0.207);"
        self.string_nodenames_flags = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)bobby:0.762_m1_)bobbybubby:0.921_m2_)robert:0.207);"
 
        with open('tests/newickFiles/printed_noflags.txt', 'r') as f:
            self.true_noflags = str(f.read())
        with open('tests/newickFiles/printed_flags.txt', 'r') as f:
            self.true_flags = str(f.read())
        with open('tests/newickFiles/printed_nodenames.txt', 'r') as f:
            self.true_nodenames = str(f.read())
        with open('tests/newickFiles/printed_nodenames_flags.txt', 'r') as f:
            self.true_nodenames_flags = str(f.read())    
    
    def tearDown(self):
        os.remove("out.txt")


       
        
    def test_newick_read_tree_file(self):
        '''
            Test newick reading in tree from file, no flags.
        '''
        t = read_tree(file = 'tests/newickFiles/test_tree.tre') 
 
        orig_stdout = sys.stdout
        f = file('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())

        self.assertMultiLineEqual(printed, self.true_noflags, msg = "Couldn't read and parse tree from file properly.")


    def test_newick_read_tree_string_noflags(self):
        '''
            Test newick reading in tree from string, no flags.
        '''
        t = read_tree(tree = self.string_noflags)
        
        orig_stdout = sys.stdout
        f = file('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        
        self.assertMultiLineEqual(printed, self.true_noflags, msg = "Couldn't read and parse tree from string properly, no flags.")



    def test_newick_read_tree_string_flags(self):
        '''
            Test newick reading in tree from string, with flags.
        '''
        t = read_tree(tree = self.string_flags)
        
        orig_stdout = sys.stdout
        f = file('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_flags, msg = "Couldn't read and parse tree from string properly, no flags.")



    def test_newick_read_tree_nodenames(self):
        '''
            Test newick reading in tree from string, with node names.
        '''
        t = read_tree(tree = self.string_nodenames)
        orig_stdout = sys.stdout
        f = file('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nodenames, msg = "Couldn't read and parse tree from string properly when node names provided.")


    def test_newick_read_tree_nodenames_flags(self):
        '''
            Test newick reading in tree from string, with node names and flags.
        '''
        t = read_tree(tree = self.string_nodenames_flags)
        orig_stdout = sys.stdout
        f = file('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nodenames_flags, msg = "Couldn't read and parse tree from string properly when node names and flags provided.")
    
        
        
        
        
        
        
        
        
        
        
        
        
        
