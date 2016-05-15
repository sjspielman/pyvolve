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
        self.string_propflags_underscore = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921_m2_):0.207);"
        self.string_propflags_hashtag   = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762#m1#):0.921#m2#):0.207);"

        self.string_nodenames = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)bobby:0.762)bobbybubby:0.921)robert:0.207);"
        self.string_nodenames_propflags = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)bobby:0.762_m1_)bobbybubby:0.921_m2_)robert:0.207);"
        self.string_nodenames_nopropflags = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)bobby:0.762_m1)bobbybubby:0.921_m2)robert:0.207);"

        self.string_nopropflags_underscore = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1):0.921_m2):0.207);"  
        self.string_nopropflags_hashtag = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762#m1):0.921#m2):0.207);" 
 
        with open('tests/newickFiles/printed_noflags.txt', 'r') as f:
            self.true_noflags = str(f.read())
        with open('tests/newickFiles/printed_noflags_scaled.txt', 'r') as f:
            self.true_noflags_scaled = str(f.read())
        with open('tests/newickFiles/printed_propflags.txt', 'r') as f:
            self.true_propflags = str(f.read())
        with open('tests/newickFiles/printed_nopropflags.txt', 'r') as f:
            self.true_nopropflags = str(f.read())    
        with open('tests/newickFiles/printed_nodenames.txt', 'r') as f:
            self.true_nodenames = str(f.read())
        with open('tests/newickFiles/printed_nodenames_propflags.txt', 'r') as f:
            self.true_nodenames_propflags = str(f.read())      

        with open('tests/newickFiles/printed_nodenames_nopropflags.txt', 'r') as f:
            self.true_nodenames_nopropflags = str(f.read())    


    
    def tearDown(self):
        try:
            os.remove("out.txt")
        except:
            pass

       
        
    def test_newick_read_tree_open(self):
        '''
            Test newick reading in tree from file, no flags.
        '''
        t = read_tree(file = 'tests/newickFiles/test_tree.tre') 
 
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())

        self.assertMultiLineEqual(printed, self.true_noflags, msg = "Couldn't read and parse tree from file properly.")


    def test_newick_read_tree_file_scaletree(self):
        '''
            Test newick reading in tree from file, no flags, with scale_tree.
        '''
        t = read_tree(file = 'tests/newickFiles/test_tree.tre', scale_tree = 10.) 
 
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())

        self.assertMultiLineEqual(printed, self.true_noflags_scaled, msg = "Couldn't read and parse tree from file properly.")



    def test_newick_read_tree_string_noflags(self):
        '''
            Test newick reading in tree from string, no flags.
        '''
        t = read_tree(tree = self.string_noflags)
        
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        
        self.assertMultiLineEqual(printed, self.true_noflags, msg = "Couldn't read and parse tree from string properly, no flags.")



    def test_newick_read_tree_string_flags_prop_underscore(self):
        '''
            Test newick reading in tree from string, with propagating underscore flags.
        '''
        t = read_tree(tree = self.string_propflags_underscore)
        
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_propflags, msg = "Couldn't read and parse tree from string properly, propagating underscore flags.")


    def test_newick_read_tree_string_flags_prop_hashtag(self):
        '''
            Test newick reading in tree from string, with propagating hashtag flags.
        '''
        t = read_tree(tree = self.string_propflags_hashtag)
        
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_propflags, msg = "Couldn't read and parse tree from string properly, propagating hashtag flags.")


 
    def test_newick_read_tree_string_flags_noprop_underscore(self):
        '''
            Test newick reading in tree from string, with non-propagating underscore flags.
        '''
        t = read_tree(tree = self.string_nopropflags_underscore)
        
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nopropflags, msg = "Couldn't read and parse tree from string properly, non-propagating underscore flags.")


    def test_newick_read_tree_string_flags_noprop_hashtag(self):
        '''
            Test newick reading in tree from string, with non-propagating hashtag flags.
        '''
        t = read_tree(tree = self.string_nopropflags_hashtag)
        
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nopropflags, msg = "Couldn't read and parse tree from string properly, non-propagating hashtag flags.")

 
 


    def test_newick_read_tree_nodenames(self):
        '''
            Test newick reading in tree from string, with node names.
        '''
        t = read_tree(tree = self.string_nodenames)
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nodenames, msg = "Couldn't read and parse tree from string properly when node names provided.")


    def test_newick_read_tree_nodenames_propflags(self):
        '''
            Test newick reading in tree from string, with node names and propagating flags.
        '''
        t = read_tree(tree = self.string_nodenames_propflags)
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nodenames_propflags, msg = "Couldn't read and parse tree from string properly when node names and propagating flags provided.")


    def test_newick_read_tree_nodenames_nopropflags(self):
        '''
            Test newick reading in tree from string, with node names and non-propagating flags.
        '''
        t = read_tree(tree = self.string_nodenames_nopropflags)
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f
        print_tree(t)
        sys.stdout = orig_stdout
        f.close()
        with open("out.txt", "r") as f:
            printed = str(f.read())
        self.assertMultiLineEqual(printed, self.true_nodenames_nopropflags, msg = "Couldn't read and parse tree from string properly when node names and non-propagating flags provided.")
        
        
        
        
        
        
        
        
        
        
        
        
        
