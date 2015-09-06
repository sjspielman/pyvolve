#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

''' 
    Test evolver under a variety of model specifications. Not quite unittests, but definitely tests.
'''

import unittest
import os
from pyvolve import *
ZERO=1e-8
DECIMAL=8




class evolver_mrca(unittest.TestCase):
    ''' 
        Suite of tests for evolver when MRCA is provided.
    '''
    
    def setUp(self):
        ''' 
            Tree and frequency set-up.
        '''
        self.tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" )        

    def test_evolver_mrca_nucleotide(self):
        rootseq = "AAATTTCCCGGG"
        m = Model("nucleotide")
        p = Partition(root_sequence = rootseq, models = m)
        evolve = Evolver(partitions = p, tree = self.tree)
        evolve(ratefile = False, infofile=False, seqfile=False)
        seqdict = evolve.get_sequences(anc = True)
        self.assertTrue(seqdict["root"] == rootseq, msg = "MRCA not preserved for nucleotide evolution.")



    def test_evolver_mrca_aminoacid(self):
        rootseq = "ARGYMMKLPQ"
        m = Model("WAG")
        p = Partition(root_sequence = rootseq, models = m)
        evolve = Evolver(partitions = p, tree = self.tree)
        evolve(ratefile = False, infofile=False, seqfile=False)
        seqdict = evolve.get_sequences(anc = True)
        self.assertTrue(seqdict["root"] == rootseq, msg = "MRCA not preserved for amino acid evolution.")


    def test_evolver_mrca_codon(self):
        rootseq = "AACCGATTTGGCCAT"
        m = Model("codon", {"beta":0.5})
        p = Partition(root_sequence = rootseq, models = m)
        evolve = Evolver(partitions = p, tree = self.tree)
        evolve(ratefile = False, infofile=False, seqfile=False)
        seqdict = evolve.get_sequences(anc = True)
        self.assertTrue(seqdict["root"] == rootseq, msg = "MRCA not preserved for codon evolution.")




class evolver_singlepart_nohet_tests(unittest.TestCase):
    ''' 
        Suite of tests for evolver under temporally homogeneous conditions (no branch heterogeneity!!)
        Single partition.
        All tests are conducted using nucleotides.
    '''
    
    def setUp(self):
        ''' 
            Tree and frequency set-up.
        '''
        self.tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" )
        fobj = EqualFrequencies( 'nucleotide' )
        f = fobj.compute_frequencies()
        
        params = {'state_freqs':f, 'mu':{'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}}
        m1 = Model("nucleotide", params)
        self.part1 = Partition(models = m1, size = 10)




    def test_evolver_singlepart_nohet_get_sequences(self):
        '''
            Test evolver get_sequences method. 
        '''
        
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(ratefile = False, infofile=False, seqfile=False)

        seqdict = evolve.get_sequences()
        self.assertTrue(type(seqdict) is dict, msg = "Evolver .get_sequences(), without ancestors, method does not return a dictionary.")
        self.assertTrue(len(seqdict) == 5, msg = "Evolver .get_sequences(), without ancestors, is wrong size.")
        
        for entry in seqdict:
            self.assertTrue(len(seqdict[entry]) == 10, msg = "Evolver .get_sequences(), without ancestors, returns sequences of the wrong length.")




    def test_evolver_singlepart_nohet_get_sequences_anc(self):
        '''
            Test evolver get_sequences method. 
        '''
        
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(ratefile = False, infofile=False, seqfile=False)

        seqdict = evolve.get_sequences(anc = True)
        self.assertTrue(type(seqdict) is dict, msg = "Evolver .get_sequences(), with ancestors, method does not return a dictionary.")
        self.assertTrue(len(seqdict) == 9, msg = "Evolver .get_sequences(), with ancestors, is wrong size.")
        
        for entry in seqdict:
            self.assertTrue(len(seqdict[entry]) == 10, msg = "Evolver .get_sequences(), with ancestors, returns sequences of the wrong length.")
            
 

    
    def test_evolver_singlepart_nohet_ratefile(self):
        '''
            Test evolver with a single partition, no heterogeneity at all.
            Ensure rate file correct.
        '''
        
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(ratefile = "rates.txt", infofile=False, seqfile=False)
        # Check ratefile
        with open('tests/evolFiles/single_part_no_het_rates.txt', 'r') as ref_h:
            ref = str(ref_h.read())
        with open('rates.txt', 'r') as test_h:
            test = str(test_h.read())
        os.remove("rates.txt")
        self.assertMultiLineEqual(test, ref, msg = "Rate file improperly written for single partition, no het.")



    def test_evolver_singlepart_nohet_seqfile(self):
        '''
            Test evolver with a single partition, no heterogeneity at all.
            Ensure leaf sequences only properly written to seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.fasta", infofile=False, ratefile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        assert(len(aln) == 5), "Wrong number of sequences were written to file when write_anc=False."
        assert(len(aln[0]) == 10), "Output alignment incorrect length."
        
        
        
    def test_evolver_singlepart_nohet_seqfile_anc(self):
        '''
            Test evolver with a single partition, no heterogeneity at all.
            Ensure ancestors properly written to seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.fasta", write_anc = True, infofile=False, ratefile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        assert(len(aln) == 9), "Wrong number of sequences were written to file when write_anc=False."
        assert(len(aln[0]) == 10), "Output alignment incorrect length."
        


    def test_evolver_singlepart_nohet_seqfile_phy(self):
        '''
            Test evolver with a single partition, no heterogeneity at all.
            Ensure can save in phylip for seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.phy", seqfmt = "phylip", infofile=False, ratefile=False)
        
        # Check seqfile, no ancestors
        try:
            aln = AlignIO.read("out.phy", "phylip")
            os.remove("out.phy")
        except:
            raise AssertionError("Couldn't overwrite FASTA default to save as PHYLIP.")
            os.remove("out.phy")


class evolver_twopart_nohet_tests(unittest.TestCase):
    ''' 
        Suite of tests for evolver under temporally homogeneous conditions (no branch heterogeneity!!)
        Two partitions.
        All tests are conducted using nucleotides.
    '''
    
    def setUp(self):
        ''' 
            Tree and frequency set-up.
        '''
        self.tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" )
        fobj = EqualFrequencies( 'nucleotide' )
        f = fobj.compute_frequencies()
        params = {'state_freqs':f, 'mu':{'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}}
       
        m1 = Model('nucleotide', params)
        self.part1 = Partition(models = m1, size = 10)
        
        m2 = Model('nucleotide', params)
        self.part2 = Partition(models = m2, size = 12)

  
        
    def test_evolver_twopart_nohet_ratefile(self):
        '''
            Test evolver with two partitions, no heterogeneity at all.
            Ensure rate file correct.
        '''
        
        evolve = Evolver(partitions = [self.part1, self.part2], tree = self.tree)
        evolve(ratefile = "rates.txt", seqfile=False, infofile=False)
        # Check ratefile
        with open('tests/evolFiles/two_part_no_het_rates.txt', 'r') as ref_h:
            ref = str(ref_h.read())
        with open('rates.txt', 'r') as test_h:
            test = str(test_h.read())
        os.remove("rates.txt")
        self.assertMultiLineEqual(test, ref, msg = "Rate file improperly written for two partitions, no het.")



    def test_evolver_twopart_nohet_seqfile(self):
        '''
            Test evolver with two partitions, no heterogeneity at all.
            Ensure leaf sequences only properly written to seqfile.
        '''
        evolve = Evolver(partitions = [self.part1, self.part2], tree = self.tree)
        evolve(seqfile = "out.fasta", ratefile=False, infofile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        self.assertTrue(len(aln) == 5, msg="Wrong number of sequences were written to file when write_anc=False.")
        self.assertTrue(len(aln[0]) == 22, msg="Output alignment incorrect length.")
        

    def test_evolver_twopart_nohet_seqfile_anc(self):
        '''
            Test evolver with two partitions, no heterogeneity at all.
            Ensure ancestors properly written to seqfile.
        '''
        evolve = Evolver(partitions = [self.part1, self.part2], tree = self.tree)
        evolve(seqfile = "out.fasta", write_anc = True, ratefile=False, infofile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        self.assertTrue(len(aln) == 9, msg="Wrong number of sequences were written to file when write_anc=True.")
        self.assertTrue(len(aln[0]) == 22, msg="Output alignment incorrect length.")
        
        
        




class evolver_sitehet_tests(unittest.TestCase):
    ''' 
        Suite of tests for evolver under temporally homogeneous conditions (no branch heterogeneity!!).
        Uses a single partition with nucleotides WITH site heterogeneity.
    '''
    
    def setUp(self):
        ''' 
            Tree and frequency set-up.
        '''
        self.tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" )
        fobj = EqualFrequencies( 'nucleotide' )
        f = fobj.compute_frequencies()
        
        
        params = {'state_freqs':f, 'mu':{'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}}
        m1 = Model('nucleotide', params, rate_factors = [2.0783848 ,  0.89073634,  0.05938242], rate_probs = [0.33, 0.33, 0.34])
        self.part1 = Partition(models = m1, size = 12)

                
    def test_evolver_sitehet_ratefile(self):
        '''
            Test evolver with one partition, site heterogeneity.
            Ensure rate file correct.
        '''
        
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(ratefile = "rates.txt", seqfile = False, infofile = False)
        # Check ratefile
        test = []
        with open('rates.txt', 'r') as test_h:
            for line in test_h:
                test.append(line)
        os.remove("rates.txt")
        self.assertTrue( len(test) == 13 , msg="Ratefile improperly written for single partition, site het (wrong num lines).")
        for i in range(1, 13):
            self.assertRegexpMatches( test[i], str(i) + "\t1\t[123]", msg = "Ratefile improperly written for single partition, site het (wrong line contents). NOTE: THIS FAILURE MAY HAVE RESULTED DUE TO PYTHON 2vs3 INCOMPATIBILITY. Please contact the author.")


    def test_evolver_sitehet_infofile(self):
        '''
            Test evolver with one partition, site heterogeneity.
            Ensure rate info file correct.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(infofile = "info.txt", seqfile = False, ratefile = False)
        test = []
        with open('info.txt', 'r') as info_h:
            for line in info_h:
                test.append(line)
        os.remove("info.txt")
        self.assertTrue( len(test) == 4, msg="Infofile improperly written for single partition, site het (wrong num lines).")
        for i in range(1, 4):
            self.assertRegexpMatches( test[i],  "1\tNone\t" + str(i) + "\t" + str(round(self.part1.models[0].rate_probs[i-1],4)) + "\t" + str(round(self.part1.models[0].rate_factors[i-1],4)), msg = "Infofile improperly written for single partition, site het (wrong line contents). NOTE: THIS FAILURE MAY HAVE RESULTED DUE TO PYTHON 2vs3 INCOMPATIBILITY. Please contact the author.")

        
        
        
        

    def test_evolver_sitehet_seqfile(self):
        '''
            Test evolver with one partition, site heterogeneity.
            Ensure leaf sequences only properly written to seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.fasta", infofile = False, ratefile = False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        self.assertTrue(len(aln) == 5, msg="Wrong number of sequences were written to file when write_anc=False.")
        self.assertTrue(len(aln[0]) == 12, msg="Output alignment incorrect length.")
        

    def test_evolver_sitefile_seqfile_anc(self):
        '''
            Test evolver with one partition, site heterogeneity.
            Ensure ancestors properly written to seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.fasta", write_anc = True, ratefile=False, infofile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        self.assertTrue(len(aln) == 9, msg="Wrong number of sequences were written to file when write_anc=True.")
        self.assertTrue(len(aln[0]) == 12, msg="Output alignment incorrect length.")







class evolver_branchhet_tests(unittest.TestCase):
    ''' 
        Suite of tests for evolver under branch heterogeneity.
        Uses a single partition with nucleotides *without* site heterogeneity.
    '''
    
    def setUp(self):
        ''' 
            Tree and frequency set-up.
        '''
        self.tree = read_tree( tree = "(((t2:0.36_m2_,t1:0.45):0.001,t3:0.77):0.44_m1_,(t5:0.77,t4:0.41):0.89);" )
        fobj = EqualFrequencies( 'nucleotide' )
        f = fobj.compute_frequencies()
        params = {'state_freqs':f, 'mu':{'AC':1, 'AG':1, 'AT':1, 'CG':1, 'CT':1, 'GT':1}}
        type = 'nucleotide'
        
        root = Model(type, params)
        root.assign_name('root_model')
   
        m1 = Model(type, params)
        m1.assign_name('m1')

        m2 = Model(type, params)
        m2.assign_name('m2')
           
        
        self.part1 = Partition(models = [ m1, m2, root ], size = 10, root_model_name = "root_model")

        
        
    def test_evolver_branchhet_ratefile(self):
        '''
            Test evolver with one partition, branch heterogeneity.
            Ensure rate file correct.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(ratefile = "rates.txt", seqfile=False, infofile=False)
        # Check ratefile
        with open('tests/evolFiles/single_part_no_het_rates.txt', 'r') as ref_h:
            ref = str(ref_h.read())
        with open('rates.txt', 'r') as test_h:
            test = str(test_h.read())
        os.remove("rates.txt")
        self.assertMultiLineEqual(test, ref, msg = "Rate file improperly written for single partition, branch het.")


    def test_evolver_sitehet_seqfile(self):
        '''
            Test evolver with one partition, branch heterogeneity.
            Ensure leaf sequences only properly written to seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.fasta", ratefile=False, infofile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        self.assertTrue(len(aln) == 5, msg="Wrong number of sequences were written to file when write_anc=False, branch het.")
        self.assertTrue(len(aln[0]) == 10, msg="Output alignment incorrect length, branch het.")
        

    def test_evolver_sitefile_seqfile_anc(self):
        '''
            Test evolver with one partition, branch heterogeneity.
            Ensure ancestors properly written to seqfile.
        '''
        evolve = Evolver(partitions = self.part1, tree = self.tree)
        evolve(seqfile = "out.fasta", write_anc = True, ratefile=False, infofile=False)
        
        # Check seqfile, no ancestors
        aln = AlignIO.read("out.fasta", "fasta")
        os.remove("out.fasta")
        self.assertTrue(len(aln) == 9, msg="Wrong number of sequences were written to file when write_anc=True, branch het.")
        self.assertTrue(len(aln[0]) == 10, msg= "Output alignment incorrect length, branch het.")
 
 
 
        



class evolver_setcode(unittest.TestCase):
    '''
        Tests to ensure that the self._code is properly setup.
    '''
   
    def setUp(self):
        ''' 
            genetics set-up.
        '''
        self.g = Genetics()
    
    



    def test_evolver_setcode_codon(self):
        '''
            Codons properly assigned as code in codon model?"
        '''
        tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" ) 
        m = Model("GY", {"omega":1.5})
        p = Partition(models = m, size = 5)
        evolve = Evolver(partitions = p, tree = tree)
        self.assertTrue(evolve._code == self.g.codons, msg = "Codons not properly assigned as code for a codon model.")


    def test_evolver_setcode_mutsel_codon(self):
        '''
            Codons properly assigned as code in mutsel model?"
        '''
        tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" ) 
        m = Model("mutsel", {"state_freqs":np.repeat(1./61, 61)})
        p = Partition(models = m, size = 5)
        evolve = Evolver(partitions = p, tree = tree)
        self.assertTrue(evolve._code == self.g.codons, msg = "Codons not properly assigned as code for a mutsel codon model.")

            
    def test_evolver_setcode_mutsel_nucleotide(self):
        '''
            Nucleotides properly assigned as code in mutsel model?"
        '''
        tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" ) 
        m = Model("mutsel", {"state_freqs":np.repeat(1./4, 4)})
        p = Partition(models = m, size = 5)
        evolve = Evolver(partitions = p, tree = tree)
        self.assertTrue(evolve._code == self.g.nucleotides, msg = "Nucleotides not properly assigned as code for a mutsel nuc model.")
  

    def test_evolver_setcode_custom(self):
        '''
            Custom code properly assigned as code for custom model?"
        '''
        tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" ) 
        matrix = np.array([ [-0.5, 0.25, 0.25], [0.25, -0.5, 0.25], [0.25, 0.25, -0.5] ]) 
        code = ["0", "1", "2"]
        m = Model("custom", {"matrix":matrix, "code":code})
        p = Partition(models = m, size = 5)
        evolve = Evolver(partitions = p, tree = tree)
        self.assertTrue(evolve._code == code, msg = "Custom not properly assigned as code for a custom model.")

    def test_evolver_setcode_nucleotide(self):
        '''
            Nucleotides properly assigned as code in nuc model?"
        '''
        tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" ) 
        m1 = Model("nucleotide")
        p = Partition(models = m1, size = 5)
        evolve = Evolver(partitions = p, tree = tree)
        self.assertTrue(evolve._code == self.g.nucleotides, msg = "Nucleotides not properly assigned as code for a nucleotide model.")

    def test_evolver_setcode_aminoacid(self):
        '''
            Amino acids properly assigned as code in aa model?"
        '''
        tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" ) 
        m = Model("JTT")
        p = Partition(models = m, size = 5)
        evolve = Evolver(partitions = p, tree = tree)
        self.assertTrue(evolve._code == self.g.amino_acids, msg = "Amino acids not properly assigned as code for an AA model.")


            
# def run_evolver_test():
# 
#     run_tests = unittest.TextTestRunner()
# 
#     print "Testing evolver no het, one partition"
#     test_suite0 = unittest.TestLoader().loadTestsFromTestCase(evolver_singlepart_nohet_tests)
#     run_tests.run(test_suite0)
# 
#     print "Testing evolver no het, two partitions"
#     test_suite1 = unittest.TestLoader().loadTestsFromTestCase(evolver_twopart_nohet_tests)
#     run_tests.run(test_suite1)
# 
#     print "Testing evolver site het, one partition"
#     test_suite2 = unittest.TestLoader().loadTestsFromTestCase(evolver_sitehet_tests)
#     run_tests.run(test_suite2)
#     
#     print "Testing evolver branch het, one partition"
#     test_suite3 = unittest.TestLoader().loadTestsFromTestCase(evolver_branchhet_tests)
#     run_tests.run(test_suite3)
            
      
            
            
            
            
            
            
            
            