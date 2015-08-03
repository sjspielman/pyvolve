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
 
 
 
        
          
          
class evolver_noisy_branch_lengths_tests(unittest.TestCase):
    ''' 
        Suite of tests for evolver to ensure that matrices are being properly up and that branch lengths are begin drawn ok.
    '''     
    
    def setUp(self):
        ''' 
            Tree and frequency set-up.
        '''
        
        self.tree = read_tree( tree = "(((t2:0.36,t1:0.45):0.001,t3:0.77):0.44,(t5:0.77,t4:0.41):0.89);" )
        fobj = EqualFrequencies( 'codon' )
        f = fobj.compute_frequencies()
        
        params = {'state_freqs':f, 'kappa':3.5, 'omega':1.5}
        self.m1 = Model("codon", params)

    
        self.Q =  np.array([[-1.349975704555383960e+00, 6.364679535939857247e-02, 1.277097380283088168e-01, 7.919074340016568625e-02, 1.439089117091822689e-01, 0.0, 0.0, 0.0, 3.227081589895426372e-01, 0.0, 0.0, 0.0, 1.471601513268189498e-01, 0.0, 0.0, 0.0, 9.823916713055715066e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.674120386114100301e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.190800144017270801e-01, -8.111253192454067618e-01, 7.662584281698529842e-02, 1.319845723336094678e-01, 0.0, 1.199587104099346102e-01, 0.0, 0.0, 0.0, 9.050053333775449904e-02, 0.0, 0.0, 0.0, 2.743545360204182118e-02, 0.0, 0.0, 0.0, 1.168601017315806651e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.363908796195048845e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.504100264982280721e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.984666906695451427e-01, 6.364679535939857247e-02, -9.238071720597268444e-01, 7.919074340016568625e-02, 0.0, 0.0, 9.037343479593649975e-02, 0.0, 0.0, 0.0, 1.173579081806063878e-01, 0.0, 0.0, 0.0, 2.466472050901183377e-02, 0.0, 0.0, 0.0, 1.726023602629593867e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.775045188821034459e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.190800144017270801e-01, 1.060779922656642782e-01, 7.662584281698529842e-02, -9.879788158155174971e-01, 0.0, 0.0, 0.0, 6.767091186769211286e-02, 0.0, 0.0, 0.0, 2.947906887808233178e-01, 0.0, 0.0, 0.0, 2.876245830508082982e-02, 0.0, 0.0, 0.0, 1.553419280862226581e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.006183285564997698e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.956714643567183720e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.190800144017270801e-01, 0.0, 0.0, 0.0, -1.211134648584501372e+00, 7.997247360662307347e-02, 1.506223913265608561e-01, 4.511394124512807524e-02, 1.290832635958170660e-01, 0.0, 0.0, 0.0, 3.679003783170473607e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.221210415561683159e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.106262563608926330e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.652382557508824845e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 6.364679535939857247e-02, 0.0, 0.0, 9.593927447278817466e-02, -9.824150157043188836e-01, 6.024895653062433548e-02, 1.127848531128201881e-01, 0.0, 3.620021333510180239e-02, 0.0, 0.0, 0.0, 6.858863400510455122e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.048405264624480276e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.504331893387535035e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.973257308727992931e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 7.662584281698529842e-02, 0.0, 2.398481861819704297e-01, 7.997247360662307347e-02, -9.810530313793327517e-01, 4.511394124512807524e-02, 0.0, 0.0, 4.694316327224255792e-02, 0.0, 0.0, 0.0, 6.166180127252958443e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.117061183535621061e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.464910238155726963e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.432259873329247979e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 7.919074340016568625e-02, 9.593927447278817466e-02, 1.999311840165576559e-01, 6.024895653062433548e-02, -8.439058287624919830e-01, 0.0, 0.0, 0.0, 1.179162755123293188e-01, 0.0, 0.0, 0.0, 7.190614576270207281e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.370515304155900116e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.103119265075475863e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.069052600098004860e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [2.977000360043177141e-01, 0.0, 0.0, 0.0, 1.439089117091822689e-01, 0.0, 0.0, 0.0, -9.967407148622288160e-01, 3.620021333510180239e-02, 7.823860545373759190e-02, 1.179162755123293188e-01, 1.471601513268189498e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.423620178618960619e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.113803197345517859e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 1.591169883984964173e-01, 0.0, 0.0, 0.0, 1.199587104099346102e-01, 0.0, 0.0, 1.290832635958170660e-01, -1.170400259162157619e+00, 4.694316327224255792e-02, 1.965271258538822119e-01, 0.0, 2.743545360204182118e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.049094320248859930e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.365443026494024326e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.498818193554544986e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.915646070424632530e-01, 0.0, 0.0, 0.0, 9.037343479593649975e-02, 0.0, 2.151387726596950822e-01, 3.620021333510180239e-02, -9.855502720423577889e-01, 1.179162755123293188e-01, 0.0, 0.0, 2.466472050901183377e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.662126838542537843e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.220946797366974423e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.097630006569721983e-02, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.979768585004142156e-01, 0.0, 0.0, 0.0, 6.767091186769211286e-02, 1.290832635958170660e-01, 6.033368889183633732e-02, 4.694316327224255792e-02, -7.856847210063072628e-01, 0.0, 0.0, 0.0, 2.876245830508082982e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.215271444449976701e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.349724737527979002e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.926441475344454812e-02, 0.0, 0.0, 0.0, 0.0], [1.190800144017270801e-01, 0.0, 0.0, 0.0, 3.597722792729556307e-01, 0.0, 0.0, 0.0, 1.290832635958170660e-01, 0.0, 0.0, 0.0, -1.296960509868487765e+00, 1.829030240136121643e-02, 6.166180127252958443e-02, 1.917497220338722219e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.207358666313057827e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.075158091099481128e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.164620097945595206e-02, 0.0, 0.0, 0.0], [0.0, 6.364679535939857247e-02, 0.0, 0.0, 0.0, 2.998967760248365533e-01, 0.0, 0.0, 0.0, 3.620021333510180239e-02, 0.0, 0.0, 9.810676755121261472e-02, -8.608616761463714351e-01, 2.466472050901183377e-02, 4.793743050846805548e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.671139981280925224e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.221088784918274023e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.186096238422058337e-03, 0.0, 0.0], [0.0, 0.0, 7.662584281698529842e-02, 0.0, 0.0, 0.0, 2.259335869898412563e-01, 0.0, 0.0, 0.0, 4.694316327224255792e-02, 0.0, 3.679003783170473607e-01, 2.743545360204182118e-02, -1.161178816982928019e+00, 2.876245830508082982e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.215878636189055673e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.540769091681253866e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.513422381496730695e-01, 0.0], [0.0, 0.0, 0.0, 7.919074340016568625e-02, 0.0, 0.0, 0.0, 1.691772796692302960e-01, 0.0, 0.0, 0.0, 1.179162755123293188e-01, 9.810676755121261472e-02, 4.572575600340304108e-02, 2.466472050901183377e-02, -9.345291863553736311e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.833975078917709250e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.483425822025813190e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.030653107182625539e-01], [1.190800144017270801e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.199751187587493551e+00, 1.168601017315806651e-01, 2.876706004382656445e-01, 1.553419280862226581e-01, 1.221210415561683159e-02, 0.0, 0.0, 0.0, 2.408857566982110232e-01, 0.0, 0.0, 0.0, 1.207358666313057827e-01, 0.0, 0.0, 0.0, 1.469648154445640231e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 6.364679535939857247e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.823916713055715066e-02, -1.315677782692045383e+00, 1.726023602629593867e-01, 2.589032134770377080e-01, 0.0, 1.048405264624480276e-01, 0.0, 0.0, 0.0, 2.622735800622149616e-01, 0.0, 0.0, 0.0, 1.671139981280925224e-01, 0.0, 0.0, 0.0, 2.545563518478019885e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.626025066245570250e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 7.662584281698529842e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.637319452175952650e-01, 1.168601017315806651e-01, -9.942207800478171098e-01, 1.553419280862226581e-01, 0.0, 0.0, 4.117061183535621061e-02, 0.0, 0.0, 0.0, 2.873297564453451969e-01, 0.0, 0.0, 0.0, 8.215878636189055673e-02, 0.0, 0.0, 0.0, 7.100180755284137002e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 7.919074340016568625e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.823916713055715066e-02, 1.947668362193011315e-01, 1.726023602629593867e-01, -1.183514773560439215e+00, 0.0, 0.0, 0.0, 1.370515304155900116e-01, 0.0, 0.0, 0.0, 2.303817861112494036e-01, 0.0, 0.0, 0.0, 4.833975078917709250e-02, 0.0, 0.0, 0.0, 2.402473314225999149e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.989178660891795791e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 1.439089117091822689e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.823916713055715066e-02, 0.0, 0.0, 0.0, -1.195781172223651545e+00, 6.989368430829868972e-02, 6.861768639226034638e-02, 9.136768694372666977e-02, 9.635430267928440928e-02, 0.0, 0.0, 0.0, 3.018396665782644428e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.425050254435705044e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.413095639377206558e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 1.199587104099346102e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.168601017315806651e-01, 0.0, 0.0, 8.141402770411221637e-03, -1.388025642626965928e+00, 2.744707455690414272e-02, 2.284192173593166952e-01, 0.0, 1.049094320248859930e-01, 0.0, 0.0, 0.0, 4.177849953202313338e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.401732757355013959e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.243314327181998302e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.037343479593649975e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.726023602629593867e-01, 0.0, 2.035350692602804976e-02, 6.989368430829868972e-02, -1.221580919578354685e+00, 9.136768694372666977e-02, 0.0, 0.0, 1.149319025781380676e-01, 0.0, 0.0, 0.0, 2.053969659047263918e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.859640952622908405e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.580649683323119947e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.767091186769211286e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.553419280862226581e-01, 8.141402770411221637e-03, 1.747342107707466896e-01, 2.744707455690414272e-02, -7.434764115321713662e-01, 0.0, 0.0, 0.0, 9.215271444449976701e-02, 0.0, 0.0, 0.0, 1.208493769729427486e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.041247706030190415e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.672631500245011282e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.605550906387803012e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.455979178263928697e-01, 0.0, 0.0, 0.0, 1.221210415561683159e-02, 0.0, 0.0, 0.0, -8.320814608475014529e-01, 6.993962134992399993e-02, 1.915531709635634461e-01, 6.143514296299984467e-02, 1.207358666313057827e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.455212789382071575e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.620021333510180239e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.921502543289516973e-01, 0.0, 0.0, 0.0, 1.048405264624480276e-01, 0.0, 0.0, 6.423620178618960619e-02, -1.364072589282105730e+00, 7.662126838542537843e-02, 1.535878574074996117e-01, 0.0, 1.671139981280925224e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.461772105976097580e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.747045483886362605e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.129544218149503398e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.315059006573984668e-01, 0.0, 0.0, 0.0, 4.117061183535621061e-02, 0.0, 1.605905044654740155e-01, 6.993962134992399993e-02, -9.943746318734599798e-01, 6.143514296299984467e-02, 0.0, 0.0, 8.215878636189055673e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.883787189467898249e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.744075016424304611e-02, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.179162755123293188e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.883548202155566176e-01, 0.0, 0.0, 0.0, 1.370515304155900116e-01, 6.423620178618960619e-02, 1.748490533748099651e-01, 7.662126838542537843e-02, -1.238928836312801263e+00, 0.0, 0.0, 0.0, 4.833975078917709250e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.339889895011191601e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.981610368836113634e-01, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.471601513268189498e-01, 0.0, 0.0, 0.0, 9.823916713055715066e-02, 0.0, 0.0, 0.0, 3.053026038904207984e-02, 0.0, 0.0, 0.0, 9.635430267928440928e-02, 0.0, 0.0, 0.0, -9.186010166834390755e-01, 1.114093320853950381e-01, 1.369313106031509464e-01, 3.222650052611806398e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.630063236439792396e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.027436682990932465e-01, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.743545360204182118e-02, 0.0, 0.0, 0.0, 1.168601017315806651e-01, 0.0, 0.0, 0.0, 2.621013161561200899e-01, 0.0, 0.0, 0.0, 1.049094320248859930e-01, 0.0, 0.0, 8.049057775420385974e-02, -7.789444488181741511e-01, 5.477252424126036884e-02, 8.056625131529515649e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.884355139673095952e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.965240596055145842e-03, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.466472050901183377e-02, 0.0, 0.0, 0.0, 1.726023602629593867e-01, 0.0, 0.0, 0.0, 1.029265295883905196e-01, 0.0, 0.0, 0.0, 1.149319025781380676e-01, 0.0, 2.012264443855096563e-01, 1.114093320853950381e-01, -1.073855617185561329e+00, 3.222650052611806398e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.163076366725015603e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.522370635827885232e-01, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.876245830508082982e-02, 0.0, 0.0, 0.0, 1.553419280862226581e-01, 0.0, 0.0, 0.0, 3.426288260389749873e-01, 0.0, 0.0, 0.0, 9.215271444449976701e-02, 8.049057775420385974e-02, 2.785233302134875744e-01, 5.477252424126036884e-02, -1.389672668760419061e+00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.933703288103253037e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.576632767956563708e-01], [2.977000360043177141e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.823916713055715066e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.223930633062055184e-01, 2.545563518478019885e-02, 1.183363459214022972e-01, 2.402473314225999149e-02, 8.425050254435705044e-02, 0.0, 0.0, 0.0, 1.113803197345517859e-01, 0.0, 0.0, 0.0, 1.630063236439792396e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 1.591169883984964173e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.168601017315806651e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.469648154445640231e-01, -1.024587067462706891e+00, 7.100180755284137002e-02, 4.004122190376665363e-02, 0.0, 1.401732757355013959e-01, 0.0, 0.0, 0.0, 2.365443026494024326e-01, 0.0, 0.0, 0.0, 4.884355139673095952e-02, 0.0, 0.0, 6.504100264982280721e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.915646070424632530e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.726023602629593867e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.449413590742733349e-01, 2.545563518478019885e-02, -1.040910547636912931e+00, 2.402473314225999149e-02, 0.0, 0.0, 9.859640952622908405e-02, 0.0, 0.0, 0.0, 2.220946797366974423e-01, 0.0, 0.0, 0.0, 6.163076366725015603e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.979768585004142156e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.553419280862226581e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.469648154445640231e-01, 4.242605864130033488e-02, 7.100180755284137002e-02, -8.965253719776286045e-01, 0.0, 0.0, 0.0, 2.041247706030190415e-02, 0.0, 0.0, 0.0, 8.349724737527979002e-02, 0.0, 0.0, 0.0, 9.933703288103253037e-02, 0.0, 7.956714643567183720e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 3.597722792729556307e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.221210415561683159e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.469648154445640231e-01, 0.0, 0.0, 0.0, -1.338925479192911050e+00, 9.344885049033425928e-02, 1.643273492103818068e-01, 1.360831804020126885e-02, 4.455212789382071575e-02, 0.0, 0.0, 0.0, 4.075158091099481128e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 9.652382557508824845e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 2.998967760248365533e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.048405264624480276e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.545563518478019885e-02, 0.0, 0.0, 5.616700169623803363e-02, -8.925708467918269662e-01, 6.573093968415272270e-02, 3.402079510050316780e-02, 0.0, 9.461772105976097580e-02, 0.0, 0.0, 0.0, 1.221088784918274023e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 8.973257308727992931e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.259335869898412563e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.117061183535621061e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.100180755284137002e-02, 0.0, 1.404175042405950979e-01, 9.344885049033425928e-02, -9.717214475448986422e-01, 1.360831804020126885e-02, 0.0, 0.0, 8.883787189467898249e-02, 0.0, 0.0, 0.0, 1.540769091681253866e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 1.432259873329247979e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.691772796692302960e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.370515304155900116e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.402473314225999149e-02, 5.616700169623803363e-02, 2.336221262258356968e-01, 6.573093968415272270e-02, -9.982056179869800427e-01, 0.0, 0.0, 0.0, 3.339889895011191601e-02, 0.0, 0.0, 0.0, 2.483425822025813190e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 3.069052600098004860e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.227081589895426372e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.635430267928440928e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.674120386114100301e-01, 0.0, 0.0, 0.0, 8.425050254435705044e-02, 0.0, 0.0, 0.0, -1.267138859632953807e+00, 6.307848070650731720e-02, 1.480631198244649893e-01, 2.226593263340794285e-02, 1.630063236439792396e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.050053333775449904e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.049094320248859930e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.363908796195048845e-02, 0.0, 0.0, 0.0, 1.401732757355013959e-01, 0.0, 0.0, 2.970141859588047370e-02, -7.425391979214640559e-01, 5.922524792978599295e-02, 5.566483158351986232e-02, 0.0, 4.884355139673095952e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.498818193554544986e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.173579081806063878e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.149319025781380676e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.775045188821034459e-01, 0.0, 0.0, 0.0, 9.859640952622908405e-02, 0.0, 7.425354648970118598e-02, 6.307848070650731720e-02, -7.405957627296407830e-01, 2.226593263340794285e-02, 0.0, 0.0, 6.163076366725015603e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.097630006569721983e-02, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.947906887808233178e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.215271444449976701e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.006183285564997698e-02, 0.0, 0.0, 0.0, 2.041247706030190415e-02, 2.970141859588047370e-02, 1.576962017662683069e-01, 5.922524792978599295e-02, -8.926420290676866376e-01, 0.0, 0.0, 0.0, 9.933703288103253037e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.926441475344454812e-02, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.679003783170473607e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.207358666313057827e-01, 0.0, 0.0, 0.0, 1.469648154445640231e-01, 0.0, 0.0, 0.0, 2.106262563608926330e-01, 0.0, 0.0, 0.0, 4.455212789382071575e-02, 0.0, 0.0, 0.0, -1.153930641257678857e+00, 3.256236759782064200e-02, 1.027179394454169198e-01, 6.622468858735502950e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.164620097945595206e-02, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.858863400510455122e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.671139981280925224e-01, 0.0, 0.0, 0.0, 2.545563518478019885e-02, 0.0, 0.0, 0.0, 3.504331893387535035e-01, 0.0, 0.0, 0.0, 9.461772105976097580e-02, 0.0, 0.0, 1.086708824293194930e-01, -1.022715053630787851e+00, 4.108717577816677069e-02, 1.655617214683875738e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.186096238422058337e-03, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.166180127252958443e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.215878636189055673e-02, 0.0, 0.0, 0.0, 7.100180755284137002e-02, 0.0, 0.0, 0.0, 2.464910238155726963e-01, 0.0, 0.0, 0.0, 8.883787189467898249e-02, 0.0, 2.716772060732987604e-01, 3.256236759782064200e-02, -1.071957791305660601e+00, 6.622468858735502950e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.513422381496730695e-01, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.190614576270207281e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.833975078917709250e-02, 0.0, 0.0, 0.0, 2.402473314225999149e-02, 0.0, 0.0, 0.0, 5.103119265075475863e-02, 0.0, 0.0, 0.0, 3.339889895011191601e-02, 1.086708824293194930e-01, 8.140591899455161540e-02, 4.108717577816677069e-02, -5.629300092153062263e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.030653107182625539e-01], [0.0, 6.364679535939857247e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.921502543289516973e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.545563518478019885e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.794878133135884735e-01, 1.326119107261197194e-01, 0.0, 8.973257308727992931e-02, 0.0, 0.0, 3.747045483886362605e-01, 0.0, 0.0, 0.0, 1.186096238422058337e-03, 0.0, 0.0], [0.0, 0.0, 0.0, 7.919074340016568625e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.883548202155566176e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.402473314225999149e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.084016710830380259e-01, -9.318888414438744050e-01, 0.0, 0.0, 0.0, 3.069052600098004860e-02, 0.0, 0.0, 1.981610368836113634e-01, 0.0, 0.0, 0.0, 1.030653107182625539e-01], [0.0, 0.0, 0.0, 0.0, 1.439089117091822689e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.053026038904207984e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.425050254435705044e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.317972220382692461e-01, 5.982171539151994594e-02, 2.387099788882079965e-01, 2.046035066732003124e-02, 0.0, 0.0, 0.0, 1.541155024486398628e-01, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 1.199587104099346102e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.621013161561200899e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.401732757355013959e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.504100264982280721e-02, 0.0, 6.434921705005883230e-02, -9.511054501765308089e-01, 9.548399155528319859e-02, 5.115087666830007984e-02, 1.498818193554544986e-01, 0.0, 0.0, 0.0, 2.965240596055145842e-03, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.037343479593649975e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.029265295883905196e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.859640952622908405e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.608730426251470669e-01, 5.982171539151994594e-02, -9.223833780344231625e-01, 2.046035066732003124e-02, 0.0, 1.097630006569721983e-02, 0.0, 0.0, 0.0, 3.783555953741826738e-01, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.767091186769211286e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.426288260389749873e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.041247706030190415e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.956714643567183720e-02, 6.434921705005883230e-02, 1.495542884787998683e-01, 9.548399155528319859e-02, -1.156594550035883628e+00, 0.0, 0.0, 7.926441475344454812e-02, 0.0, 0.0, 0.0, 2.576632767956563708e-01], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.620021333510180239e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.622735800622149616e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.461772105976097580e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.626025066245570250e-01, 0.0, 0.0, 8.973257308727992931e-02, 0.0, 0.0, -7.896963483954415608e-01, 1.097630006569721983e-02, 1.321073579224075756e-01, 0.0, 1.186096238422058337e-03, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.694316327224255792e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.873297564453451969e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.883787189467898249e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.432259873329247979e-01, 0.0, 1.498818193554544986e-01, -9.468252512037638180e-01, 7.926441475344454812e-02, 0.0, 0.0, 1.513422381496730695e-01, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.179162755123293188e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.303817861112494036e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.339889895011191601e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 1.989178660891795791e-01, 0.0, 0.0, 0.0, 3.069052600098004860e-02, 2.498030322590908681e-01, 1.097630006569721983e-02, -9.751499957069009739e-01, 0.0, 0.0, 0.0, 1.030653107182625539e-01], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.471601513268189498e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.012264443855096563e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.630063236439792396e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 2.413095639377206558e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.109190953833501636e+00, 1.186096238422058337e-03, 2.522370635827885232e-01, 1.030653107182625539e-01], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.743545360204182118e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.177849953202313338e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.884355139673095952e-02, 0.0, 0.0, 6.504100264982280721e-02, 0.0, 0.0, 2.243314327181998302e-01, 0.0, 0.0, 1.498818193554544986e-01, 0.0, 0.0, 6.164620097945595206e-02, -1.318082212035381229e+00, 1.513422381496730695e-01, 1.717755178637709323e-01], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.466472050901183377e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.369313106031509464e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.163076366725015603e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 3.580649683323119947e-01, 0.0, 0.0, 1.097630006569721983e-02, 0.0, 1.027436682990932465e-01, 1.186096238422058337e-03, -7.992631384332000710e-01, 1.030653107182625539e-01], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.876245830508082982e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.208493769729427486e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.933703288103253037e-02, 0.0, 7.956714643567183720e-02, 0.0, 0.0, 0.0, 7.672631500245011282e-02, 0.0, 0.0, 7.926441475344454812e-02, 6.164620097945595206e-02, 1.976827064036764184e-03, 1.513422381496730695e-01, -6.994720105437883806e-01]])
        self.t = 0.05

    def tearDown(self):
        ''' 
            The tearDown deletes the output evolver files.
            
        '''
        try:
            os.remove("simulated_alignment.fasta")   
        except:
            pass
        try:
            os.remove("site_rates.txt")   
        except:
            pass
        try:
            os.remove("site_rates_info.txt")   
        except:
            pass


    def test_evolver_noisy_branch_lengths_false(self):
        
        part1 = Partition(models = self.m1, size = 50)
        evolve = Evolver(partitions = part1, tree = self.tree)
        evolve(seqfile = False, ratefile=False, infofile=False)
        mat, map = evolve._generate_transition_matrices(self.Q, self.t)
        

        self.assertTrue(len(mat) ==  1, msg = "Too many matrices made without noisy branch lengths. There should just be 1!")
        self.assertTrue(len(map) == 50, msg = "Incorrect map length without noisy branch lengths.")
        np.testing.assert_array_almost_equal(map, np.zeros(50), decimal = DECIMAL, err_msg = "Mapping does not all map to 0 without noisy branch lengths.")



    def test_evolver_noisy_branch_lengths_matrix_mapping(self):
        '''
            Test matrix, mapping using normal branch lengths with default number of categories.
        '''
            
        
        part1 = Partition(models = self.m1, size = 50)
        evolve = Evolver(partitions = part1, tree = self.tree, branch_lengths = {"dist":"normal", "sd":0.1})
        evolve(seqfile = False, ratefile=False, infofile=False)
        mat, map = evolve._generate_transition_matrices(self.Q, self.t)
        
        for entry in map:
            try:
                dummy = mat[map]
            except IndexError:
                # yeah, yeah, yeah.
                self.assertTrue(1 == 7, msg = "Mapping values do not have corresponding matrices, for default noisy branch lengths.")



    def test_evolver_noisy_branch_lengths_num_categories(self):
        '''
            Test num_categories correctly setup using normal branch lengths with default number of categories.
        '''
            
        
        part1 = Partition(models = self.m1, size = 50)
        evolve = Evolver(partitions = part1, tree = self.tree, branch_lengths = {"dist":"normal", "sd":0.1})
        self.assertTrue(evolve.bl_noise["num_categories"] ==  5, msg = "Incorrect number of bl categories when default 10%.")


    def test_evolver_noisy_branch_lengths_fulln(self):
        
        part1 = Partition(models = self.m1, size = 50)
        evolve = Evolver(partitions = part1, tree = self.tree, branch_lengths = {"dist":"normal", "sd":0.1, "num_categories":"full"})
        mat, map = evolve._generate_transition_matrices(self.Q, self.t)
        
        np.testing.assert_array_almost_equal(map, np.arange(50), decimal = DECIMAL, err_msg = "Mapping incorrect for 'full' noisy branch lengths.")
        for entry in map:
            self.assertTrue( entry in range(len(mat)), msg = "Mapping values do not have corresponding matrices, for noisy branch lengths with num_categories='full'.")



    def test_evolver_noisy_branch_lengths_true_custom_n(self):
        
        part1 = Partition(models = self.m1, size = 50)
        n=20
        evolve = Evolver(partitions = part1, tree = self.tree, branch_lengths = {"dist":"normal", "sd":0.1, "num_categories":n})
        mat, map = evolve._generate_transition_matrices(self.Q, self.t)
        
        for entry in map:
            self.assertTrue( entry in range(len(mat)), msg = "Mapping values do not have corresponding matrices, for noisy branch lengths with n=custom value (20, here).")

          



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

    def test_evolver_setcode_aaaaanucleotide(self):
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
            
      
            
            
            
            
            
            
            
            