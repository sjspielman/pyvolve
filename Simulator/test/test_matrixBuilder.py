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


class matrixBuilder_baseClass_tests(unittest.TestCase):
    ''' 
        Set of unittests for simple base class functions in matrixBuilder.
        Functions tested here include isTI, orderNucleotidePair. 
        All other functions require full model specification, so they are tested elsewhere.
        Note: stateFreqs specification is required. To test this, all frequencies must be unique!.
    '''
    
    def setUp(self):
        # Do not rely on misc for codons in case something happens to it!
        self.codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        basicModel = Model()
        codonFreqs = [0.01617666, 0.00291771, 0.02664918, 0.02999061, 0.00717921, 0.00700012, 0.01435559, 0.0231568, 0.02403056, 0.00737008, 0.03185765, 0.0193576, 0.03277142, 0.02141258, 0.0127537, 0.00298803, 0.0256333, 0.02312437, 0.01861465, 0.01586447, 0.00373147, 0.02662654, 0.00082524, 0.00048916, 0.01191673, 0.00512658, 0.00050502, 0.01688169, 0.01843001, 0.00215437, 0.02659356, 0.02377742, 0.01169375, 0.00097256, 0.02937344, 0.00268204, 0.01414414, 0.02781933, 0.00070877, 0.02370841, 0.02984617, 0.01828081, 0.01002825, 0.00870788, 0.00728006, 0.02179328, 0.00379049, 0.01978996, 0.00443774, 0.01201798, 0.02030269, 0.01238501, 0.01279963, 0.02094385, 0.02810987, 0.00918507, 0.02880549, 0.0029311, 0.0237658, 0.03194712, 0.06148723]
        basicModel.params = {'stateFreqs': codonFreqs}
        self.baseObject = MatrixBuilder(basicModel)
        self.zero = 1e-8
        
    def test_matrixBuilder_baseClass_isTI(self):    
        ''' Test that transitions can be properly identified. '''
                
        self.assertTrue( self.baseObject.isTI('A', 'G'), msg = "matrixBuilder.isTI() does not think A -> G is a transition.")
        self.assertTrue( self.baseObject.isTI('G', 'A'), msg = "matrixBuilder.isTI() does not think G -> A is a transition.")
        self.assertTrue( self.baseObject.isTI('C', 'T'), msg = "matrixBuilder.isTI() does not think C -> T is a transition.")
        self.assertTrue( self.baseObject.isTI('T', 'C'), msg = "matrixBuilder.isTI() does not think C -> T is a transition.")
        self.assertFalse( self.baseObject.isTI('A', 'C'), msg = "matrixBuilder.isTI() mistakenly thinks A -> C is a transition.")
        self.assertFalse( self.baseObject.isTI('C', 'A'), msg = "matrixBuilder.isTI() mistakenly thinks C -> A is a transition.")
        self.assertFalse( self.baseObject.isTI('A', 'T'), msg = "matrixBuilder.isTI() mistakenly thinks A -> T is a transition.")
        self.assertFalse( self.baseObject.isTI('T', 'A'), msg = "matrixBuilder.isTI() mistakenly thinks T -> A is a transition.")
        self.assertFalse( self.baseObject.isTI('G', 'C'), msg = "matrixBuilder.isTI() mistakenly thinks G -> C is a transition.")
        self.assertFalse( self.baseObject.isTI('C', 'G'), msg = "matrixBuilder.isTI() mistakenly thinks C -> G is a transition.")
        self.assertFalse( self.baseObject.isTI('G', 'T'), msg = "matrixBuilder.isTI() mistakenly thinks G -> T is a transition.")
        self.assertFalse( self.baseObject.isTI('T', 'G'), msg = "matrixBuilder.isTI() mistakenly thinks T -> G is a transition.")

    def test_mechCodon_MatrixBuilder_isSyn(self):    
        ''' Test that synonymous vs nonsynymous changes can be properly identified. 
            Assumes that biopython is not broken. This is (theoretically...) a very safe assumption.
        '''
        for source in range(61):
            for target in range(61):
                sourceCodon = self.codons[source]
                targetCodon = self.codons[target]
                source_aa = str( Seq.Seq(sourceCodon, generic_dna).translate() )
                target_aa = str( Seq.Seq(targetCodon, generic_dna).translate() )
                if source_aa == target_aa:
                    self.assertTrue( self.baseObject.isSyn(source, target), msg = ("matrixBuilder.isSyn() does not think", source, " -> ", target, " is synonymous.") )
                else:
                    self.assertFalse( self.baseObject.isSyn(source, target), msg = ("matrixBuilder.isSyn() mistakenly thinks", source, " -> ", target, " is synonymous.") )



    def test_matrixBuilder_baseClass_orderNucleotidePair(self):    
        ''' Test that nucleotides can properly be ordered. ''' 
        self.assertEqual( self.baseObject.orderNucleotidePair('A', 'G'), 'AG', msg = "matrixBuilder.orderNucleotidePair can't order 'A', 'G' .")
        self.assertEqual( self.baseObject.orderNucleotidePair('G', 'A'), 'AG', msg = "matrixBuilder.orderNucleotidePair can't order 'G', 'A' .")
        self.assertEqual( self.baseObject.orderNucleotidePair('C', 'T'), 'CT', msg = "matrixBuilder.orderNucleotidePair can't order 'C', 'T' .")
        self.assertEqual( self.baseObject.orderNucleotidePair('T', 'C'), 'CT', msg = "matrixBuilder.orderNucleotidePair can't order 'T', 'C' .")
 


    def test_matrixBuilder_baseClass_getNucleotideDiff(self):
        ''' Test that nucleotide differences between codons can be identified properly. '''
        
        # No difference. Pos=F,T, Mul=F,T
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAA') ),              False       , msg = "matrixBuilder.getNucleotideDiff can't do same, pos=F." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAA'), True ),       (False, None), msg = "matrixBuilder.getNucleotideDiff can't do same, pos=T." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAA'), False, True ), '', msg = "matrixBuilder.getNucleotideDiff can't do same, pos=T." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAA'), True,  True ), '', msg = "matrixBuilder.getNucleotideDiff can't do same, pos=T." )
        
        # 1 difference, Pos=F,T, Mul=F,T
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('CAA')),        'AC',  msg = "matrixBuilder.getNucleotideDiff can't do one difference, position=F." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('ATA'), True), ('AT', 1),  msg = "matrixBuilder.getNucleotideDiff can't do one difference, pos=T, pos 2 change." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAT'), True), ('AT', 2),  msg = "matrixBuilder.getNucleotideDiff can't do one difference, pos=T, pos 3 change." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAT'), False, True),'AT',  msg = "matrixBuilder.getNucleotideDiff can't do one difference, pos=F, mul=T." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('AAT'), True, True), 'AT',  msg = "matrixBuilder.getNucleotideDiff can't do one difference, pos=T, mul=T." )
       
        # 2 difference. Pos=F,T, Mul=F,T
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('ATT')),        False, msg = "matrixBuilder.getNucleotideDiff can't do two differences, pos=F." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('ATT'), True), (False, None), msg = "matrixBuilder.getNucleotideDiff can't do two differences, pos=T." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('ATT'), False, True), 'ATAT', msg = "matrixBuilder.getNucleotideDiff can't do two differences, pos=F, mul=T." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('ATT'), True, True),  'ATAT', msg = "matrixBuilder.getNucleotideDiff can't do two differences, pos=T, mul=T." )


        # 3 difference. Pos=F,T, Mul=F,T.
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('TTT') ),        False       , msg = "matrixBuilder.getNucleotideDiff can't do fully distinct." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('TTT'), True ), (False, None), msg = "matrixBuilder.getNucleotideDiff can't do fully distinct, pos=T, mul=F." )
        self.assertEqual( self.baseObject.getNucleotideDiff( self.codons.index('AAA'), self.codons.index('TTT'), False, True ), 'ATATAT', msg = "matrixBuilder.getNucleotideDiff can't do fully distinct, pos=F, mul=T." )








class matrixBuilder_buildQ_tests(unittest.TestCase):
    ''' Test that the instantaneous matrix is being properly built for codon, mutsel, amino acid, and nucleotide models.
        The scaleMatrix function is implicitly tested within. This is not ideal and I intend to come back to this and separately test it.
    '''
    
    def setUp(self):
        self.dec = 8
    
    def test_matrixBuilder_buildQ_codon(self):    
        ''' Test proper Q construction for codon models.'''
        correctMatrix = np.loadtxt('matrixFiles/codonMatrix.txt')
        myFrequencies = np.loadtxt('matrixFiles/codonFreqs.txt')
        codonParams = {'stateFreqs': myFrequencies, 'alpha':1.0, 'beta':1.5, 'mu': {'AC': 1., 'AG': 2.5, 'AT': 1., 'CG': 1., 'CT': 2.5, 'GT': 1.}}
        codonModel = Model()
        codonModel.params = codonParams
        m = mechCodon_MatrixBuilder(codonModel)
        testMatrix = m.buildQ()
        np.testing.assert_array_almost_equal(correctMatrix, testMatrix, decimal = self.dec, err_msg = "Q improperly constructed for codon model.")












class matrixBuilder_mechCodon_MatrixBuilder_tests(unittest.TestCase):
    ''' 
        Set of unittests for the mechCodon_MatrixBuilder subclass of matrixBuilder.
        Functions tested here include isSyn, getCodonFreq, calcSynProb, calcNonsynProb, calcInstProb.
    '''
    
    def setUp(self):
        
        ################### DO NOT CHANGE ANY OF THESE EVER. #######################
        # Do not rely on misc for codons in case something happens to it!
        self.codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.nucleotides = ["A", "C", "G", "T"]
        self.codonFreqs = np.array( [0.01617666, 0.00291771, 0.02664918, 0.02999061, 0.00717921, 0.00700012, 0.01435559, 0.0231568, 0.02403056, 0.00737008, 0.03185765, 0.0193576, 0.03277142, 0.02141258, 0.0127537, 0.00298803, 0.0256333, 0.02312437, 0.01861465, 0.01586447, 0.00373147, 0.02662654, 0.00082524, 0.00048916, 0.01191673, 0.00512658, 0.00050502, 0.01688169, 0.01843001, 0.00215437, 0.02659356, 0.02377742, 0.01169375, 0.00097256, 0.02937344, 0.00268204, 0.01414414, 0.02781933, 0.00070877, 0.02370841, 0.02984617, 0.01828081, 0.01002825, 0.00870788, 0.00728006, 0.02179328, 0.00379049, 0.01978996, 0.00443774, 0.01201798, 0.02030269, 0.01238501, 0.01279963, 0.02094385, 0.02810987, 0.00918507, 0.02880549, 0.0029311, 0.0237658, 0.03194712, 0.06148723] )
        self.nucFreqs = np.array( [ 0.25, 0.20, 0.15, 0.40 ] )
        self.posNucFreqs = np.array([[ 0.25, 0.20, 0.15, 0.40 ],
                                     [ 0.15, 0.30, 0.05, 0.50 ],
                                     [ 0.35, 0.05, 0.30, 0.30 ]] )
                                
        muCodonParams = {'AG': 4.0, 'CT': 2.0, 'AC': 1.75, 'AT': 1.5, 'CG': 1.56, 'GT': 4.65}
        self.mycodon = Model()
        self.mycodon.params = { 'alpha':1.83, 'beta':5.7, 'mu': muCodonParams}
        self.zero = 1e-8
        ############################################################################
        
    
    def test_mechCodon_MatrixBuilder_GY94_getTargetFreq(self):
        ''' Test getTargetFreqs for GY94-style models '''
        self.mycodon.params['stateFreqs'] = self.codonFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )
        
        self.assertEqual( codonMatrix.getTargetFreq(0, None), 0.01617666, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for GY94, with position=None"))
        self.assertEqual( codonMatrix.getTargetFreq(1, 2), 0.00291771, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for GY94, with position=2"))
        
        

    def test_mechCodon_MatrixBuilder_MG94_nucFreqs_getTargetFreq(self):
        ''' Test getTargetFreqs for MG94-style models given 1D global nucleotide frequencies '''
        self.mycodon.params['stateFreqs'] = self.nucFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )

        self.assertEqual( codonMatrix.getTargetFreq(0, 0), 0.25, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with position=0"))
        self.assertEqual( codonMatrix.getTargetFreq(0, 1), 0.25, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with position=1"))
        self.assertEqual( codonMatrix.getTargetFreq(0, 2), 0.25, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with position=2"))




    def test_mechCodon_MatrixBuilder_MG94_posNucFreqs_getTargetFreq(self):
        ''' Test getTargetFreqs for MG94-style models given positional nucleotide frequencies '''
        self.mycodon.params['stateFreqs'] = self.posNucFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )

        self.assertEqual( codonMatrix.getTargetFreq(0, 0), 0.25, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with target=0, position=0"))
        self.assertEqual( codonMatrix.getTargetFreq(0, 1), 0.15, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with target=0, position=1"))
        self.assertEqual( codonMatrix.getTargetFreq(0, 2), 0.35, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with target=0, position=2"))
        self.assertEqual( codonMatrix.getTargetFreq(3, 2), 0.30, msg = ("mechCodon_MatrixBuilder.getTargetFreq() doesn't work for MG94, with target=0, position=2"))




    def test_mechCodon_MatrixBuilder_calcSynProb(self):
        ''' Test synonymous calculation. Note that, since target frequencies are already assigned, this function encompasses both GY94 and MG94. '''
        
        self.mycodon.params['stateFreqs'] = self.codonFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )

        # GCA -> GCT
        correctProb1 =  0.02370841 * 1.5 * 1.83
        self.assertTrue( abs(codonMatrix.calcSynProb(0.02370841, "AT") - correctProb1) < self.zero, msg = "mechCodon_MatrixBuiler.calcSynProb can't do GCA -> GCT.")
        # TTT -> TTC
        correctProb2 = 0.0237658 * 2.0 * 1.83
        self.assertTrue( abs(codonMatrix.calcSynProb(0.0237658, "CT") - correctProb2) < self.zero, msg = "mechCodon_MatrixBuiler.calcSynProb can't do TTT -> TTC.")
        # CAA -> CAG
        correctProb3 = 0.01861465 * 4.0 * 1.83 
        self.assertTrue( abs(codonMatrix.calcSynProb(0.01861465, "AG") - correctProb3) < self.zero, msg = "mechCodon_MatrixBuiler.calcSynProb can't do CAA -> CAG.")
        # CAG -> CAA (reverse of above.)
        correctProb4 = 0.0256333 * 4.0 * 1.83 
        self.assertTrue( abs(codonMatrix.calcSynProb(0.0256333, "AG") - correctProb4) < self.zero, msg = "mechCodon_MatrixBuiler.calcSynProb can't do CAG -> CAA.")



 
    def test_mechCodon_MatrixBuilder_calcNonsynProb(self):
        ''' Test that instantaneous substitution probabilities are properly calculated for nonsynonymous codons. 
            For tractability, just test a few nonsynonymous changes.
        '''
        self.mycodon.params['stateFreqs'] = self.codonFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )
        
        # TTA -> ATA
        correctProb1 =  0.03277142 * 1.5 * 5.7
        self.assertTrue( abs(codonMatrix.calcNonsynProb(0.03277142, "AT") - correctProb1) < self.zero, msg = "mechCodon_MatrixBuiler.calcNonsynProb can't do TTA -> ATA.")
        # CGT -> AGT
        correctProb2 =  0.0193576 * 1.75 * 5.7
        self.assertTrue( abs(codonMatrix.calcNonsynProb(0.0193576, "AC") - correctProb2) < self.zero, msg = "mechCodon_MatrixBuiler.calcNonsynProb can't do CGT -> AGT.")
        # TCC -> TGC
        correctProb3 =  0.02810987 * 1.56 * 5.7
        self.assertTrue( abs(codonMatrix.calcNonsynProb(0.02810987, "CG") - correctProb3) < self.zero, msg = "mechCodon_MatrixBuiler.calcNonsynProb can't do TCC -> TGC.")
        # TGC -> TCC, reverse of above.
        correctProb4 =  0.01238501 * 1.56 * 5.7
        self.assertTrue( abs(codonMatrix.calcNonsynProb(0.01238501, "CG") - correctProb4) < self.zero, msg = "mechCodon_MatrixBuiler.calcNonsynProb can't do TGC -> TCC.")



    def test_mechCodon_MatrixBuilder_calcInstProb_GY94(self):    
        ''' Test that substitution probabilities are properly calculated for GY94-style models
            Conduct tests for - no change, two changes, three changes, synonymous, nonsynonymous.
        '''
        self.mycodon.params['stateFreqs'] = self.codonFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )
        
        # Test no change, two changes, three changes. All should be 0
        self.assertTrue( abs(codonMatrix.calcInstProb(7, 7) - 0.) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb doesn't return 0 for same codon substitution when GY94.")
        self.assertTrue( abs(codonMatrix.calcInstProb(7, 8) - 0.) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb doesn't return 0 for two nucleotide changes when GY94.")
        self.assertTrue( abs(codonMatrix.calcInstProb(7, 24) - 0.) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb doesn't return 0 for three nucleotide changes when GY94.")
        
        # Synonymous. GAG -> GAA
        correctProbSyn = 0.01169375 * 1.83 * 4.0
        self.assertTrue( abs(codonMatrix.calcInstProb(34, 32) - correctProbSyn) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb wrong for GAG -> GAA (synonymous) when GY94.")

        # Nonsynonymous. TCG -> ACG
        correctProbNonsyn = 0.01435559 * 5.7 * 1.5
        self.assertTrue( abs(codonMatrix.calcInstProb(52, 6) - correctProbNonsyn) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb wrong for TCG -> ACG (nonsynonymous) when GY94.")
        
    
    
    
    def test_mechCodon_MatrixBuilder_calcInstProb_MG94(self):    
        ''' Test that substitution probabilities are properly calculated for MG94-style models
            Conduct tests for synonymous, nonsynonymous (the test function for GY94 has already tested for 0,2,3 nuc changes).
        '''
        self.mycodon.params['stateFreqs'] = self.posNucFreqs
        codonMatrix = mechCodon_MatrixBuilder( self.mycodon )

        # Synonymous, position 3. GAG -> GAA
        correctProbSyn = 0.35 * 1.83 * 4.0
        self.assertTrue( abs(codonMatrix.calcInstProb(34, 32) - correctProbSyn) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb wrong for GAG -> GAA (synonymous, pos 3) when MG94.")

        # Nonsynonymous, position 1. ACG -> TCG
        correctProbNonsyn = 0.4 * 5.7 * 1.5
        self.assertTrue( abs(codonMatrix.calcInstProb(6, 52) - correctProbNonsyn) < self.zero, msg = "mechCodon_MatrixBuilder.calcInstProb wrong for TCG -> ACG (nonsynonymous, pos 1) when MG94.")     







class matrixBuilder_aminoAcid_MatrixBuilder_tests(unittest.TestCase):
    ''' 
        Set of unittests for the aminoAcid_MatrixBuilder subclass of matrixBuilder, which deals with empirical amino acid models.
        Functions tested here include initEmpiricalMatrix, calcInstProb.
        We are going to just test everything with the LG matrix.
    '''
    
    def setUp(self):
        ################### DO NOT CHANGE ANY OF THESE EVER. #######################
        self.amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        aminoFreqs = [0.03467759, 0.02767874, 0.05246165, 0.03385361, 0.07774578, 0.06417143, 0.01284807, 0.09590641, 0.04063569, 0.04324441, 0.05862815, 0.0198158, 0.07119442, 0.08725553, 0.0105108, 0.05877746, 0.03735058, 0.04630849, 0.05067367, 0.0762617]
        self.lg_mat = np.array([[0.000000000000000000e+00, 2.489084000000000074e+00, 3.951439999999999952e-01, 1.038545000000000051e+00, 2.537010000000000098e-01, 2.066040000000000099e+00, 3.588580000000000103e-01, 1.498299999999999910e-01, 5.365180000000000504e-01, 3.953369999999999940e-01, 1.124034999999999895e+00, 2.768180000000000085e-01, 1.177651000000000003e+00, 9.698940000000000339e-01, 4.250929999999999986e-01, 4.727181999999999995e+00, 2.139501000000000097e+00, 2.547870000000000079e+00, 1.807169999999999888e-01, 2.189589999999999870e-01],[2.489084000000000074e+00, 0.000000000000000000e+00, 6.255600000000000049e-02, 3.498999999999999940e-03, 1.105250999999999983e+00, 5.692650000000000210e-01, 6.405429999999999735e-01, 3.206269999999999953e-01, 1.326600000000000001e-02, 5.940069999999999517e-01, 8.936800000000000299e-01, 5.287680000000000158e-01, 7.538200000000000456e-02, 8.480799999999999450e-02, 5.345509999999999984e-01, 2.784478000000000009e+00, 1.143480000000000052e+00, 1.959290999999999894e+00, 6.701279999999999459e-01, 1.165532000000000012e+00],[3.951439999999999952e-01, 6.255600000000000049e-02, 0.000000000000000000e+00, 5.243870000000000253e+00, 1.741600000000000092e-02, 8.449259999999999549e-01, 9.271139999999999937e-01, 1.068999999999999985e-02, 2.829590000000000161e-01, 1.507600000000000086e-02, 2.554800000000000126e-02, 5.076149000000000022e+00, 3.944559999999999733e-01, 5.233860000000000179e-01, 1.239539999999999947e-01, 1.240275000000000016e+00, 4.258600000000000163e-01, 3.796700000000000075e-02, 2.988999999999999990e-02, 1.351070000000000049e-01],[1.038545000000000051e+00, 3.498999999999999940e-03, 5.243870000000000253e+00, 0.000000000000000000e+00, 1.881100000000000133e-02, 3.488470000000000182e-01, 4.238810000000000078e-01, 4.426499999999999879e-02, 1.807177000000000033e+00, 6.967299999999999882e-02, 1.737350000000000005e-01, 5.417119999999999713e-01, 4.194089999999999763e-01, 4.128591000000000122e+00, 3.639700000000000157e-01, 6.119729999999999892e-01, 6.045449999999999990e-01, 2.450340000000000018e-01, 7.785200000000000453e-02, 1.200370000000000048e-01],[2.537010000000000098e-01, 1.105250999999999983e+00, 1.741600000000000092e-02, 1.881100000000000133e-02, 0.000000000000000000e+00, 8.958599999999999897e-02, 6.821390000000000509e-01, 1.112727000000000022e+00, 2.391799999999999829e-02, 2.592691999999999997e+00, 1.798853000000000035e+00, 8.952499999999999347e-02, 9.446400000000000630e-02, 3.585499999999999798e-02, 5.272199999999999803e-02, 3.618190000000000017e-01, 1.650010000000000088e-01, 6.546830000000000149e-01, 2.457120999999999889e+00, 7.803901999999999894e+00],[2.066040000000000099e+00, 5.692650000000000210e-01, 8.449259999999999549e-01, 3.488470000000000182e-01, 8.958599999999999897e-02, 0.000000000000000000e+00, 3.114839999999999831e-01, 8.704999999999999197e-03, 2.966360000000000108e-01, 4.426100000000000173e-02, 1.395379999999999954e-01, 1.437645000000000062e+00, 1.969609999999999972e-01, 2.679590000000000027e-01, 3.901919999999999833e-01, 1.739989999999999926e+00, 1.298360000000000070e-01, 7.670100000000000529e-02, 2.684909999999999797e-01, 5.467899999999999844e-02],[3.588580000000000103e-01, 6.405429999999999735e-01, 9.271139999999999937e-01, 4.238810000000000078e-01, 6.821390000000000509e-01, 3.114839999999999831e-01, 0.000000000000000000e+00, 1.088820000000000066e-01, 6.972639999999999949e-01, 3.663170000000000037e-01, 4.424719999999999764e-01, 4.509237999999999857e+00, 5.088510000000000533e-01, 4.813505000000000145e+00, 2.426600999999999786e+00, 9.900120000000000031e-01, 5.842619999999999481e-01, 1.190129999999999938e-01, 5.970539999999999736e-01, 5.306834000000000273e+00],[1.498299999999999910e-01, 3.206269999999999953e-01, 1.068999999999999985e-02, 4.426499999999999879e-02, 1.112727000000000022e+00, 8.704999999999999197e-03, 1.088820000000000066e-01, 0.000000000000000000e+00, 1.590689999999999882e-01, 4.145067000000000057e+00, 4.273607000000000156e+00, 1.915030000000000066e-01, 7.828100000000000336e-02, 7.285400000000000209e-02, 1.269909999999999928e-01, 6.410499999999999532e-02, 1.033738999999999963e+00, 6.491069999999999895e-01, 1.116599999999999954e-01, 2.325230000000000075e-01],[5.365180000000000504e-01, 1.326600000000000001e-02, 2.829590000000000161e-01, 1.807177000000000033e+00, 2.391799999999999829e-02, 2.966360000000000108e-01, 6.972639999999999949e-01, 1.590689999999999882e-01, 0.000000000000000000e+00, 1.375000000000000111e-01, 6.566039999999999655e-01, 2.145077999999999818e+00, 3.903220000000000023e-01, 3.234293999999999780e+00, 6.326067000000000107e+00, 7.486829999999999874e-01, 1.136862999999999957e+00, 1.852020000000000055e-01, 4.990599999999999897e-02, 1.319319999999999937e-01],[3.953369999999999940e-01, 5.940069999999999517e-01, 1.507600000000000086e-02, 6.967299999999999882e-02, 2.592691999999999997e+00, 4.426100000000000173e-02, 3.663170000000000037e-01, 4.145067000000000057e+00, 1.375000000000000111e-01, 0.000000000000000000e+00, 6.312357999999999691e+00, 6.842700000000000171e-02, 2.490600000000000036e-01, 5.824570000000000025e-01, 3.018480000000000052e-01, 1.822870000000000046e-01, 3.029359999999999831e-01, 1.702744999999999953e+00, 6.196319999999999606e-01, 2.996480000000000254e-01],[1.124034999999999895e+00, 8.936800000000000299e-01, 2.554800000000000126e-02, 1.737350000000000005e-01, 1.798853000000000035e+00, 1.395379999999999954e-01, 4.424719999999999764e-01, 4.273607000000000156e+00, 6.566039999999999655e-01, 6.312357999999999691e+00, 0.000000000000000000e+00, 3.710040000000000004e-01, 9.984899999999999332e-02, 1.672568999999999972e+00, 4.841329999999999800e-01, 3.469599999999999906e-01, 2.020366000000000106e+00, 1.898717999999999906e+00, 6.961749999999999883e-01, 4.813060000000000116e-01],[2.768180000000000085e-01, 5.287680000000000158e-01, 5.076149000000000022e+00, 5.417119999999999713e-01, 8.952499999999999347e-02, 1.437645000000000062e+00, 4.509237999999999857e+00, 1.915030000000000066e-01, 2.145077999999999818e+00, 6.842700000000000171e-02, 3.710040000000000004e-01, 0.000000000000000000e+00, 1.617869999999999864e-01, 1.695751999999999926e+00, 7.518780000000000463e-01, 4.008358000000000310e+00, 2.000678999999999874e+00, 8.368799999999999850e-02, 4.537599999999999967e-02, 6.120250000000000412e-01],[1.177651000000000003e+00, 7.538200000000000456e-02, 3.944559999999999733e-01, 4.194089999999999763e-01, 9.446400000000000630e-02, 1.969609999999999972e-01, 5.088510000000000533e-01, 7.828100000000000336e-02, 3.903220000000000023e-01, 2.490600000000000036e-01, 9.984899999999999332e-02, 1.617869999999999864e-01, 0.000000000000000000e+00, 6.242940000000000156e-01, 3.325330000000000230e-01, 1.338132000000000099e+00, 5.714679999999999760e-01, 2.965010000000000145e-01, 9.513099999999999334e-02, 8.961299999999999821e-02],[9.698940000000000339e-01, 8.480799999999999450e-02, 5.233860000000000179e-01, 4.128591000000000122e+00, 3.585499999999999798e-02, 2.679590000000000027e-01, 4.813505000000000145e+00, 7.285400000000000209e-02, 3.234293999999999780e+00, 5.824570000000000025e-01, 1.672568999999999972e+00, 1.695751999999999926e+00, 6.242940000000000156e-01, 0.000000000000000000e+00, 2.807907999999999848e+00, 1.223827999999999916e+00, 1.080135999999999985e+00, 2.103319999999999912e-01, 2.361989999999999923e-01, 2.573360000000000092e-01],[4.250929999999999986e-01, 5.345509999999999984e-01, 1.239539999999999947e-01, 3.639700000000000157e-01, 5.272199999999999803e-02, 3.901919999999999833e-01, 2.426600999999999786e+00, 1.269909999999999928e-01, 6.326067000000000107e+00, 3.018480000000000052e-01, 4.841329999999999800e-01, 7.518780000000000463e-01, 3.325330000000000230e-01, 2.807907999999999848e+00, 0.000000000000000000e+00, 8.581509999999999971e-01, 5.789870000000000294e-01, 1.708870000000000111e-01, 5.936069999999999958e-01, 3.144399999999999973e-01],[4.727181999999999995e+00, 2.784478000000000009e+00, 1.240275000000000016e+00, 6.119729999999999892e-01, 3.618190000000000017e-01, 1.739989999999999926e+00, 9.900120000000000031e-01, 6.410499999999999532e-02, 7.486829999999999874e-01, 1.822870000000000046e-01, 3.469599999999999906e-01, 4.008358000000000310e+00, 1.338132000000000099e+00, 1.223827999999999916e+00, 8.581509999999999971e-01, 0.000000000000000000e+00, 6.472279000000000337e+00, 9.836899999999999811e-02, 2.488619999999999999e-01, 4.005469999999999864e-01],[2.139501000000000097e+00, 1.143480000000000052e+00, 4.258600000000000163e-01, 6.045449999999999990e-01, 1.650010000000000088e-01, 1.298360000000000070e-01, 5.842619999999999481e-01, 1.033738999999999963e+00, 1.136862999999999957e+00, 3.029359999999999831e-01, 2.020366000000000106e+00, 2.000678999999999874e+00, 5.714679999999999760e-01, 1.080135999999999985e+00, 5.789870000000000294e-01, 6.472279000000000337e+00, 0.000000000000000000e+00, 2.188158000000000047e+00, 1.408250000000000057e-01, 2.458410000000000040e-01],[2.547870000000000079e+00, 1.959290999999999894e+00, 3.796700000000000075e-02, 2.450340000000000018e-01, 6.546830000000000149e-01, 7.670100000000000529e-02, 1.190129999999999938e-01, 6.491069999999999895e-01, 1.852020000000000055e-01, 1.702744999999999953e+00, 1.898717999999999906e+00, 8.368799999999999850e-02, 2.965010000000000145e-01, 2.103319999999999912e-01, 1.708870000000000111e-01, 9.836899999999999811e-02, 2.188158000000000047e+00, 0.000000000000000000e+00, 1.895100000000000118e-01, 2.493130000000000068e-01],[1.807169999999999888e-01, 6.701279999999999459e-01, 2.988999999999999990e-02, 7.785200000000000453e-02, 2.457120999999999889e+00, 2.684909999999999797e-01, 5.970539999999999736e-01, 1.116599999999999954e-01, 4.990599999999999897e-02, 6.196319999999999606e-01, 6.961749999999999883e-01, 4.537599999999999967e-02, 9.513099999999999334e-02, 2.361989999999999923e-01, 5.936069999999999958e-01, 2.488619999999999999e-01, 1.408250000000000057e-01, 1.895100000000000118e-01, 0.000000000000000000e+00, 3.151815000000000033e+00],[2.189589999999999870e-01, 1.165532000000000012e+00, 1.351070000000000049e-01, 1.200370000000000048e-01, 7.803901999999999894e+00, 5.467899999999999844e-02, 5.306834000000000273e+00, 2.325230000000000075e-01, 1.319319999999999937e-01, 2.996480000000000254e-01, 4.813060000000000116e-01, 6.120250000000000412e-01, 8.961299999999999821e-02, 2.573360000000000092e-01, 3.144399999999999973e-01, 4.005469999999999864e-01, 2.458410000000000040e-01, 2.493130000000000068e-01, 3.151815000000000033e+00, 0.000000000000000000e+00]])
        model = Model()
        model.params = {'stateFreqs': aminoFreqs, 'aaModel': 'LG'}
        self.aaMatrix = aminoAcid_MatrixBuilder(model)
        self.dec = 8
        ############################################################################
    
   
    def test_aminoAcid_MatrixBuilder_initEmpiricalMatrix(self):
        ''' Tests that class initialization properly imported the empirical replacement matrix. '''
        np.testing.assert_array_almost_equal(self.aaMatrix.empMat, self.lg_mat, decimal = self.dec, err_msg = "aminoAcid_MatrixBuilder.initEmpiricalMatrix doesn't return empirical matrix properly.")

        
    def test_aminoAcid_MatrixBuilder_calcInstProb(self):
        ''' Test calcInstProb function for amino acid class. '''
        correctProb = 3.499e-03 * 0.03385361
        self.assertEqual(self.aaMatrix.calcInstProb(1, 3), correctProb, msg = "aminoAcid_MatrixBuilder.calcInstProb fail.")





class matrixBuilder_nucleotide_MatrixBuilder_tests(unittest.TestCase):
    ''' 
        Set of unittests for the nucleotide_MatrixBuilder subclass of matrixBuilder.
        Functions tested here include getNucleotideFreq, calcInstProb.
    '''
    
    def setUp(self):
        ################### DO NOT CHANGE ANY OF THESE EVER. #######################
        self.nucleotides = ['A', 'C', 'G', 'T']
        muParams = {'AG':0.1, 'GA':0.15, 'CT':0.125, 'TC':0.125, 'AC': 0.08, 'CA':0.17, 'AT':0.05, 'TA':0.075, 'CG':0.125, 'GC':0.125, 'GT':0.13, 'TG':0.12}
        model = Model()
        model.params = {'stateFreqs': [0.34, 0.21, 0.27, 0.18], 'mu': muParams}
        self.nucMatrix = nucleotide_MatrixBuilder(model)
        ############################################################################
    
   
    def test_nucleotide_MatrixBuilder_calcInstProb(self):
        ''' Test function to retrieve instantaneous substitution probability between nucleotides. Just test a few. '''
        correctAT = 0.18 * 0.05
        self.assertEqual(self.nucMatrix.calcInstProb(0, 3), correctAT, msg = "nucleotideMatrix.calcInstProb doesn't properly work for A->T.")
        correctTA = 0.34 * 0.075
        self.assertEqual(self.nucMatrix.calcInstProb(3, 0), correctTA, msg = "nucleotideMatrix.calcInstProb doesn't properly work for T->A.")
        correctCA = 0.34 * 0.17
        self.assertEqual(self.nucMatrix.calcInstProb(1, 0), correctCA, msg = "nucleotideMatrix.calcInstProb doesn't properly work for C->A.")









class matrixBuilder_mutSel_MatrixBuilder_tests(unittest.TestCase):
    ''' 
        Set of unittests for the mutSel_MatrixBuilder subclass of matrixBuilder.
        Functions tested here include isSyn, getCodonFreq, calcSynProb, calcNonsynProb, calcInstProb.
    '''
    def setUp(self):
        ################### DO NOT CHANGE ANY OF THESE EVER. #######################
        self.codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        codonFreqs = np.array( [0, 0.04028377, 0.02664918, 0, 0.00717921, 0.00700012, 0.0231568, 0.0231568, 0.02403056, 0.00737008, 0.03185765, 0.0193576, 0.03277142, 0.02141258, 0.0127537, 0.00298803, 0.0256333, 0.02312437, 0.01861465, 0.01586447, 0.00373147, 0.02662654, 0.00082524, 0.00048916, 0.01191673, 0.00512658, 0.00050502, 0.01688169, 0.01843001, 0.00215437, 0.02659356, 0.02377742, 0.01169375, 0.00097256, 0.02937344, 0.00268204, 0.01414414, 0.02781933, 0.00070877, 0.02370841, 0.02984617, 0.01828081, 0.01002825, 0.00870788, 0.00728006, 0.02179328, 0.00379049, 0.01978996, 0.00443774, 0.01201798, 0.02030269, 0.01238501, 0.01279963, 0.02094385, 0.02810987, 0.00918507, 0.02880549, 0.0029311, 0.0237658, 0.03194712, 0.06148723] )
        muParams = {'AG':0.125, 'GA':0.125, 'CT':0.125, 'TC':0.125, 'AC': 0.125, 'CA':0.125, 'AT':0.125, 'TA':0.125, 'CG':0.125, 'GC':0.125, 'GT':0.13, 'TG':0.12}
        model = Model()
        model.params = {'stateFreqs': codonFreqs, 'mu': muParams}
        self.mutSelMatrix = mutSel_MatrixBuilder( model )
        self.zero = 1e-8
        ############################################################################

    def test_mutSel_MatrixBuilder_calcSubstitutionProb(self):    
        ''' Test function calcSubstitutionProb for mutation-selection model subclass.
            Test where target has 0 freq, source has 0 freq, both have 0 freq, they have equal freq, they have different freq.
        '''

        # Target and/or source have 0 or equal frequency. Assertions should be raised.
        self.assertRaises(AssertionError, self.mutSelMatrix.calcSubstitutionProb, 0., 0., 0.65, 0.32)
        self.assertRaises(AssertionError, self.mutSelMatrix.calcSubstitutionProb, 0., 0.084572, 0.82, 0.71)    
        self.assertRaises(AssertionError, self.mutSelMatrix.calcSubstitutionProb, 0.10599277, 0., 0.111, 0.0099)
        self.assertRaises(AssertionError, self.mutSelMatrix.calcSubstitutionProb, 0.0756, 0.0756, 0.982, 0.00234)

        # Target and source have different frequencies
        self.assertTrue( abs(self.mutSelMatrix.calcSubstitutionProb(0.367, 0.02345, 0.09, 0.06) - 0.02237253623) < self.zero, msg = "mutSel_MatrixBuilder.calcSubstitutionProb fails when target and source have different frequencies.")

    def test_mutSel_MatrixBuilder_calcInstProb(self):    
        ''' Test function calcInstProb for mutation-selection model subclass.
            Test for one or both have freq 0, have equal freq, have no changes, have multiple changes, and finally, have different freq.
        '''

        # Target and/or source have 0 or equal frequency should return 0.
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(0, 3) - 0.) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when both codons have 0 frequency.")
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(0, 4) - 0.) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when source codon only has 0 frequency.")
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(7, 3) - 0.) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when target codon only has 0 frequency.")
        
        # Equal frequency should return forward mutation rate
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(6, 7) - 0.13) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when codons have equal frequency.")
        
        # Too few or too many changes
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(6, 6) - 0.) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when source and target are the same.")
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(6, 53) - 0.) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when two changes between codons.")
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(0, 60) - 0.) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when three changes between codons.")

        # Different frequencies.
        self.assertTrue( abs(self.mutSelMatrix.calcInstProb(56, 55) - 0.0612161452749) < self.zero, msg = "mutSel_MatrixBuilder.calcInstProb wrong when codons have equal frequency.")





if __name__ == '__main__':
    run_tests = unittest.TextTestRunner()
    
    print "Testing aminoAcids_MatrixBuilder, a subclass of the parent matrixBuilder"
    test_suite_aminoAcidMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_aminoAcid_MatrixBuilder_tests)
    run_tests.run(test_suite_aminoAcidMatrix)

    print "Testing the simple functions in the base class matrixBuilder"
    test_suite_baseMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_baseClass_tests)
    run_tests.run(test_suite_baseMatrix)
    
    print "Testing mechCodon_MatrixBuilder, a subclass of the parent matrixBuilder"
    test_suite_codonMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_mechCodon_MatrixBuilder_tests)
    run_tests.run(test_suite_codonMatrix)

    print "Testing mutSel_MatrixBuilder, a subclass of the parent matrixBuilder"
    test_suite_mutSelMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_mutSel_MatrixBuilder_tests)
    run_tests.run(test_suite_mutSelMatrix)
    
    print "Testing nucleotide_MatrixBuilder, a subclass of the parent matrixBuilder"
    test_suite_nucleotideMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_nucleotide_MatrixBuilder_tests)
    run_tests.run(test_suite_nucleotideMatrix)

    print "Testing buildQ function of matrixBuilder for codon model"
    test_suite_buildQ = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_buildQ_tests)
    run_tests.run(test_suite_buildQ)
