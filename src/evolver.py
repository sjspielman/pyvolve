#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Evolve sequences along a phylogeny.
'''



import numpy as np
from scipy import linalg
import random as rn
from time import strftime
from misc import ZERO, Genetics, Model
MOLECULES = Genetics()
        
class Evolver(object):
    ''' 
        Class to evolve sequences along a phylogeny.
    '''
    def __init__(self, partitions, root_model):
        
        
        self._partitions    = partitions   # this should be a list of tuples. Each tuple is (length, {flag:model, ...}). 
        self._number_partitions = len(self._partitions)
        assert ( self._number_partitions >= 1), "You have nothing to evolve. Partitions, please!"
        self._root_model = root_model # some checking should probably be done here.
        self.alndict = {}
        self._seq_length = 0
        for i in range(self._number_partitions):
            self._seq_length += self._partitions[i][0]   
        self._set_code()

        
    def _set_code(self):
        ''' 
            Assign genetic code.
        '''   
        
        dim = self._partitions[0][1].values()[0].Q.shape[0] # Go python, go!!
        if dim == 4:
            self._code = MOLECULES.nucleotides
        elif dim == 20:
            self._code = MOLECULES.amino_acids
        elif dim == 61:
            self._code = MOLECULES.codons
        else:
            raise AssertionError("wut.")
            
                        
    def _sequence_to_integer(self, entry):
        ''' 
            Convert a dna/protein character to its appropriate integer (index in self._code).
        '''
        return self._code.index(entry)
    
    def _integer_to_sequence(self, index):
        '''
            Convert an integer (index in self._code) to its appropriate dna/protein character.
        '''
        return self._code[index]
 
 
    def _intseq_to_string(self, intseq):
        ''' 
            Convert a full sequence coded as integers (indices in self._code) and convert to a dna/protein string.
        '''
        stringseq = ''
        for i in intseq:
            stringseq += self._integer_to_sequence(i)
        return stringseq   


    def write_sequences(self, **kwargs):
        ''' 
            Write resulting sequence alignment (self.alndict) to a file in fasta format.
            NOTE: THIS FUNCTION WILL BE REPLACED BY A MUCH MORE GENERAL ONE IN THE NEXT AND/OR FIRST RELEASE.
        '''
        
        outfile  = kwargs.get("outfile", "seqs_"+strftime("%m.%d.%H.%M.%S")+".fasta")  
        out_handle=open(outfile, 'w')
        for entry in self.alndict:
            seq = self._intseq_to_string(self.alndict[entry])
            out_handle.write(">"+entry+"\n"+seq+"\n")
        out_handle.close()  
  
        
            
    def _generate_prob_from_unif(self, prob_array):
        ''' 
            Sample a sequence (nuc,aa,or codon), and return an integer for the sequence chosen from a uniform distribution.
            Arugment "prob_array_ is any list and/or numpy array of probabilities which sum to 1.
        '''
        
        assert ( abs(np.sum(prob_array) - 1.) < ZERO), "Probabilities do not sum to 1. Cannot generate a new sequence."
        r = rn.uniform(0,1)
        i = 0
        sum = prob_array[i]
        while sum < r:
            i += 1
            sum += prob_array[i]
        return i     

        
    def _generate_root_seq(self):
        ''' 
            Generate a root sequence based on the stationary frequencies, for each partition and corresponding model (if they differ).
            Return a complete root sequence (again, coded in integers).
        '''
        
        root_sequence = np.empty(self._seq_length, dtype=int)
        index = 0
        for n in range(self._number_partitions):
            partlen = self._partitions[n][0]
            freqs = self._partitions[n][1][self._root_model].params['state_freqs']
            for j in range(partlen):
                root_sequence[index] = self._generate_prob_from_unif(freqs)
                index += 1
        return root_sequence

      
        
    def _check_parent_branch(self, node, parent_seq, parent_model):
        ''' 
            Function ensures that, for a given node we'd like to evolve to, an appropriate branch length exists. 
            If the branch length is acceptable, an evolutionary model is then assigned to the node.
        '''
        assert (parent_seq != None), "\n\nThere is no parent sequence from which to evolve!"

        branch_model = node.model_flag
        if branch_model is None:
            node.model_flag = parent_model
        assert (branch_length > 0), "\n\nYour tree has a negative branch length. I'm going to quit now."
        return float(node.branch_length), node.model_flag        
        
        
        
    def simulate(self, current_node, parent_node = None):
        ''' 
            Function to traverse a Tree object recursively and simulate sequences.
            Arguments:
                1. "current_node" is the node (either internal node or leaf) we are evolving TO
                2. "parent_node" is the node we are evolving FROM. Default of None is only called when the root sequence is not yet made.
        '''

        # We are at the base and must generate root sequence
        if (parent_node is None):
            current_node.seq = self._generate_root_seq() 
            current_node.model_flag = self._root_model 
        else:
            self.evolve_branch(current_node, parent_node)
            
            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                self.simulate(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[current_node.name]=current_node.seq
            
            
            
            
            
    def evolve_branch(self, current_node, parent_node):
        ''' 
            Function to evolve a given sequence during tree traversal.
            Arguments:
                1. "current_node" is the node (either internal node or leaf) we are evolving TO
                2. "parent_node" is the node we are evolving FROM.
        '''
    
        # Ensure brank length ok, parent sequence exists, and model is assigned.
        parent_seq = parent_node.seq
        parent_model = parent_node.model_flag
        branch_length, branch_model = self._check_parent_branch(current_node, parent_seq, parent_model)

        # Evolve only if branch length is greater than 0. 
        if branch_length <= ZERO:
            new_seq = parent_seq
        else:
            # Evolve each partition, i, and add onto new_seq as we go
            new_seq = np.empty(self._seq_length, dtype=int)
            index = 0
            for i in range(self._number_partitions):
            
                # set the length and the instantaneous rate matrix for this partition at this node
                seqlen  = self._partitions[i][0]
                inst_matrix = self._partitions[i][1][branch_model].Q
                
                # Generate probability matrix for evolution along this branch and assert correct
                Qt = np.multiply(inst_matrix, branch_length)
                prob_matrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
                for i in range(len(self._code)):
                    assert( abs(np.sum(prob_matrix[i]) - 1.) < ZERO ), "Row in P(t) matrix does not sum to 1."
    
                # Move along parentSeq and evolve. 
                for j in range(seqlen):
                    new_seq[index] = self._generate_prob_from_unif( prob_matrix[parent_seq[index]] )
                    index+=1
                             
        # Attach final sequence to node
        current_node.seq = new_seq 