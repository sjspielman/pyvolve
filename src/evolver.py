# SJS. NOTE - AS STANDS, DOES NOT SUPPORT INDELS.
import numpy as np
from scipy import linalg
import random as rn
from time import strftime
from misc import ZERO, Genetics, Model
MOLECULES = Genetics()
        
class Evolver(object):
    def __init__(self, partitions, root_model):
        
        
        self._partitions    = partitions   # this should be a list of tuples. Each tuple is (length, {flag:model, ...}). 
        self._numparts      = len(self._partitions)
        assert ( self._number_partitions >= 1), "You have nothing to evolve. Partitions, please!"
        self._root_model = root_model # some checking should probably be done here.
        self.alndict = {}
        self._seq_length = 0
        for i in range(self._number_partitions):
            self._seq_length += self._partitions[i][0]   
        self.set_code()

        
    def _set_code(self):
        ''' Determine whether we have nuc, amino, or codon alphabet and assign a genetic code. '''   
        # Go python, go!!
        dim = self._partitions[0][1].values()[0].Q.shape[0]
        if dim == 4:
            self._code = MOLECULES.nucleotides
        elif dim == 20:
            self._code = MOLECULES.amino_acids
        elif dim == 61:
            self._code = MOLECULES.codons
        else:
            raise AssertionError("wut.")
            
                        
    def _sequence_to_integer(self, entry):
        ''' Take a DNA/PROT value and return its integer index '''
        return self._code.index(entry)
    
    def _integer_to_sequence(self, index):
        ''' Take a DNA/PROT index and return its corresponding molecule ''' 
        return self._code[index]
 
    def _intseq_to_string(self, intseq):
        ''' Take a sequence coded as ints and turn to actual molecule string '''
        stringseq = ''
        for i in intseq:
            stringseq += self._integer_to_sequence(i)
        return stringseq   

    def _check_parent_branch(self, node, parent_seq, parent_model):
        ''' Check that the parent sequence exists, branch length is reasonable, and assign a model. ''' 
        assert (parent_seq != None), "There is no parent sequence."

        branch_model = node.model_flag
        if branch_model is None:
            node.model_flag = parent_model
        branch_length = float( node.branch_length )
        assert (branch_length >= ZERO), "Branch length is negative. Must be >= 0."
        return branch_length, node.model_flag    
        
            
    def _generate_prob_from_unif(self, prob_array):
        ''' Sample a sequence letter (nuc,aa,or codon). prob_array can be any list/numpy array of probabilities that sum to 1.'''
        assert ( abs(np.sum(prob_array) - 1.) < ZERO), "Probabilities do not sum to 1. Cannot generate a new sequence."
        r = rn.uniform(0,1)
        i=0
        sum=prob_array[i]
        while sum < r:
            i+=1
            sum+=prob_array[i]
        return i     

        
    def _generate_root_seq(self):
        ''' Select starting sequence based on state frequencies, for each partition, and return full root sequence. '''
        
        root_sequence = np.empty(self._seq_length, dtype=int)
        index = 0
        for n in range(self._number_partitions):
            partlen = self._partitions[n][0]
            freqs = self._partitions[n][1][self._root_model].params['state_freqs']
            for j in range(partlen):
                root_sequence[index] = self._generate_prob_from_unif(freqs)
                index += 1
       return root_sequence



    def write_sequences(self, **kwargs):
        ''' Write resulting sequences to a file, currently only in fasta format.'''
        
        outfile  = kwargs.get("outfile", "seqs_"+strftime("%m.%d.%H.%M.%S")+".fasta")  
        out_handle=open(outfile, 'w')
        for entry in self.alndict:
            seq = self._intseq_to_string(self.alndict[entry])
            out_handle.write(">"+entry+"\n"+seq+"\n")
        out_handle.close()    
        
        
        
        
    def simulate(self, current_node, parent_node = None):
        ''' Traverse the tree and simulate. '''

        # We are at the base and must generate root sequence
        if (parent_node is None):
            current_node.seq = self._generate_root_seq() 
            current_node.model_flag = self.root_model 
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
        ''' Crux function to evolve sequences along a branch. INDELS ARE NOT SUPPORTED.'''
    
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
                seqlen  = self.parts[i][0]
                inst_matrix = self.parts[i][1][branch_bodel].Q
                
                # Generate probability matrix for evolution along this branch and assert correct
                Qt = np.multiply(inst_matrix, branch_length)
                prob_matrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
                for i in range(len(self._code)):
                    assert( abs(np.sum(prob_matrix[i]) - 1.) < ZERO ), "Row in P(t) matrix does not sum to 1."
    
                # Move along parentSeq and evolve. 
                for j in range(seqlen):
                    new_sequence[index] = self._generate_seq( prob_matrix[parent_seq[index]] )
                    index+=1
                             
        # Attach final sequence to node
        current_node.seq = new_seq 