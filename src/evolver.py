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
from misc import ZERO, Genetics, Model, Partition
MOLECULES = Genetics()
        
class Evolver(object):
    ''' 
        Class to evolve sequences along a phylogeny. 
        Currently supported:
            1. Site heterogeneity (via partitions)
            2. Branch heterogeneity (via model flags in phylogeny, similar to approach used by Indelible)
        
        Coming soon:
            1. Indels
    '''
    
    def __init__(self, partitions):
        # The first argument should be a list of the partitions. *required*
        # The second argument should be the name of the evolutionary model at the root of the tree. This argument MUST be provided when there is branch heterogeneity, but if the process is time-homogeneous then it is not needed.
        
        self._partitions = partitions
        self._setup_partitions()
        self.alndict = {} # Will store final alignment (TIPS ONLY)


    def _setup_partitions(self):
        '''
            Setup and some sanity checks on partitions and root_model. Also determine  full_seq_length.
        '''
        if isinstance(self._partitions, Partition):
            self._partitions = [self._partitions]
        else:
            assert(type(self._partitions) is list), "\n\nMust provide either a single Partition object or list of Partition objects to evolver."
        self._full_seq_length = 0
        
        for part in self._partitions:
            
            # Full sequence length [THIS WILL BE REMOVED WHEN INDELS ARE INCORPORATED]
            if type(part.size) is int:
                part.size = [part.size]
            self._full_seq_length += sum( part.size )
            
            # Branch *homogeneity*
            if isinstance(part.model, Model):
                part.root_model = None
                dim = part.model.params['state_freqs'].shape[0]
                 
            # Branch *heterogeneity*
            elif type(part.model) is list:
                dim = part.model[0].params['state_freqs'].shape[0]
                found_root = False
                for m in part.model:
                    if m.name == part.root_model:
                        found_root = True
                        break
                assert(found_root is True), "\n\n Your root_model does not correspond to any of the Model() objects provided to your Partition() objects."
                
        assert(self._full_seq_length > 0), "Partitions have no size!" 
        self._set_code(dim)



    def _set_code(self, dim):
        ''' 
            Assign genetic code.
        '''    
        if dim == 4:
            self._code = MOLECULES.nucleotides
        elif dim == 20:
            self._code = MOLECULES.amino_acids
        elif dim == 61:
            self._code = MOLECULES.codons
        else:
            raise AssertionError("This should never be reached.")
            
                        
                        
                        
                        
    def _sequence_to_integer(self, entry):
        ''' 
            Convert a dna/protein character to its appropriate integer (index in self._code).
            Argument *entry* is the character to convert.
        '''
        return self._code.index(entry)
    
    
    
    def _integer_to_sequence(self, index):
        '''
            Convert an integer (index in self._code) to its appropriate dna/protein character.
            Argument *index* is the integer to convert.
        '''
        return self._code[index]
 
 
    def _intseq_to_string(self, intseq):
        ''' 
            Convert a full sequence coded as integers (indices in self._code) and convert to a dna/protein string.
            Argument *intseq* is the sequence to convert.
        '''
        stringseq = ''
        for i in intseq:
            stringseq += self._integer_to_sequence(i)
        return stringseq   



    def write_sequences(self, **kwargs):
        ''' 
            Write resulting sequence alignment (self.alndict) to a file in specified format.
            Arguments:
                1. "outfile" is the name of the file for saving the alignment
                2. "format" is the alignment output file format (either fasta, nexus, phylip, phylip-relaxed, stockholm, etc. Anything that Biopython can accept!!) If not provided, will output in fasta format.
        '''
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import generic_alphabet
        from Bio import SeqIO

        format   = kwargs.get("format", "fasta").lower()
        outfile  = kwargs.get("outfile", "simulated_alignment.txt") 
        alignment = [] 
        for entry in self.alndict:
            seq_object = SeqRecord( Seq( self._intseq_to_string( self.alndict[entry] ) , generic_alphabet ), id = entry, description = "")
            alignment.append(seq_object)
        try:
            SeqIO.write(alignment, outfile, format)
        except:
            raise AssertionError("\n Output file format is unknown. Consult with Biopython manual to see which I/O formats are accepted.")


####################################
    def _shuffle_sites(self):
        ''' 
            Shuffle evolved sequences. Can either shuffle within each partition or shuffle the entire alignment. Column integrity is maintained.
        '''
####################################



    def _generate_prob_from_unif(self, prob_array):
        ''' 
            Sample a sequence (nuc,aa,or codon), and return an integer for the sequence chosen from a uniform distribution.
            Arugment *prob_array* is any list and/or numpy array of probabilities which sum to 1.
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
        
        root_sequence = np.empty(self._full_seq_length, dtype=int)
        index = 0
        
        for part in self._partitions:
        
            # Obtain root frequencies
            if isinstance(part.model, Model):  #if part.root_model is None:
                freqs = part.model.params['state_freqs']
            else:
                for m in part.model:
                    if m.name == part.root_model:
                        freqs = m.params['state_freqs']
                        break
    
            
            # Simulate root          
            for j in range( sum(part.size) ):
                root_sequence[index] = self._generate_prob_from_unif(freqs)
                index += 1
                
        return root_sequence, part.root_model

        
        
    def simulate(self, current_node, parent_node = None):
        ''' 
            Function to traverse a Tree object recursively and simulate sequences.
            Arguments:
                1. *current_node* is the node (either internal node or leaf) TO WHICH we evolving
                2. *parent_node* is the node we are evolving FROM. Default of None is only called when the root sequence is not yet made.
        '''
        
        # We are at the base and must generate root sequence
        if (parent_node is None):
            current_node.seq, current_node.model_flag = self._generate_root_seq() 
        else:
            self.evolve_branch(current_node, parent_node) 
            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                self.simulate(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[current_node.name]=current_node.seq
            
            
            
            
            
    def _check_parent_branch(self, parent_node, current_node):
        ''' 
            Function ensures that, for a given node we'd like to evolve to, an appropriate branch length exists. 
            If the branch length is acceptable, an evolutionary model is then assigned to the node if one has yet been assigned. This will typically be the case of 
            
            Arguments:
                1. *parent_node* is node FROM which we evolve
                2. *current_node* is the node (either internal node or leaf) TO WHICH we evolve
        '''
        assert (parent_node.seq != None), "\n\nThere is no parent sequence from which to evolve!"
        assert (current_node.branch_length > 0), "\n\n Your tree has a negative branch length. I'm going to quit now."
        if current_node.model_flag is None:
            current_node.model_flag = parent_node.model_flag

            
            
    def evolve_branch(self, current_node, parent_node):
        ''' 
            Function to evolve a given sequence during tree traversal.
            Arguments:
                1. *current_node* is the node (either internal node or leaf) we are evolving TO
                2. *parent_node* is the node we are evolving FROM.
        '''
    
        # Ensure parent sequence exists and branch length is acceptable. Return the model flag to use here.
        self._check_parent_branch(parent_node, current_node)
 
        # Evolve only if branch length is greater than 0 (1e-8). 
        if current_node.branch_length <= ZERO:
            new_seq = parent_seq
        
        else:
            new_seq = np.empty(self._full_seq_length, dtype=int)
            current_model = None
            index = 0
            for part in self._partitions:
                
                # Obtain current model, inst_matrix
                if isinstance(part.model, Model):
                    current_model = part.model
                else:
                    for m in part.model:
                        if m.name == current_node.model_flag:
                            current_model = m
                            break
                assert( current_model is not None ), "\n\nCould not retrieve model for partition in evolve_branch."

                # Incorporate rate heterogeneity (for nuc, amino acid, or ECM models) if specified. If homogeneous, part.size will be len=1 anyways.
                for x in range(len(part.size)):
                
                    # Generate probability matrix and assert correct
                    inst_matrix = current_model.matrix * current_model.rates[x]
                    prob_matrix = linalg.expm( np.multiply(inst_matrix, float(current_node.branch_length) ) )
                    assert( np.allclose( np.sum(prob_matrix, axis = 1), np.ones(len(self._code))) ), "Rows in transition matrix do not each sum to 1."
                
                    # Evolve branch
                    for j in range( part.size[x] ):
                        new_seq[index] = self._generate_prob_from_unif( prob_matrix[parent_node.seq[index]] )
                        index+=1
                             
        # Attach final sequence to node
        current_node.seq = new_seq 
