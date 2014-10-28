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
from misc import *
MOLECULES = Genetics()
        
class Evolver(object):
    ''' 
        Class to evolve sequences along a phylogeny. 
        Currently supported:
            1. Site heterogeneity, among and within partitions
            2. Branch heterogeneity (via model flags in phylogeny, similar to approach used by Indelible)
        
        Coming soon:
            1. Indels
    '''
    
    def __init__(self, **kwargs):
        '''             
            Required arguments:
                1. *tree* is the phylogeny along which we will evolve
                2. *partitions* is a list of Partition() objects to evolve 
            
            Required arguments for branch heterogeneity [homogenous models need not apply]:
                1. *root_model* is the name of the model at root of the tree. This argument is unneeded and entirely useless if the process is time-homogenous, although it won't hurt to provide something... it's just very silly.
            
            Optional arguments:    
                1. *seqfile* is an output file for saving final simulated sequences
                2. *seqfmt* is the format for seqfile (either fasta, nexus, phylip, phylip-relaxed, stockholm, etc. Anything that Biopython can accept!!) Default is FASTA.
                3. *write_anc* is a bool for whether ancestral sequences should be output. If so, they are output with the tip sequences in seqfile. Default is False.
            
            TODO arguments:
                1. *ratefile* is an output file for saving rate information about each simulated column. For codon sequences, saves dN and dS. For nucleotide and amino acid sequences, saves heterogeneity info (if applicable)
        '''
                
        self.partitions = kwargs.get('partitions', None)
        self.full_tree  = kwargs.get('tree', Tree())
        self.root_seq   = kwargs.get('root_seq', None)
        self.seqfile    = kwargs.get('seqfile', None)
        self.seqfmt     = kwargs.get('seqfmt', 'fasta').lower()
        self.write_anc  = kwargs.get('write_anc', False)
        
        self.leaf_seqs = {} # Store final tip sequences only
        self.evolved_seqs = {} # Stores sequences from all nodes, including internal and tips
        
        # Setup partitions, with some sanity checking, before evolution begins. 
        self._setup_partitions()
        

    def _setup_partitions(self):
        '''
            Setup and some sanity checks on partitions and root_model. Also determine  full_seq_length.
        '''
        if isinstance(self.partitions, Partition):
            self.partitions = [self.partitions]
        else:
            assert(type(self.partitions) is list), "\n\nMust provide either a single Partition object or list of Partition objects to evolver."
        self._full_seq_length = 0
        
        for part in self.partitions:
            
            # Full sequence length and root_seq sanity. [THIS CHUNK WILL BE REMOVED WHEN INDELS ARE INCORPORATED]
            #####
            if type(part.size) is int:
                part.size = [part.size]
            if part.root_seq:
                assert( len(part.root_seq) == sum( part.size ) ), "\n\nProvided root sequence is wrong length. Since we don't have indels yet, a given partition's root sequence *must* be the same size as the partition!"
            self._full_seq_length += sum( part.size )
            #####
            
            # Branch *homogeneity*
            if isinstance(part.model, Model):
                part.root_model = None
                dim = part.model.params['state_freqs'].shape[0]
                 
            # Branch *heterogeneity*. Ensure that a appropriate model flags have been specified.
            elif type(part.model) is list:
                dim = part.model[0].params['state_freqs'].shape[0]
                for m in part.model:
                    if m.name == part.root_model:
                        self.full_tree.model_flag = part.root_model
                        break
                assert(self.full_tree.model_flag is not None), "\n\n Your root_model does not correspond to any of the Model() objects provided to your Partition() objects."
                
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
        
            
      
    ####################### CRUX SIMULATION AND SAVING FUNCTION #############################
    def simulate(self):
        '''
            Crux function for sequence simulation, post-processing, and file saving.
        '''

        # Simulate recursively
        self.sim_subtree(self.full_tree)
        
        self.write_sequences("preshuf.fasta", "fasta", self.leaf_seqs)
        
        # Shuffle sequences?
        #if self._shuffle:
        shuffle_list = self._shuffle_sites() # returns a list of how things were shuffled, to deal with rates and such.
        
        self.write_sequences("postshuf.fasta", "fasta", self.leaf_seqs)
        
        assert 1==5
        # Save sequences?
        if self.seqfile is not None:
            if self.write_anc:
                self.write_sequences(self.seqfile, self.seqfmt, self.evolved_seqs)
            else:
                self.write_sequences(self.seqfile, self.seqfmt, self.leaf_seqs)
    #########################################################################################                      
                        
                        
                        
                        
                        
    ######################## FUNCTIONS TO PROCESS SIMULATED SEQUENCES #######################              
  
    def _sequence_to_integer(self, entry):
        ''' 
            Convert a dna/protein character to its appropriate integer (index in self._code).
            Argument *entry* is the character to convert.
        '''
        return self._code.index(entry)
        
        
    def _sequence_to_intseq(self, seqstring):
        '''
            Convert a full sequence string into numpy array of integers.
            Argument *seqstring* is the sequence to convert.
        '''
        intseq = np.empty(len(seqstring))
        for i in range(len(seqstring)):
            intseq[i] = self._sequence_to_integer(seqstring[i])
        return intseq
        
        
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



    def _shuffle_sites(self):
        ''' 
            10/26/14 DEBUGGED BUT NEEDS TESTING!!
            Shuffle evolved sequences. Can either shuffle within each partition or shuffle the entire alignment. Column integrity is maintained.
            In particular, we shuffle sequences in the self.evolved_seqs dictionary, and then we copy over to the self.leaf_seqs dictionary.
            
        '''
        group_shuffle = [] # Indices from mul. partitions which should be shuffled as a single unit (partitions whose .part==2)
        shuffle_map   = {} # Key is final site index, value is FORMER site index (where that column came from)
        
               
        # Shuffle within individual partitions, as needed
        start = 0
        for part in self.partitions:
            size = sum(part.size)
            
            # No shuffle
            if part.shuffle == 0:
                # Map 1:1
                for entry in ( np.arange(size) + start ):
                    shuffle_map[entry] = entry
            
            # Shuffle partition
            elif part.shuffle == 1:
                part_pos = np.arange( size ) + start
                np.random.shuffle(part_pos)       
                for record in self.evolved_seqs:
                    self.evolved_seqs[record][start:start + size] = self.evolved_seqs[record][part_pos]
                # Map the shuffled sites
                for i in range(start, start + size):
                    shuffle_map[i] = part_pos[ i - start ]
         
            # Add partition to global shuffle list
            elif part.shuffle == 2:
                group_shuffle.extend( np.arange( size ) + start )
            start += size
    
        # Shuffle remaining partition groups together, as needed
        if len(group_shuffle) > 0:
            orig = np.array( group_shuffle )
            np.random.shuffle( group_shuffle )
            orig = np.sort(orig)
            
            for record in self.evolved_seqs:
                self.evolved_seqs[record][ orig ] = self.evolved_seqs[record][ group_shuffle ]
            for i in range(len( orig )):
                shuffle_map[ orig[i] ] = group_shuffle[i] 
                     
        # Apply shuffling to self.leaf_seqs
        for record in self.leaf_seqs:
            self.leaf_seqs[record] = self.evolved_seqs[record]
        
        print shuffle_map
        return shuffle_map


               
                    


    def write_sequences(self, outfile, fmt, seqdict):
        ''' 
            Write resulting sequences (seqdict - this is either just the tips or all nodes) to a file in specified format.
            Arguments:
                1. "outfile" is the name of the file for saving the alignment
                2. "fmt" is the alignment output file format (either fasta, nexus, phylip, phylip-relaxed, stockholm, etc. Anything that Biopython can accept!!) If not provided, will output in fasta format.
        '''
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import generic_alphabet
        from Bio import SeqIO

        alignment = [] 
        for entry in seqdict:
            seq_object = SeqRecord( Seq( self._intseq_to_string( seqdict[entry] ) , generic_alphabet ), id = entry, description = "")
            alignment.append(seq_object)
        try:
            SeqIO.write(alignment, outfile, fmt)
        except:
            raise AssertionError("\n Output file format is unknown. Consult with Biopython manual to see which I/O formats are accepted.")




    ######################### FUNCTIONS INVOLVED IN SEQUENCE EVOLUTION ############################
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
        
        root_sequence = np.empty(self._full_seq_length, dtype = int)
        
        index = 0
        for part in self.partitions:
        
            # Add a specified root sequence to the array
            if part.root_seq:
                root_sequence[index : sum(part.size) ] = self._sequence_to_intseq(part.root_seq)
                index += sum(part.size)
            
            # Generate a root sequence
            else:
                
                # Determine root model for selecting freq vector. The root model name has been stored in self.full_tree.model_flag
                if isinstance(part.model, Model):
                    freqs = part.model.params['state_freqs']
                else:
                    for m in part.model:
                        if m.name == self.full_tree.model_flag:
                            freqs = m.params['state_freqs']
                            break
    
                # Generate
                for j in range( sum(part.size) ):
                    root_sequence[index] = self._generate_prob_from_unif(freqs)
                    index += 1
            
        assert( np.all(root_sequence) >= 0 and np.all(root_sequence) < len(self._code) ), "\n\n Root sequence improperly generated, evolution cannot proceed."
        return root_sequence

        
        
    def sim_subtree(self, current_node, parent_node = None):
        ''' 
            Function to traverse a Tree object recursively and simulate sequences.
            Arguments:
                1. *current_node* is the node (either internal node or leaf) TO WHICH we evolving
                2. *parent_node* is the node we are evolving FROM. Default of None is only called when the root sequence is not yet made.
        '''
        
        # We are at the base and must generate root sequence
        if (parent_node is None):
            current_node.seq = self._generate_root_seq()
            self.evolved_seqs['root'] = current_node.seq
        else:
            self.evolve_branch(current_node, parent_node) 
            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                self.sim_subtree(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.leaf_seqs[current_node.name]=current_node.seq
            
            
            
            
            
    def _check_parent_branch(self, parent_node, current_node):
        ''' 
            Function ensures that, for a given node we'd like to evolve to, an appropriate branch length exists. 
            If the branch length is acceptable, an evolutionary model is then assigned to the node if one has yet been assigned.
            
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
            for part in self.partitions:
                
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
                             
        # Attach final sequence to node and save to self.evolved_dict
        current_node.seq = new_seq
        self.evolved_seqs[ current_node.name ] = current_node.seq
        
        
        
        
        
         
