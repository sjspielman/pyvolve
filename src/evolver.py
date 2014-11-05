#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Evolve sequences along a phylogeny.

TODO: 
1. Provide output file for gamma rate categories. In cases of branch het, the model for each set of categories should be indicated. 
'''

import itertools
from copy import deepcopy
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
            
            Optional arguments:    
                1. *seqfile* is an output file for saving final simulated sequences
                2. *seqfmt* is the format for seqfile (either fasta, nexus, phylip, phylip-relaxed, stockholm, etc. Anything that Biopython can accept!!) Default is FASTA.
                3. *write_anc* is a bool for whether ancestral sequences should be output. If so, they are output with the tip sequences in seqfile. Default is False.
                4. *ratefile* is an output file for saving rate information about each simulated column. Gives the partition, rate cateogy for each site in final simulated alignment.
                5. *infofile* is an output file for saving the actual site-rate heterogeneity values (either the rate scaling factor or dN and dS).
        
                Note: for all output files, users may specify False (e.g. seqfile = False) to suppress output. Otherwise, files are automatically output with default names, given below.
        '''
        
                
        self.partitions = kwargs.get('partitions', None)
        self.full_tree  = kwargs.get('tree', Tree())
        self.seqfile    = kwargs.get('seqfile', 'simulated_alignment.fasta')
        self.seqfmt     = kwargs.get('seqfmt', 'fasta').lower()
        self.write_anc  = kwargs.get('write_anc', False)
        self.ratefile   = kwargs.get('ratefile', 'site_rates.txt')
        self.infofile   = kwargs.get('infofile', 'site_rates_info.txt')
        
        # These dictionaries enable convenient post-processing of the simulated alignment. Otherwise we'd have to always loop over full tree, which would be very slow.
        self.leaf_seqs = {} # Store final tip sequences only
        self.evolved_seqs = {} # Stores sequences from all nodes, including internal and tips
        
        # Setup and sanity checks 
        self._setup_partitions()
        self._set_code()

            


    def _setup_partitions(self):
        '''
            Setup and various sanity checks. 
        '''
        if isinstance(self.partitions, Partition):
            self.partitions = [self.partitions]
        else:
            assert(type(self.partitions) is list), "\n\nYou must provide either a single Partition object or list of Partition objects to evolver."
            for p in self.partitions:
                assert(isinstance(p, Partition)), "\n\nYou must provide either a single Partition object or list of Partition objects to evolver." 
        
        self._root_seq_length = 0
        for part in self.partitions:
        
            ################ Set up branch heterogeneity, if specified ################
            # Yes branch heterogeneity -> sanity check the (hopefully) specified root_model, and yell if not assigned or assigned incorrectly.
            if part.branch_het():
                for m in part.models:
                    if m.name == part.root_model:
                        self.full_tree.model_flag = part.root_model
                assert(self.full_tree.model_flag is not None), "\n\n Your root_model does not correspond to any of the Model()/CodonModel() objects provided to your Partition() objects."
            else:
                part.models = [part.models]
                part.root_model = None
                        
            ################ Set up rate heterogeneity, if specified ################
            # Site-rate classes will not change, even with branch heterogeneity, so can simply pick first model.
            m = part.models[0]
            full = part.size
            # No rate heterogeneity. Turn part.size into list of length 1. No need to shuffle.
            if m.num_classes() == 1:
                part.size = [part.size]
            
            # Yes rate heterogeneity. 
            else:
                # Divvy up part.size into rate het chunks, and set shuffle to True
                part.shuffle = True
                remaining = part.size
                part.size = []
                for i in range(len(m.rates) - 1): # don't fill in last one yet since rounding issues will occur.
                    section = int( m.probs[i] * full )
                    part.size.append( section )
                    remaining -= section
                part.size.append(remaining)  
                
                # Ensure that all the models have properly normalized rates, or fix them accordingly
                for model in part.models:
                    model.probs, model.rates = self._setup_rates(model.probs, model.rates)
                  
            assert( sum(part.size) ==  full ), "\n\nImproperly divvied up rate heterogeneity."
            self._root_seq_length += full


        ################ Final check on size ################      
        assert(self._root_seq_length > 0), "\n\nPartitions have no size!"


    def _setup_rates(self, probs, rates):
        ''' 
            Sanity check rate categories for site-heterogeneity.
            Arguments *probs* and *rates* are the rate category probabilities and scalars, respectively, for the model we are checking
        '''
        probs = np.array(probs)
        rates = np.array(rates)
        if abs( 1. - np.sum(probs)) > ZERO:
            probs /= np.sum(probs)
        if abs( 1. - np.sum(probs * rates)) > ZERO:
            rates /= np.sum(rates * probs)
        return rates, probs
    
    
    
    def _set_code(self):
        ''' 
            Assign genetic code.
        '''    
        dim = self.partitions[0].models[0].params['state_freqs'].shape[0] 
        if dim == 4:
            self._code = MOLECULES.nucleotides
        elif dim == 20:
            self._code = MOLECULES.amino_acids
        elif dim == 61:
            self._code = MOLECULES.codons
        else:
            raise AssertionError("This should never be reached.")
        
            
            
            
            
            
    def __call__(self):
        '''
            Run evolver! Perform simulation, any necessary post-processing, and save sequences and/or other info to appropriate files.
        '''

        # Simulate recursively
        self.sim_subtree(self.full_tree)

        # Shuffle sequences?
        self._shuffle_sites()

        # Save rate info        
        if self.ratefile:
            self.write_ratefile()
        if self.infofile:
            self.write_infofile()
        
        
        # Save sequences
        if self.seqfile:
            if self.write_anc:
                self.write_sequences(self.seqfile, self.seqfmt, self.evolved_seqs)
            else:
                self.write_sequences(self.seqfile, self.seqfmt, self.leaf_seqs)
    #########################################################################################                      
                        
                        
                        
                        
                        
    ######################## FUNCTIONS TO PROCESS SIMULATED SEQUENCES #######################              
    def _site_to_sequence(self, site):
        '''
            Convert a single Site() object or list of Site() objects into a sequence string.
            Argument *site* is either a Site() object or a list of them.
        '''
        if type(site) is list:
            sequence = ''
            for s in site:
                sequence += self._code[s.int_seq]
        else:
            sequence = self._code[site.int_seq]
        return sequence




    def _shuffle_sites(self):
        ''' 
            Shuffle evolved sequences within partitions, if specified.
            In particular, we shuffle sequences in the self.evolved_seqs dictionary, and then we copy over to the self.leaf_seqs dictionary.
            
            LATER: only shuffle sites which have Site.origin == "root". Inserted sites don't need to be shuffled.
            
        ''' 
        start = 0
        for part_index in range( len(self.partitions) ):            
            part = self.partitions[part_index]
            if part.shuffle:
                size = sum( part.size )
                part_pos = np.arange( size ) + start
                np.random.shuffle(part_pos)     
                for record in self.evolved_seqs:
                    temp = []
                    for pp in part_pos:
                        temp.append( self.evolved_seqs[record][part_index][pp] )
                    self.evolved_seqs[record][part_index][start:start + size] = temp

        # Apply shuffling to self.leaf_seqs
        for record in self.leaf_seqs:
            self.leaf_seqs[record] = self.evolved_seqs[record]

               
                    


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
            merged_entry = list( itertools.chain.from_iterable( seqdict[entry] ) )
            sequence = self._site_to_sequence( merged_entry )
            seq_object = SeqRecord( Seq( sequence , generic_alphabet ), id = entry, description = "")
            alignment.append(seq_object)
        try:
            SeqIO.write(alignment, outfile, fmt)
        except:
            raise AssertionError("\n Output file format is unknown. Consult with Biopython manual to see which I/O formats are accepted.")



    def write_ratefile(self):
        '''
            Write ratefile, a tab-delimited file containing site-specific rate information. Considers leaf sequences only.
            Writes -   Site_Index    Partition_Index     Rate_Category
            All indexing is from *1*.
        '''
        refseq = self.leaf_seqs.values()[0]
        with open(self.ratefile, 'w') as ratef:
            ratef.write("Site_Index\tPartition_Index\tRate_Category")
            site_index = 1
            for p in range(len(refseq)):
                part = refseq[p]
                for i in range(len(part)):
                    w = "\n" + str(site_index) + "\t" + str(p +  1) + "\t" + str(part[i].rate + 1)
                    ratef.write(w)
                    site_index += 1
        

    def write_infofile(self):
        '''
            infofile.
        '''
        with open(self.infofile, 'w') as infof:
            infof.write("Partition\tModel_Name\tRate_Category\tRate_Probability\tRate_Factor")
            
            for p in range( len(self.partitions) ):
                part = self.partitions[p]                
                for m in part.models:
                    if m.num_classes() == 1:
                        infof.write("\n" + str(p+1) + "\t" + str(m.name) + "\t1\t1")
                    if m.num_classes() > 1:
                        for r in range(len(m.rates)):
                            outstr = "\n" + str(p+1) + "\t" + str(m.name) + "\t" + str(r+1) + "\t" + str(round(m.probs[r],4)) + "\t"
                            if isinstance(m, CodonModel):
                                infof.write(outstr + str(round(m.rates[0],4)) + "\t" + str(round(m.rates[1],4)) )
                            else:
                                infof.write(outstr + str(round(m.rates[r],4)) )
                                
                            
                        
        
        
        
        
    ######################### FUNCTIONS INVOLVED IN SEQUENCE EVOLUTION ############################
    def _obtain_model(self, part, flag):
        '''
            Obtain the appropriate Model()/CodonModel() for evolution along a particular branch.
        '''
        my_model = None
        if part.branch_het():
            my_model = None
            for m in part.models:
                if m.name == flag:
                    my_model = m
                    break
        else:
            my_model = part.models[0]
        assert( my_model is not None ), "\n\nCould not retrieve model a particular branch's evolution."
        return my_model


    
    
    
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
            Generate a root sequence based on the stationary frequencies or using a provided sequence
            Return a complete root sequence list of Site objects.
        '''
        
        root_sequence = [] # This will contain a list for each partition's sequence (which is itself a list of Site() objects)

        for part in self.partitions:
            
            # Grab model info for this partition to get frequency vector for root simulation
            root_model = self._obtain_model(part, self.full_tree.model_flag)

            # Generate root_sequence and assign the Site a rate class
            part_root = []
            for i in range( root_model.num_classes() ):
                for j in range( part.size[i] ):
                    new_site = Site()
                    new_site.rate = i
                    new_site.int_seq = self._generate_prob_from_unif( root_model.params['state_freqs'] )
                    part_root.append( new_site )
            assert( len(part_root) == sum(part.size) ), "\n\nRoot sequence improperly generated for a partition, evolution cannot happen."
            root_sequence.append(part_root)
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
            current_node.seq = self._generate_root_seq() # the .seq attribute is actually a list of Site() objects.
            self.evolved_seqs['root'] = current_node.seq
        else:
            current_node.seq = self.evolve_branch(current_node, parent_node) 
            self.evolved_seqs[current_node.name] = current_node.seq

            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                self.sim_subtree(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.leaf_seqs[current_node.name] = current_node.seq

        
            
            
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
            new_seq = deepcopy(parent_node.seq)
        
        else:
            new_seq = []
            for p in range( len(self.partitions) ):

                # Obtain current model
                part = self.partitions[p]
                current_model = self._obtain_model(part, current_node.model_flag)
                
                index = 0
                temp_new_seq = [] 
                for i in range( current_model.num_classes() ):
                    # Grab instantaneous rate matrix, which is done differently depending if codon (dN/dS) model or not. This is the rate het in the partition.
                    inst_matrix = None
                    if part.codon_model():
                        inst_matrix = current_model.matrices[i]
                    else:
                        inst_matrix = current_model.matrix * current_model.rates[i] # note that rates = [1.] if no site heterogeneity, so matrix unchanged
                    assert( inst_matrix is not None ), "\n\nCouldn't retrieve instantaneous rate matrix."
                    
                    # Generate transition matrix and assert correct
                    prob_matrix = linalg.expm( np.multiply(inst_matrix, float(current_node.branch_length) ) )
                    assert( np.allclose( np.sum(prob_matrix, axis = 1), np.ones(len(self._code))) ), "Rows in transition matrix do not each sum to 1."
                
                    # Evolve branch
                    part_parent_seq = parent_node.seq[p][index : index + part.size[i]]
                    for j in range( part.size[i] ):
                        new_site = deepcopy( part_parent_seq[j] )
                        new_site.int_seq = self._generate_prob_from_unif( prob_matrix[ new_site.int_seq ] )
                        temp_new_seq.append( new_site )
                        index += 1
                new_seq.append( temp_new_seq )
  
        return new_seq

        
        
        
         
