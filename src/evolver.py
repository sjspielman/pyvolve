#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
The module will evolve sequences along a phylogeny.
'''

import itertools
from copy import deepcopy
import numpy as np
from scipy import linalg
import random as rn
from model import *
from newick import *
from genetics import *
from partition import *
ZERO      = 1e-8
MOLECULES = Genetics()
        
        
class Site():
    '''
        Defines the Site class, which holds information for each evolved site.
        Currently the sequence, but this will be expanded.
    '''
    def __init__(self):
        self.int_seq      = None # integer sequence at a site
        #self.position     = None # location of site in full alignment size. <- eventually, this will be what gets shuffled.



class Evolver(object):
    ''' 
        This callable class evolves sequences along a phylogeny. By default, Evolver will evolve sequences and create several output files:
            1. simulated_alignment.fasta
                - The resulting simulated alignment
            2. site_rates.txt 
                - File providing rate information about each simulated column. Gives the partition and rate cateogy for each site in final simulated alignment.
                - Tab-delimited file with fields, Site_Index    Partition_Index     Rate_Category . All indexing is from *1*.
            3. site_rates_info.txt
                - File providing the true site-rate heterogeneity values (either the rate scaling factor or dN and dS) for each rate category.
                - Tab-delimited file with fields, Partition_Index    Model_Name    Rate_Category    Rate_Probability    Rate_Factor
          
          Note that file creation may be suppressed or files may be renamed using optional arguments given below.
    
    '''    
    def __init__(self, **kwargs):
        '''             
            Required keyword arguments include,
                1. **tree** is the phylogeny (parsed with the ``newick.read_tree`` function) along which sequences are evolved
                2. **partitions** is a list of Partition instances to evolve
    
            Optional keyword arguments include,
                1. **noisy_branch_lengths** is a boolean argument (True or False) for whether noise should be added to the branch lengths given in the newick tree. Default is False. \n If True, the branch lengths for each site along a branch will be sampled from a uniform distribution with center equal to the provided branch length. The range is, by default, 10% of the center. This 10% factor can be customized with the argument "noisy_branch_lengths_scale".
                2. **noisy_branch_lengths_n** is the number of noisy branch lengths. By default, when noisy_branch_lengths is True, 10 branch lengths are drawn per branch and randomly applied to sites. This value may be customized using this argument. Note: the argument "full" means a unique branch length at each site (this will likely be quite slow!!).
                3. **noisy_branch_lengths_scale** is a scaling factor to determine the range of the uniform distribution for drawing branch lengths with some noise. This option is only used if noisy_branch_lengths is True. Default 0.1.
        
        '''
        
                
        self.partitions = kwargs.get('partitions', None)
        self.full_tree  = kwargs.get('tree', Tree())
        self.bl_noise   = kwargs.get('noisy_branch_lengths', False)
        self.bl_noise_n = kwargs.get('noisy_branch_lengths_n', 10)
        self.bl_noise_scale = kwargs.get('noisy_branch_lengths_scale', 0.1)
        
        # ATTRIBUTE FOR THE sitewise_dnds_mutsel PROJECT
        self.select_root_type = kwargs.get('select_root_type', 'random') # other options are min, max to select the lowest prob and highest prob state, respectively, for the root sequence.
                
        # These dictionaries enable convenient post-processing of the simulated alignment. 
        self._leaf_sites = {} # Store final tip Site lists only
        self._evolved_sites = {} # Stores Site lists from all nodes, including internal and tips
        
        # Setup and sanity checks 
        self._root_seq_length = 0
        self._setup_partitions()
        self._set_code()
        self._bl_noise_sanity()


           

    def _bl_noise_sanity(self):
        ''' 
            Perform some sanity checks on noisy branch length preferences.
            Conditions for each keyword argument:
                1. noisy_branch_lengths must be boolean
                2. noisy_branch_lengths_n must be a positive integer <= than the full sequence length or the string "full"
                3. noisy_branch_lengths_scale must be a positive float
            Note that last 2 conditions are only checked if 1 is True
        '''
        try:
            self.bl_noise = bool(self.bl_noise)
        except:
            raise AssertionError("Argument noisy_branch_lengths must be True/False (or 1/0).")
        
        if self.bl_noise:
            if type(self.bl_noise_n) is int:
                assert( self.bl_noise_n > 0 and self.bl_noise_n <= self._root_seq_length ), "Value for noisy_branch_lengths_n should be either a postive integer (in range [1,partition size]) or the string 'full' (for each site has own branch length)." 
            else:
                assert( self.bl_noise_n == "full"), "Value for noisy_branch_lengths_n should be either a postive integer or the string 'full' (for each site has own branch length)." 
            if self.bl_noise_n == self._root_seq_length:
                self.bl_noise_n = "full"
        
            assert(self.bl_noise_scale > ZERO), "Value for noisy_branch_lengths_scale must be positive."
 
 

    def _setup_partitions(self):
        '''
            Setup and various sanity checks. 
        '''
        # If partitions is not a list but indeed a Partition, turn into a list. If not a partition, assert.
        if isinstance(self.partitions, Partition):
            self.partitions = [self.partitions]
        else:
            assert(type(self.partitions) is list), "\n\nYou must provide either a single Partition object or list of Partition objects to evolver."
            for p in self.partitions:
                assert(isinstance(p, Partition)), "\n\nYou must provide either a single Partition object or list of Partition objects to evolver." 
        
        # Assign root model flag to the full tree and determine length of root sequence
        for part in self.partitions:        
            self.full_tree.model_flag = part.root_model_name
            if part.branch_het():
                assert(self.full_tree.model_flag is not None), "\n\n Your root model name does not correspond to any of the Model()/CodonModel() objects' names provided to your Partition() objects."
            self._root_seq_length += sum( part.size )

        # Final check on size
        assert(self._root_seq_length > 0), "\n\nPartitions have no size!"
    
    
    
    
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
        
            
            
            
            
            
    def __call__(self, **kwargs):
        '''
            Simulate sequences, perform any necessary post-processing, and save sequences and/or other info to appropriate files.
 
            Optional keyword arguments:
                1. **seqfile** is a custom name for the output simulated alignment. Provide None or False to suppress file creation.
                2. **seqfmt**  is the format for seqfile (either fasta, nexus, phylip, phylip-relaxed, stockholm, etc. Anything that Biopython can accept!!) Default is FASTA.
                3. **ratefile** is a custom name for the "site_rates.txt" file. Provide None or False to suppress file creation.
                4. **infofile** is a custom name for the "site_rates_info.txt" file. Provide None or False to suppress file creation.
                5. **write_anc** is a boolean argument (True or False) for whether ancestral sequences should be output along with the tip sequences. Default is False.
       
            Examples:
                .. code-block:: python
                   
                   >>> # Evolve according to default settings
                   >>> evolve = Evolver(tree = my_tree, partitions = my_partition_list)()
        
                   >>> # Include ancestral sequences in output file
                   >>> evolve = Evolver(tree = my_tree, partitions = my_partition_list, write_anc = True)()

                   >>> # Custom sequence file name and format, and suppress rate information
                   >>> evolve = Evolver(tree = my_tree, partitions = my_partition_list, seqfile = "my_seqs.phy", seqfmt = "phylip", ratefile = None, infofile = None)()
      
        '''
        # Input arguments
        self.seqfile    = kwargs.get('seqfile', 'simulated_alignment.fasta')
        self.seqfmt     = kwargs.get('seqfmt', 'fasta').lower()
        self.write_anc  = kwargs.get('write_anc', False)
        self.ratefile   = kwargs.get('ratefile', 'site_rates.txt')
        self.infofile   = kwargs.get('infofile', 'site_rates_info.txt')


        # Simulate recursively
        self._sim_subtree(self.full_tree)

        # Shuffle sequences?
        self._shuffle_sites()

        # Convert Site dictionaries to sequence dictionaries
        self.leaf_seqs = self._convert_site_to_seq_dict(self._leaf_sites)
        self.evolved_seqs = self._convert_site_to_seq_dict(self._evolved_sites)

        # Save rate info, as needed       
        if self.ratefile:
            self._write_ratefile()
        if self.infofile:
            self._write_infofile()
        
        
        # Save sequences, as needed
        if self.seqfile:
            if self.write_anc:
                self._write_sequences(self.evolved_seqs)
            else:
                self._write_sequences(self.leaf_seqs)
    #########################################################################################                      
                        
                        
                        
                        
                        
    ######################## FUNCTIONS TO PROCESS SIMULATED SEQUENCES #######################              
    def _convert_site_to_seq_dict(self, seqdict):
        '''
            Return dictionary with key:value pairs of ID:sequence string from the self._leaf_sites or self._evolved_sites dictionaries.
        '''
        new_dict = {}
        for entry in seqdict:
            merged_entry = list( itertools.chain.from_iterable( seqdict[entry] ) )
            sequence = self._site_to_sequence( merged_entry )
            new_dict[entry] = sequence
        return new_dict


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
            In particular, we shuffle sequences in the self._evolved_sites dictionary, and then we copy over to the self._leaf_sites dictionary.            
        ''' 
        start = 0
        for part_index in range( len(self.partitions) ):            
            part = self.partitions[part_index]
            if part.shuffle:
                size = sum( part.size )
                part_pos = np.arange( size ) + start
                np.random.shuffle(part_pos)     
                for record in self._evolved_sites:
                    temp = []
                    for pp in part_pos:
                        temp.append( self._evolved_sites[record][part_index][pp] )
                    self._evolved_sites[record][part_index][start:start + size] = temp

        # Apply shuffling to self._leaf_sites
        for record in self._leaf_sites:
            self._leaf_sites[record] = self._evolved_sites[record]

               
                    


    def _write_sequences(self, seqdict):
        ''' 
            Write resulting sequences to a file in specified format.
        '''
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import generic_alphabet
        from Bio import SeqIO

        alignment = [] 
        for entry in seqdict:
            seq_object = SeqRecord( Seq( seqdict[entry] , generic_alphabet ), id = entry, description = "")
            alignment.append(seq_object)
        try:
            SeqIO.write(alignment, self.seqfile, self.seqfmt)
        except:
            raise AssertionError("\n Output file format is unknown. Consult with Biopython manual to see which I/O formats are accepted.\n NOTE: If you are attempting to save as phylip and are receiving this error, try seqfmt = 'phylip-relaxed' instead.")



    def _write_ratefile(self):
        '''
            Write ratefile, a tab-delimited file containing site-specific rate information. Considers leaf sequences only.
            Writes -   Site_Index    Partition_Index     Rate_Category
            All indexing is from *1*.
        '''
        refseq = self._leaf_sites.values()[0]
        with open(self.ratefile, 'w') as ratef:
            ratef.write("Site_Index\tPartition_Index\tRate_Category")
            site_index = 1
            for p in range(len(refseq)):
                part = refseq[p]
                for i in range(len(part)):
                    w = "\n" + str(site_index) + "\t" + str(p +  1) + "\t" + str(part[i].rate + 1)
                    ratef.write(w)
                    site_index += 1
        

    def _write_infofile(self):
        '''
            Write infofile, a tab-delimited file which maps ratefile to actual rate values. Considers leaf sequences only.
            Writes -   Partition_Index    Model_Name    Rate_Category    Rate_Probability    Rate_Factor
            All indexing is from *1*.
        '''
        with open(self.infofile, 'w') as infof:
            infof.write("Partition_Index\tModel_Name\tRate_Category\tRate_Probability\tRate_Factor")
            for p in range( len(self.partitions) ):
                part = self.partitions[p]  
                prob_list = part._root_model.rate_probs      
                        
                for m in part.models:
                    for r in range(len(prob_list)):
                        outstr = "\n" + str(p+1) + "\t" + str(m.name) + "\t" + str(r+1) + "\t" + str(round(prob_list[r], 4)) + "\t"
                        if m.codon_model():
                            infof.write(outstr + str(round(m.params['beta'][r],4)) + "," + str(round(m.params['alpha'][r],4)) )
                        else:
                            infof.write(outstr + str(round(m.rate_factors[r],4)) )
                  
                  
                                
    def _get_sequences(self, anc = False):
        '''
            Method to return the dictionary of simulated sequences.
            Default anc = False will return the leaf_sites dictionary.
            If anc == True, then will return the evolved_sites dictionary.
        '''
        if anc:
            return self.evolved_seqs
        else:
            return self.leaf_seqs
                            


        
        
    ######################### FUNCTIONS INVOLVED IN SEQUENCE EVOLUTION ############################
    def _draw_noisy_branch_lengths(self, center):
        '''
            Draw an array of length self._bl_noise_n, from a uniform distribution with range ( center*(1.-self.bl_noise_scale), center*(1.+self.bl_noise_scale) ).
            Note that if the minimum bound of this distribution is negative, it is changed to be ZERO.
            Randomly assign lengths to sites.
        '''
        
        # Determine range of branch lengths distribution
        min = center * (1. - self.bl_noise_scale)
        if min <= ZERO:
            min = ZERO
        max = center * (1. + self.bl_noise_scale)
        
        # Assign each site own branch length if "full" specified, or assign random mapping
        if self.bl_noise_n == "full":
            bls = np.random.uniform( min, max, size = self._root_seq_length)
            mapping = np.arange(0, self._root_seq_length)
        
        else:
            bls = np.random.uniform( min, max, size = self.bl_noise_n)
            mapping = np.random.randint(0, self.bl_noise_n, size = self._root_seq_length)
        
        return bls, mapping


    def _exponentiate_matrix(self, Q, t):
        '''
            Perform exponentiation on instantaneous matrix to produce produce transition matrix, P = exp(Qt)
            Assert that all rows sum to 1.
            Return P
        '''
        P = linalg.expm( np.multiply(Q, float(t) ) )
        assert( np.allclose( np.sum(P, axis = 1), np.ones(len(self._code))) ), "Rows in transition matrix do not each sum to 1."
        return P
        
        

    def _generate_transition_matrices(self, Q, t):
        '''
            Generate the transition matrix/ces for this branch, P = exp(Qt).
            Two options are possible:
                + Single transition matrix for entire branch (self._gamma_branch_lengths == False)
                + Multiple transition matrices along branch, to account for stochasticity in number of substitutions per site (as bl are expected values of this). (self._gamma_branch_lengths == True)
            Returns:
                + a numpy array of matrices (or the single matrix, as the case may be) to be used along the branch 
                + an array mapping each site to the matrix it will use. The array index is the site position and the value is the index of the matrix in the array of matrices.
        '''
                
        if self.bl_noise:
            if self.bl_noise_n == "full":
                matrices = np.zeros([self._root_seq_length, len(self._code), len(self._code)])
            else:
                matrices = np.zeros([self.bl_noise_n, len(self._code), len(self._code)])
            bls, mapping = self._draw_noisy_branch_lengths(t)
            for i in range(len(bls)):
                matrices[i] = self._exponentiate_matrix(Q, bls[i])
        
        else:
            matrices = np.array([self._exponentiate_matrix(Q, t)])
            mapping = np.zeros( self._root_seq_length )
            
        
        return matrices, mapping
            
                
    
    
    
    
    def _obtain_model(self, part, flag):
        '''
            Obtain the appropriate Model()/CodonModel() for evolution along a particular branch.
        '''
        my_model = None
        if part.branch_het():
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
            Generate a root sequence based on the stationary frequencies.
            Return a complete root sequence list of Site objects.
            
            NOTE: The select_root_type attribute is for the sitewise_dnds_mutsel project and was created on 4/30/15.
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
                    ########### SECTION EDITED FOR sitewise_dnds_mutsel PROJECT ############
                    if self.select_root_type == "min":
                        new_site.int_seq = np.argmin(root_model.params['state_freqs'])
                        
                    elif self.select_root_type == "max":
                        new_site.int_seq = np.argmax(root_model.params['state_freqs'])
                        
                    elif self.select_root_type == "random": 
                        new_site.int_seq = self._generate_prob_from_unif( root_model.params['state_freqs'] )
                    #########################################################################
                    part_root.append( new_site )
                    del new_site
            assert( len(part_root) == sum(part.size) ), "\n\nRoot sequence improperly generated for a partition, evolution cannot happen."
            root_sequence.append(part_root)
        return root_sequence

        
        
    def _sim_subtree(self, current_node, parent_node = None):
        ''' 
            Function to traverse a Tree object recursively and simulate sequences.
            Required positional arguments include,
                1. **current_node** is the node (either internal node or leaf) TO WHICH we evolving
                2. **parent_node** is the node we are evolving FROM. Default of None is only called when the root sequence is not yet made.
        '''
        
        # We are at the base and must generate root sequence
        if (parent_node is None):
            current_node.seq = self._generate_root_seq() # the .seq attribute is actually a list of Site() objects.
            #self._evolved_sites['root'] = current_node.seq
        else:
            current_node.seq = self._evolve_branch(current_node, parent_node) 
        self._evolved_sites[current_node.name] = current_node.seq

            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                self._sim_subtree(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self._leaf_sites[current_node.name] = current_node.seq

        
            
            
    def _check_parent_branch(self, parent_node, current_node):
        ''' 
            Function ensures that, for a given node we'd like to evolve to, an appropriate branch length exists. 
            If the branch length is acceptable, an evolutionary model is then assigned to the node if one has yet been assigned.
            
            Required positional arguments include,
                1. **parent_node** is node FROM which we evolve
                2. **current_node** is the node (either internal node or leaf) TO WHICH we evolve
        '''
        assert (parent_node.seq != None), "\n\nThere is no parent sequence from which to evolve!"
        assert (current_node.branch_length >= 0.), "\n\n Your tree has a negative branch length. I'm going to quit now."
        if current_node.model_flag is None:
            current_node.model_flag = parent_node.model_flag
            
            
            
            
            
    def _evolve_branch(self, current_node, parent_node):
        ''' 
            Function to evolve a given sequence during tree traversal.
            
            Required positional arguments include, 
                1. **current_node** is the node (either internal node or leaf) we are evolving TO
                2. **parent_node** is the node we are evolving FROM.
        '''

        # Ensure parent sequence exists and branch length is acceptable. Return the model flag to use here.
        self._check_parent_branch(parent_node, current_node)
 
        # Evolve only if branch length is greater than 0 (1e-8). 
        if current_node.branch_length <= ZERO:
            new_seq = deepcopy(parent_node.seq)
        
        else:
            new_seq = []            
            
            for p in range( len(self.partitions) ):
                # Obtain current model for this partition at this branch
                part = self.partitions[p]
                current_model = self._obtain_model(part, current_node.model_flag)
                index = 0
                part_new_seq = []  # will temporarily store this partition's new sequence
                
                for i in range( current_model.num_classes() ):
                    # Grab instantaneous rate matrix, which is done differently depending if codon (dN/dS) model or not. This is the rate het in the partition.
                    inst_matrix = None
                    if part.codon_model():
                        inst_matrix = current_model.matrices[i]
                    else:
                        inst_matrix = current_model.matrix * current_model.rate_factors[i] # note that rate_factors = [1.] if no site heterogeneity, so matrix unchanged
                    assert( inst_matrix is not None ), "\n\nCouldn't retrieve instantaneous rate matrix."
                    
                    # Generate transition matrix/ces. Matrices are all matrices for this branch, and mappings maps sites (position) to matrix (value)
                    matrices, mappings = self._generate_transition_matrices(inst_matrix, float(current_node.branch_length))
                
                    # Evolve branch
                    part_parent_seq = parent_node.seq[p][index : index + part.size[i]]
                    for j in range( part.size[i] ):
                        this_matrix = matrices[ mappings[index] ] # determine which matrix to use
                        new_site = deepcopy( part_parent_seq[j] )
                        new_site.int_seq = self._generate_prob_from_unif( this_matrix[ new_site.int_seq ] )
                        part_new_seq.append( new_site )
                        index += 1
                new_seq.append( part_new_seq )
        return new_seq

        
        
        
         
