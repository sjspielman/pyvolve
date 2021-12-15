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
from .model import *
from .newick import *
from .genetics import *
from .partition import *
ZERO      = 1e-8
MOLECULES = Genetics()


class Site():
    '''
        Defines the Site class, which holds information for each evolved site.
    '''
    def __init__(self):
        self.int_seq      = None # integer sequence at a site
        self.rate         = None # site rate category


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
                2. **partitions** (or **partition**) is a list of Partition instances to evolve.
        '''

        self.partitions = kwargs.get('partitions', None)
        if self.partitions is None:
            self.partitions = kwargs.get('partition', None)
        self.full_tree  = kwargs.get('tree', Node())



        # These dictionaries enable convenient post-processing of the simulated alignment.
        self._leaf_sites = {} # Store final tip Site lists only
        self._evolved_sites = {} # Stores Site lists from all nodes, including internal and tips

        # Setup and sanity checks
        self._root_seq_length = 0
        self._setup_partitions()
        self._code = self.partitions[0]._root_model.code

        ######### In-house projects #######
        # start root with certain fitness #
        self.select_root_type = kwargs.get('select_root_type', 'random').lower() # other options are min, max to select the lowest prob and highest prob state, respectively, for the root sequence.
        assert(self.select_root_type in ["random", "min", "max"]), "\nValue for keyword argument select_root_type argument must be either 'random', 'min', or 'max'. Default behavior is random."



    def _setup_subcounts(self):

        if self._code == MOLECULES.nucleotides:
            self.branch_substitution_counts["nucleotide"] = {}
        elif self._code == MOLECULES.amino_acids:
            self.branch_substitution_counts["amino_acid"] = {}
        elif self._code == MOLECULES.codons:
            self.branch_substitution_counts["nucleotide"] = {}
            self.branch_substitution_counts["amino_acid"] = {} ## aka also covers nonsynonymous
            self.branch_substitution_counts["synonymous"] = {}
            self.branch_substitution_counts["codon"]      = {}
        else:
            self.branch_substitution_counts["custom"]     = {}



    def _setup_partitions(self):
        '''
            Setup and various sanity checks.
        '''
        # If partitions is not a list but indeed a Partition, turn into a list. If not a partition, assert.
        assert(self.partitions is not None), "\n\nNo partitions were provided to Evolver. Please specify partition(s) with the keyword argument 'partitions' or 'partition' (both are accepted)."
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
                assert(self.full_tree.model_flag is not None), "\n\n Your root model name does not correspond to any of the Model objects' names provided to your Partition object(s)."
            self._root_seq_length += sum( part.size )

        # Final check on size
        assert(self._root_seq_length > 0), "\n\nPartitions have no size!"




    def __call__(self, **kwargs):
        '''
            Simulate sequences, perform any necessary post-processing, and save sequences and/or other info to appropriate files.

            Optional keyword arguments:
                1. **seqfile** is a custom name for the output simulated alignment. Provide None or False to suppress file creation.
                2. **seqfmt**  is the format for seqfile (either fasta, nexus, phylip, phylip-relaxed, stockholm, etc. Anything that Biopython can accept!!) Default is FASTA.
                3. **ratefile** is a custom name for the "site_rates.txt" file. Provide None or False to suppress file creation.
                4. **infofile** is a custom name for the "site_rates_info.txt" file. Provide None or False to suppress file creation.
                5. **write_anc** is a boolean argument (True or False) for whether ancestral sequences should be output along with the tip sequences. Default is False.
                6. **scale_tree** is a float argument for scaling the entire tree by a certain factor. Note that this argument can alternatively be used in the newick module (with `read_tree`) function, but it is included here for ease in replicates (e.g. lots of sims along same tree w/ varied branch lengths). Default: 1.
                7. **countfile** is a file to which you can save the total substitution counts per branch as a CSV. This will be fairly inaccurate unless one uses `algorithm=1` in which case it is guaranteed to be precisely accurate. Default: None (ie not exported)
                8. **seed**, a float to set your own random seed
                9. **algorithm**, 0 for exponentiation, 1 for Gillespie

            Examples:
                .. code-block:: python

                   >>> # Set up Evolver instance, and use below in various ways
                   evolve = Evolver(tree = my_tree, partitions = my_partition_list)

                   >>> # Evolve according to default settings
                   >>> evolve()

                   >>> # Include ancestral sequences in output file
                   >>> evolve(write_anc = True)

                   >>> # Custom sequence file name and format, and suppress rate information
                   >>> evolve(seqfile = "my_seqs.phy", seqfmt = "phylip", ratefile = None, infofile = None)

        '''
        # Input arguments
        self.seqfile    = kwargs.get('seqfile', 'simulated_alignment.fasta')
        self.seqfmt     = kwargs.get('seqfmt', 'fasta').lower()
        self.write_anc  = kwargs.get('write_anc', False)
        self.ratefile   = kwargs.get('ratefile', 'site_rates.txt')
        self.infofile   = kwargs.get('infofile', 'site_rates_info.txt')
        self.countfile  = kwargs.get('countfile', None)
        self.scale_tree = kwargs.get('scale_tree', 1.)
        self.seed       = kwargs.get('seed', None)
        self.algorithm  = kwargs.get('algorithm', 0) ## 0 is exponentiation; 1 is gillespie/jump-chain

        #### SET SEED ANEW ####
        if self.seed is not None:
            assert(isinstance(self.seed, int)), "\nERROR: Provided seed must be an integer."
            self.rng = np.random.default_rng(seed = self.seed)
        else:
            self.rng = np.random.default_rng()

        ### Algorithm check
        try:
            self.algorithm = int(self.algorithm)
        except:
            raise AssertionError("[ERROR]\n algorithm must be 0 (exponentiation) or 1 (Gillespie).")
        assert(self.algorithm in [0,1]), "[ERROR]\n algorithm must be 0 (Q matrix exponentiation) or 1 (simulate waiting times with Gillespie i.e. jump chain algorithm)."

        ## Setup substitution counting dictionaries. MUST happen here and not in __init__
        self.branch_substitution_counts = {}
        self._setup_subcounts()

        ### Warn about algorithm/countfile
        if self.countfile is not None and self.countfile is not False and self.algorithm == 0:
            print("\nWARNING: You are exporting substitution counts but using matrix exponentiation as the simulating algorithm. This means counts will most often be an _underestimate_. Consider setting `algorithm=1` when calling your Evolver instance to get _accurate_ counts.")


        # Simulate recursively
        self._sim_subtree(self.full_tree)


        # Shuffle sequences?
        self._shuffle_sites()

        # Convert Site dictionaries to sequence dictionaries
        self.leaf_seqs = self._convert_site_to_seq_dict(self._leaf_sites)
        self.evolved_seqs = self._convert_site_to_seq_dict(self._evolved_sites)

        # Save rate, count info, as needed
        if self.ratefile:
            self._write_ratefile()
        if self.infofile:
            self._write_infofile()
        if self.countfile:
            self._write_countfile()


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
            if part._shuffle:
                size = sum( part.size )
                part_pos = np.arange( size ) + start
                self.rng.shuffle(part_pos)
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
        from Bio import SeqIO

        alignment = [] # list of seqobjects
        ## Biopython 1.77 --> 1.78 deprecated the alphabet module. See pyvolve issue #21
        ## Hell hack ensues:
        from Bio import __version__ as bio_version
        bio_version = float(bio_version)
        if bio_version >= 1.77:
            for entry in seqdict:
                seq_entry = Seq( seqdict[entry])
                seq_object = SeqRecord(seq_entry, id = entry, description = "", annotations={"molecule_type": "DNA"})
                # It truly doesn't matter what type since will be written to file. molecule_type is only one of DNA, RNA, protein.
                alignment.append(seq_object)
        else:
            for entry in seqdict:
                seq_entry = Seq( seqdict[entry] )
                seq_object = SeqRecord(seq_entry, id = entry, description = "")
                alignment.append(seq_object)

        try:
            SeqIO.write(alignment, self.seqfile, self.seqfmt)
        except:
            raise TypeError("\n Output file format is unknown. Consult with Biopython manual to see which I/O formats are accepted.\n NOTE: If you are attempting to save as phylip and are receiving this error, try seqfmt = 'phylip-relaxed' instead.")



    def _write_ratefile(self):
        '''
            Write ratefile, a tab-delimited file containing site-specific rate information. Considers leaf sequences only.
            Writes -   Site_Index    Partition_Index     Rate_Category
            All indexing is from *1*.
        '''
        refseq = list(self._leaf_sites.values())[0]
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
                        if m.model_type.lower() in ["mg", "gy"]:
                            if m.is_hetcodon_model():
                                infof.write(outstr + str(round(m.params['beta'][r],4)) + "," + str(round(m.params['alpha'][r],4)) )
                            else:
                               infof.write(outstr + str(round(m.params['beta'],4)) + "," + str(round(m.params['alpha'],4)) )
                        else:
                            infof.write(outstr + str(round(m.rate_factors[r],4)) )


    def _write_countfile(self):
            '''
                Write countfile, a comma-delimited file which shows the number of substitutions per branch, INCLUDING internal branches.
            '''
            outstring = "substitution_type,branch_name,count\n"
            for substitution_type in self.branch_substitution_counts:
                for branch_name in self.branch_substitution_counts[substitution_type]:
                    n = str( self.branch_substitution_counts[substitution_type][branch_name] )
                    outstring += substitution_type + "," + branch_name + "," + n + "\n"
            with open(self.countfile, 'w') as countf:
                countf.write(outstring.strip())




    def get_sequences(self, anc = False):
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
    def _count_substitutions_branch(self, parent, child, branch_name):
        """
            Tabulate the total number of substitutions between parent and child along a given branch.
            Required keyword arguments:
                1. **parent** Site object of parent
                2. **child**  Site object of child
                3. **branch_name** name of branch where counting occurs. Will be used as key in final dictionary
        """
        parent_seq = self._site_to_sequence(parent)
        child_seq = self._site_to_sequence(child)

        ## nucleotide, amino_acid, or custom ONLY
        if len( list(self.branch_substitution_counts.keys()) ) == 1:

            codekey = list(self.branch_substitution_counts.keys())[0]
            total_changes = sum( [1 for i in range(len(parent_seq)) if parent_seq[i] != child_seq[i]] )
            try:
                self.branch_substitution_counts[codekey][branch_name] += total_changes
            except:
                self.branch_substitution_counts[codekey][branch_name] = total_changes

        ## codons
        else:
            code_changes = {"synonymous": 0, "codon": 0, "nucleotide": 0, "amino_acid": 0}

            for i in range(0, len(parent_seq), 3):
                parent_codon = parent_seq[i:i+3]
                child_codon  = child_seq[i:i+3]
                parent_aa    = MOLECULES.codon_dict[parent_codon] #MOLECULES.amino_acids.index( MOLECULES.codon_dict[parent_codon] )
                child_aa     = MOLECULES.codon_dict[child_codon] # MOLECULES.amino_acids.index( MOLECULES.codon_dict[child_codon] )
                nuc_diff = sum([1 for x in range(3) if parent_codon[x] != child_codon[x]])

                code_changes["nucleotide"] += nuc_diff
                if parent_codon != child_codon:
                    code_changes["codon"] += 1  ## keep track still for sanity
                    if parent_aa == child_aa:
                        code_changes["synonymous"] +=1
                    else:
                       code_changes["amino_acid"] += 1

            assert(code_changes["amino_acid"] + code_changes["synonymous"] == code_changes["codon"]), "\n[ERROR] Improper tabulation of branch substitution codon counts."

            for code_key in code_changes:
                try:
                    self.branch_substitution_counts[code_key][branch_name] += code_changes[code_key]
                except:
                    self.branch_substitution_counts[code_key][branch_name] = code_changes[code_key]



    def _exponentiate_matrix(self, Q, t):
        '''
            Perform exponentiation on instantaneous matrix to produce transition matrix, from a given instantaneous matrix Q and a given branch length t.
            Assert that all rows sum to 1.
            Return P
        '''
        P = linalg.expm( np.multiply(Q, float(t) ) )
        assert( np.allclose( np.sum(P, axis = 1), np.ones(len(self._code))) ), "Rows in transition matrix do not each sum to 1."
        return P





    def _obtain_model(self, part, flag):
        '''
            Obtain the appropriate Model object for evolution along a particular branch.
        '''
        my_model = None
        if part.branch_het():
            for m in part.models:
                if m.name == flag:
                    my_model = m
                    break
        else:
            my_model = part.models[0]
        assert( my_model is not None ), "\n\nCould not retrieve model a particular branch's evolution. Double-check that all models used in simulation have been properly given to a partition.\n If you are using branch heterogeneity, make sure that your model names match the model flags in the tree."
        return my_model





    def _generate_prob_from_unif(self, prob_array):
        '''
            Sample a sequence (nuc,aa,or codon), and return an integer for the sequence chosen from a uniform distribution.
            Arugment *prob_array* is any list and/or numpy array of probabilities which sum to 1.
        '''
        assert ( abs(np.sum(prob_array) - 1.) < ZERO), "Probabilities do not sum to 1. Cannot generate a new sequence."
        r = self.rng.uniform(0,1)
        i = 0
        sum = prob_array[i]
        while sum < r:
            i += 1
            sum += prob_array[i]
        return i


    def _assign_root_seq_from_MRCA(self, raw_MRCA):
        '''
            Assign a root sequence from provided MRCA. This function converts a provided MRCA into a list of Site objects.
        '''
        MRCA_sites = []
        step = len(self._code[0])
        for i in range(0, len(raw_MRCA), step):

            # Initialize a Site object
            new_site = Site()
            new_site.rate = 0

            # Convert to int_seq and append to MRCA_sites
            try:
                new_site.int_seq = self._code.index( raw_MRCA[i:i+step] )
            except:
                raise ValueError("\n\nProvided root sequence does not have the same code (alphabet) as model. Remove all noncanonical and/or wrong letters from provided root sequences. Further, if you are specifying codons, ensure that the length of your root sequence is divisible by 3.")
            MRCA_sites.append(new_site)

        assert( len(MRCA_sites)*step == len(raw_MRCA)), "\n\nRoot sequence improperly converted."
        return MRCA_sites







    def _generate_root_seq(self):
        '''
            Generate a root sequence based on the stationary frequencies.
            Return a complete root sequence list of Site objects.

            NOTE: The select_root_type attribute is for the sitewise_dnds_mutsel project and was created on 4/30/15.
        '''

        root_sequence = [] # This will contain a list for each partition's sequence (which is itself a list of Site() objects)

        for part in self.partitions:

            part_root = [] # reset to empty list

            # Is there a root sequence?
            if part.MRCA is not None:
                part_root = self._assign_root_seq_from_MRCA(part.MRCA)

            # No root sequence provided. Must generate one.
            else:

                # Grab model info for this partition to get frequency vector for root simulation
                root_model = self._obtain_model(part, self.full_tree.model_flag)

                # Generate root_sequence and assign the Site a rate class
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
            Function to traverse a Node (tree) object recursively and simulate sequences.
            Required positional arguments include,
                1. **current_node** is the node (either internal node or leaf) TO WHICH we evolving
                2. **parent_node** is the node we are evolving FROM. Default of None is only called when the root sequence is not yet made.
        '''

        # We are at the base and must generate root sequence
        if (parent_node is None and current_node.root is True):
            current_node.seq = self._generate_root_seq() # the .seq attribute is actually a list of Site() objects.

        else:
            assert(current_node.root is False), "\n\n [ERROR]: Non-root node interpreted as root."
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
        assert (parent_node.seq is not None), "\n\nThere is no parent sequence from which to evolve!"
        assert (current_node.branch_length >= 0.), "\n\n Your tree has a negative branch length. I'm going to quit now."
        if current_node.model_flag is None:
            current_node.model_flag = parent_node.model_flag


    def _make_jump_transition_matrix(self, Q_matrix):
        '''
            Return transition matrix for jump chain (Gillespie algorithm) from initial Q matrix
        '''
        jump_matrix = np.copy(Q_matrix)
        for i in range(len(Q_matrix[0])):
            denominator = -1. * Q_matrix[i][i]
            jump_matrix[i] /= denominator
            jump_matrix[i][i] = 0.
            assert( (1. - np.sum(jump_matrix[i])) <= ZERO), "\n\n[ERROR] Bad jump chain transition matrix calculation."
        return(jump_matrix)



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

            for new_seq_part in new_seq:
                self._count_substitutions_branch(new_seq_part, new_seq_part, current_node.name)

        else:
            branch_length = self.scale_tree * float(current_node.branch_length)
            new_seq = []

            for p in range( len(self.partitions) ):
                # Obtain current model for this partition at this branch
                part = self.partitions[p]
                current_model = self._obtain_model(part, current_node.model_flag)
                index = 0
                part_new_seq = []  # will temporarily store this partition's new sequence



                for i in range( current_model.num_classes() ):
                    # Grab instantaneous rate matrix, which is done differently depending if codon (dN/dS) model or not. This is the rate het in the partition.
                    Q_matrix = None
                    if current_model.is_hetcodon_model():
                        Q_matrix = current_model.matrix[i]
                    else:
                        Q_matrix = current_model.matrix * current_model.rate_factors[i] # note that rate_factors = [1.] if no site heterogeneity, so matrix unchanged
                    assert( Q_matrix is not None ), "\n\nCouldn't retrieve instantaneous rate matrix."



                    ########################## SWAP ALGORITHMS HERE ############################
                    part_parent_seq = parent_node.seq[p][index : index + part.size[i]]
                    if self.algorithm == 0:

                        # Generate transition matrix
                        P_matrix = self._exponentiate_matrix(Q_matrix, branch_length)

                        # Evolve branch
                        for j in range( part.size[i] ):
                            new_site = deepcopy( part_parent_seq[j] )
                            new_site.int_seq = self._generate_prob_from_unif( P_matrix[ new_site.int_seq ] )
                            part_new_seq.append( new_site )
                            index += 1
                        self._count_substitutions_branch(part_parent_seq, part_new_seq, current_node.name)

                    elif self.algorithm == 1:
                        ## Each site w/ this Q has the same rate matrix (we evolve each rate chunk separately)
                        ## So we can simply loop over sites and apply each matrix accordingly. Will probably be inefficient but i am not a computer scientist.

                        jump_transition_matrix = self._make_jump_transition_matrix(Q_matrix)


                        for j in range( part.size[i] ):
                            entered_while = False
                            last_jump_site = deepcopy( part_parent_seq[j] )
                            new_site       = deepcopy( part_parent_seq[j] )

                            remaining_branch_length = deepcopy(branch_length)
                            waiting_scale = -1 * Q_matrix[last_jump_site.int_seq][last_jump_site.int_seq]
                            waiting_time  = self.rng.exponential(scale = waiting_scale)

                            while waiting_time < remaining_branch_length:
                                entered_while = True
                                new_site.int_seq = self._generate_prob_from_unif( jump_transition_matrix[ last_jump_site.int_seq ] )
                                remaining_branch_length -= waiting_time

                                waiting_scale = -1 * Q_matrix[new_site.int_seq][new_site.int_seq]
                                waiting_time  = self.rng.exponential(scale = waiting_scale)

                                self._count_substitutions_branch(last_jump_site, new_site, current_node.name)
                                last_jump_site = deepcopy(new_site)

                            if not entered_while:
                                assert(last_jump_site.int_seq == new_site.int_seq), "\n[ERROR] Gillespie algorithm bug; please file an issue with reproducible example."
                                self._count_substitutions_branch(last_jump_site, new_site, current_node.name)

                            part_new_seq.append( new_site )
                            index += 1

                new_seq.append( part_new_seq )

        return new_seq
