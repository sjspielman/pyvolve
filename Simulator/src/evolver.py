
from indel import *
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
        
        #####################################################################
        # Note that these 3 lines only make sense when there are no indels. #
        self._root_seq_length = 0
        for i in range(self._number_partitions):
            self.root_seq_length += self._partitions[i][0]  
        #####################################################################     
        self.set_code()

        
    def _set_code(self):
        ''' Determine whether we have nuc, amino, or codon alphabet and assign a genetic code. '''   
        # An incredible hack! Go python, go!!
        dim = self._partitions[0][1].values()[0].Q.shape[0]
        
        if dim == 4:
            self._code = MOLECULES.nucleotides
        elif dim == 20:
            self._code = MOLECULES.amino_acids
        elif dim == 61:
            self._code = MOLECULES.codons
        else:
            raise AssertionError("wut.")

    ## definitely deprecated, or if not already will be soon.
    def write_sequences(self, **kwargs):
        ''' Write resulting sequences to a file, currently only in fasta format. '''
        outfile  = kwargs.get("outfile", "seqs_"+strftime("%m.%d.%H.%M.%S")+".fasta")  
        out_handle=open(outfile, 'w')
        for entry in self.alndict:
            seq = self._intseq_to_string(self.alndict[entry])
            out_handle.write(">"+entry+"\n"+seq+"\n")
        out_handle.close()               
            
    def _sequence_to_integer(self, entry):
        ''' Take a DNA/PROT value and return its integer index '''
        return self._code.index(entry)
    
    def _integer_to_sequence(self, index):
        ''' Take a DNA/PROT index and return its corresponding molecule ''' 
        return self._code[index]
    
    ######### FUNCTION MAY NEED OVERHAUL TO ACCOUNT FOR SITE CLASS ###############
    def _intseq_to_string(self, intseq):
        ''' Take a sequence coded as ints and turn to actual molecule string '''
        stringseq = ''
        for i in intseq:
            stringseq += self._integer_to_sequence(i)
        return stringseq   
    ######################################################################################   

    def _check_parent_branch(self, node, parent_seq, parent_model):
        ''' Check that the parent sequence exists, branch length is reasonable, and assign a model. ''' 
        assert (parent_seq != None), "There is no parent sequence."

        branch_model = node.modeletion_flag
        if branch_model is None:
            node.modeletion_flag = parent_model

        branch_length = float( node.branch_length )
        assert (branch_length >= 0.), "Branch length is negative. Must be >= 0."
        return branch_length, node.modeletion_flag           
        

    def _generate_site(self, node_name, state):
        ''' Generate a new Site object when generating the root (state 0) and insertion (state 1) sites. '''
        site        = misc.Site()
        site.origin = node_name
        site.state  = state
        return site
           
  
    def _generate_root_seq(self):
        ''' Select starting sequence based on state frequencies, for each partition, and return full root sequence. '''
        
        current_node.seq = []
        for n in range(self._number_partitions):
            partlen = self._partitions[n][0]
            partition_sites = []
            freqs = self._partitions[n][1][self._root_model].subst_params['stateFreqs']
            index = 0
            for j in range(partlen):
                newsite = _generate_site(current_node.name, 0)
                newsite.int_seq = self._generate_prob_from_unif(freqs)
                partition_sites[index] = newsite
                index += 1
            current_node.seq.append(partition_sites)


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
        
        
    def simulate(self, current_node, parent_node = None):
        ''' Traverse the tree and simulate. '''

        # We are at the base and must generate root sequence
        if (parent_node is None):
            current_node.seq = self._generate_root_seq() 
            current_node.modeletion_flag = self._root_model 
        else:
            self._simulate_substitution(current_node, parent_node)
            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                self.simulate(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[current_node.name]=current_node.seq
            
            
    def _simulate_substitution(self, node, parent_node):
        ''' Crux function to evolve sequences along a branch, without indels.'''
    
        # Ensure brank length ok, parent sequence exists, and model is assigned.
        parent_seq = parent_node.seq
        parent_model = parent_node.modeletion_flag
        branch_length, branch_model = self._check_parent_branch(node, parent_seq, parent_model)

        # Evolve only if branch length is greater than 0.
        if branch_length <= ZERO:
            new_sequence = parent_seq
        else:
            # Evolve each partition, i, and add onto new_sequence as we go
            new_sequence = np.empty(self.root_seq_length, dtype=int)
            index = 0
            for i in range(self._number_partitions):
            
                # set the length and the instantaneous rate matrix for this partition at this node
                seqlen  = self._partitions[i][0]
                inst_matrix = self._partitions[i][1][branch_model].Q
                #print branch_model, node.modeletion_flag, self._partitions[i][1][branch_model].subst_params['beta']
                
                # Generate probability matrix for evolution along this branch and assert correct
                Qt = np.multiply(inst_matrix, branch_length)
                prob_matrix = linalg.expm( Qt ) # Generate P(t) = exp(Qt)
                for i in range(len(self._code)):
                    assert( abs(np.sum(prob_matrix[i]) - 1.) < ZERO ), "Row in P(t) matrix does not sum to 1."
    
                # Move along parent_seq and evolve. 
                for j in range(seqlen):
                    new_sequence[index] = self._generate_seq( prob_matrix[parent_seq[index]] )
                    index+=1
                             
        # Attach final sequence to node
        node.seq = new_sequence 


            
            
            
            
            
            
            
            
            
            
            
            
class IndelEvolver(Evolver):
    def __init__(self, *args):
    
              
    def simulate(self, current_node, parent_node = None, wait_time = None):
        ''' Traverse the tree and simulate, with indels. '''

        # We are at the base and must generate root sequence
        if (parent_node is None):
            self._generate_root_seq( current_node )
            current_node.modeletion_flag = self._root_model 
        else:
            wait_time = self._evolve_branch(current_node, parent_node, wait_time)
            
        # We are at an internal node. Keep evolving
        if len(current_node.children)>0:
            for child_node in current_node.children:
                # only give wait_time to one child.
                if child_node is current_node.children[0]:            
                    self.simulate(child_node, current_node, wait_time)
                else:
                    self.simulate(child_node, current_node)
                
        # We are at a leaf. Save the final sequence
        else: 
            self.alndict[current_node.name]=current_node.seq          
            
    def _simulate_substitution(self, parent_seq, model, inserted_sites, branch_length):
        ''' Simulate substitution process along a branch for a single partition.
            Here, we will only be tweaking the Site.int_seq attributes (no others!!).
         '''
        ###################### NEEDS OVERHAUL TO EVOLVE NEWLY INSERTED BITS SEPARATELY FROM PREVIOUSLY PRESENT BITS #############################
                                    
        # Generate probability matrix for evolution along this branch and assert correct
        Q = model.Q
        Qt = np.multiply(Q, branch_length)
        prob_matrix = linalg.expm( Qt )
        for i in range(len(self._code)):
            assert( abs(np.sum(prob_matrix[i]) - 1.) < ZERO ), "Row in P(t) matrix does not sum to 1."
     
        
        # March along parent_seq and substitute
        partition_sites = []
        for j in range( parent_seq.numsites ):
            
            # Assign parent position to current node.
            partition_sites[j] = parent_seq[j]
            
            # If position is not a gap, evolve it
                partition_sites[j].int_seq = self._generate_prob_from_unif( prob_matrix[ partition_sites[j].int_seq ] )     

        return partition_sites   
        
        
    def _evolve_branch(self, node, parent_node, wait_time = None):
        ''' Crux function to evolve sequences along a branch. '''
    
        # Ensure brank length ok, parent sequence exists, and model is assigned.
        parent_seq = parent_node.seq
        parent_model = parent_node.modeletion_flag
        branch_length, branch_model = self._check_parent_branch(node, parent_seq, parent_model)

        # Evolve only if branch length is greater than 0.
        if branch_length <= ZERO:
            new_sequence = parent_seq
        else:
        
            new_sequence = [] # Since indels, must be dynamic list.
            for n in range(self._number_partitions):
                model = self._partitions[n][1][branch_model] # variable for clarity
                
                # Simulate indels and create the partition_sites for this partition
                partition_sites, inserted_sites, next_wait_time = self._simulate_indel(parent_seq, node.name, model, branch_length, wait_time)
                
                # Simulate substitutions.
                partition_sites = self._simulate_substitution(parent_seq, model, inserted_sites, branch_length)
            
                
                new_sequence.append(partition_sites)
                
        # Attach final new_sequence list of evolved partitions to node
        node.seq = new_sequence 
        
        return next_wait_time
        
        
        
    ############################### INDEL-SPECIFIC CODE ##########################################
            
    def _calc_prob_indeletion_event(insert_space, deletion_space, insertion_rate, deletion_rate, mean_deletion_size):
        ''' calculate probability of indel event occuring, which is a fxn of the current sequence length.
            NOTE: as length will stay the same from a parent to child branch, it's ok that things get carried over.
            # current scheme, as of 7/7/14: P_event = rate_ins*(length+1) + rate_del*(length) . <- iSG uses this.
            # ME: only real sequence bits can be deleted, as in cannot delete a gap. Therefore, the probability of deletions depends not on a naive length, but number of sequence sites that there are. Otherwise will overestimate deletion rate.
        '''
        insert_space = float(insert_space)
        deletion_space = float(deletion_space)
        p_event = insertion_rate*(insert_space+1.) + deletion_rate*(deletion_space + mean_deletion_size + 1.) 
        prob_insertion = (insertion_rate*(insert_space+1.)) / p_event
        prob_deletion = deletion_rate*deletion_space / p_event
        
        return p_event, prob_insertion, prob_deletion
        
            
    def _should_I_insert(self, prob_insertion, prob_deletion):
        ''' random number draw to determine if insertion or deletion will occur 
        '''
        r = rn.uniform(0,1)
        if r <= prob_insertion:
            return True
        else:
            return False
    
          
    def _generate_insertion(self, model, numsites, node_name):
        ''' 1. Get a length for the insertion
            2. Generate the inserted sequence using base frequencies
            3. Decide on a location for the inserted sequence
        '''
        # Generate a length
        length = self._generate_prob_from_unif(model.indeletion_params['insDist']) + 1 # +1 since lengths aren't indexed from zero, and indel length distributions begin at 1.

        # Generated new inserted sites and hold in list called insertion
        freqs = model.subst_params['stateFreqs']
        insertion = []
        for j in range(length):
            newsite = _generate_site(node_name, 1)
            newsite.int_seq = self._generate_prob_from_unif(freqs)
            insertion.append(newsite)
            
        # Generate a location
        location = rn.uniform(0, numsites)
        
        return location, insertion
        
        
        
        
    
    def _generate_deletion(self, model, deletable_sites):
        ''' 1. generate length for deletion
            2. generate starting position for deletion, conditioned on length
        '''
         # Generate a length
        length = self._generate_prob_from_unif(model.indeletion_params['delDist']) + 1 # +1 since lengths aren't indexed from zero, and indel length distributions begin at 1.

        # Generate a location
        location = rn.randint(0, len(deletable_sites) - length)
        return location, length
        
        
    def _update_deletable(self, deletable_sites, seq):
        ''' edit list which contains sites that are ok to delete. this needs to be updated as indels are simulated.
        '''
        for n in range(len(seq)):
            if seq[n].int_seq != -1
                deletable_sites.append(n)
        return deletable_sites
        
            
    def _simulate_indel(self, parent_seq, node_name, model, bl, wait_time)
        ''' ABOVE ARGUMENTS NEED OVERHAUL.
            Go, indels, go go go!! Welcome to Gillespieeeeeeeeee!!!! Wooahhhooahhh hey indels! Yeah, yeah, yeah! YEAH YEAH YEAH!     
        '''
        # Tracking variables. Not for whole simulation, just for this branch.
        inserted_sites = [] # Will contain list of tuples of len=2. In each tuple, first entry are positions that were inserted for a given insertion event, and second entry is their remaining time (to later use for p=e^(Qt) for those positions).
        deletable_sites = [] # Make sure that we don't delete gaps
        deletable_sites = self._update_deletable(deletable_sites, parent_seq.int_seq)
    
        # Scale insertion and deletion by branch_length
        insertion_rate = model.indeletion_params['insertion_rate'] * bl
        deletion_rate = model.indeletion_params['deletion_rate'] * bl

  
        # Generate event probabilities and waiting time (if the latter doesn't already exist from parent branch) 
        p_event, prob_insertion, prob_deletion = self._calc_prob_indeletion_event(seqlen, numsites, insertion_rate, deletion_rate, model.indeletion_params['meanDelLen'])
        if first_wait_time is None:
            wait_time = np.random.exponential(scale = 1./p_event), p_event            
        else:
            wait_time = first_wait_time
        
        ## SIMULATE ##
        remaining_time = branch_length - wait_time
        while remaining_time >= ZERO:
        
            # Perform insertion. track it and add it to sequence.
            if self._should_I_insert(prob_insertion, prob_deletion):
                insertion_location, insertion = self._generate_insertion( model, numsites, node_name )
                inserted_sites.append( (range(insertion_location, len(insertion)), remaining_time) ) # We can always update this in case the inserted sequences are deleted, but I suspect we don't have to because the _simulate_substitution function will skip these sites.
                for entry in insertion:
                    partition_sites.insert(insertion_location, entry)
                partition_sites.insert(insertion_location, insertion) 
              
            # Perform deletion and edit relevant site states
            else:
                deletion_length, deletion_location = self._generate_deletion( model, deletable_sites )
                deleteme = deletable_sites[deletion_location: deletion_location + deletion_length] # indices of partition_sites that we need to switch to deletions
                for x in deleteme:
                    # Change intSeq attr to None (gap). Change state attr from 0->2 (core to deleted core) or from 1->3 (insertion to deleted insertion)
                    assert(partition_sites[x].int_seq is not None), "You've just deleted a gap! WOAH NOW!"
                    assert(partition_sites[x].state == 1 or partition_sites[x].state == 3), "intSeq is not a deletion, but its state is a deletion. WOAH NOW AGAIN!"
                    partition_sites[x].int_seq = None
                    partition_sites[x].state += 2
            
            # update deletable_sites list
            deletable_sites = self._update_deletable(deletable_sites, partition_sites)          
         
            # next leap
            wait_time = self._generate_wait_time( len(partition_sites), insertion_rate, deletion_rate)
            remaining_time -= wait_time

    return partition_sites, inserted_sites, abs(remaining_time)
            
            
       