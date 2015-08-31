'''
    SJS.
    This script implements a module counting site-specific dN/dS values, as simulated by pyvolve. The counting method implemented is similar to the SLAC method [Kosakovsky Pond & Frost (2005) MBE], although as the actual ancestral sequences are known, no reconstruction is performed.
    NOTE: this module is meant to be used ONLY with sequences simulated with pyvolve and thus can't handle gaps or ambiguities.
    
    Usage:
        from count_simulated_dnds import *
    
        c = dNdS_Counter(<alnfile>, <treefile>, <mutation_dictionary>)
        # alnfile: a FASTA-formatted sequence alignment file produced by a pyvolve simulation. This file **MUST** contain ancestral sequences (get this alignment by including the argument `write_anc=True` when calling a pyvolve Evolver class).
        # treefile: a file containing the *exact* newick tree used in the pyvolve simulation which produced the alnfile.
        # mutation_dictionary: **optional** argument indicating the mutation rates between nucleotides. This argument is analogous to the "mu" dictionary provided to pyvolve when simulating with custom mutation rates. If this argument is not provided, this module assumes equal mutation rates (e.g. a JC69 situation).
        
        c.calculate_dnds()
        # This method will compute site-specific dN/dS values and output a tab-delimited file of this format:
        # site  ns_changes  s_changes   ns_sites    s_sites
        #  ...      ...         ...        ...         ... 
          
        # You can then calculate dN/dS from this file with this calculation: (ns_changes/ns_sites) / (s_changes/s_sites) . This is straight-forward in something like R or the python pandas package. Beware your parentheses placement, and also keep a look-out for NA or Inf values (these can happen when no synonymous changes occurred).
        # By default, the output file will be called "counted_dnds.txt". To change this name, include the argument savefile, e.g. c.calculate_dnds(savefile = "my_preferred_filename.txt")

    
    
    Please post all questions, bugs, etc. to https://github.com/sjspielman/pyvolve/Issues
    
'''
from copy import deepcopy
from pyvolve import newick
from Bio import AlignIO
import numpy as np
import sys


class CalcOverBranch(object):
    '''
        Compute quantities for dN/dS over a branch on a phylogeny, using a SLAC approach.
        Takes 4 ordered input arguments:
            source        = the source codon (string)
            target        =  the target codon (string
            branch_length = the branch length (float)
            nuc_sub_probs = a dictionary of normalized (sum to 1) nucleotide substitution probabilities
    '''
    def __init__(self, source, target, branch_length, nuc_sub_probs):
        
        # Genetics variables
        self.nucleotides = ["A", "C", "G", "T"]
        self.stop_codons = ["TAA", "TAG", "TGA"] 
        self.codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.translation_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
        
        # Input arguments
        self.source_codon = source # originating codon
        self.target_codon = target # final codon
        self.bl = branch_length  # branch length
        self.B = nuc_sub_probs   # dictionary of nucleotide substitution probabilities
        assert(self.source_codon in self.codons and self.target_codon in self.codons), "\n\nImproper codons provided for dN/dS calculation. They are either stop codons, or straight-up not DNA."

        
        self.path_codon_counts = {} # dictionary to contain the counts for all codons appearing in the paths between source, target
        self.n_changes = 0.
        self.s_changes = 0.
        self.n_sites = 0.
        self.s_sites = 0.
        

    def compute_branch_dnds(self):
        '''
            *This* is the function users should call to compute nonsyn, syn changes and sites over a branch.
        '''
        
        self.compute_paths_and_changes()
        self.count_sites_over_branch()
    
        assert(self.n_sites + self.s_sites > 0.), "\n\nNo sites were counted for this branch. Bad news."
        
        return self.n_changes, self.s_changes, self.n_sites, self.s_sites


    def num_nuc_diff(self, s, t):
        '''
            Count the number of nucleotide differences between two codons: s,t.
        '''
        return sum([1 for i in range(3) if s[i] != t[i]])
       
       
       
    def eval_path(self, l, n_count, s_count):
        '''
            Enumerate and evaluate the changes within a path.
        '''
        stop = False
        temp_s_count = 0.
        temp_n_count = 0.
        path = [self.source_codon]
        source = deepcopy(self.source_codon)

        for p in l:
            new_source = source[0:p] + self.target_codon[p] + source[p+1:3]

            # If our path takes us through a stop, we will disregard the whole path
            if new_source in self.stop_codons:
                stop = True
                break

            # Evaluate the change as syn, nonsyn
            if self.translation_dict[new_source] == self.translation_dict[source]:
                temp_s_count += 1.
            else:
                temp_n_count += 1.
            path.append(new_source)
            source = new_source
    
        # Save the path to count dictionary if accessible
        if not stop:
            path.append(self.target_codon) # Append the final target codon to the list
            for entry in path:
                if entry in self.path_codon_counts:
                    self.path_codon_counts[entry] += 1
                else:
                    self.path_codon_counts[entry] = 1

        # Tally and return
        s_count.append(temp_s_count)
        n_count.append(temp_n_count)
        return n_count, s_count
       
       
    def compute_paths_and_changes(self):
        '''
            Assess the shortest paths between self.source and self.target which don`t pass through stop codons.
            Builds self.path_codon_counts and computes self.n_changes, self.s_changes
        '''
    
        s_count = [] # Total number of synonymous changes. Use list for each averaging.
        n_count = [] # Total number of nonsynonymous changes. Use list for each averaging.
    
        # Which sites had a mutation?
        diff_sites = [i for i in range(3) if self.source_codon[i] != self.target_codon[i]]
        n = len(diff_sites)
        
        # No change?
        if n == 0:
            self.path_codon_counts[self.source_codon] =  1.
        
        # Single change?
        if n == 1:
            if self.translation_dict[self.source_codon] == self.translation_dict[self.target_codon]:
                self.s_changes = 1.
            else:
                self.n_changes = 1.
            self.path_codon_counts[self.source_codon] = 1.
            self.path_codon_counts[self.target_codon] = 1.
    
        # Multiple changes?
        elif n > 1:
            for i in range(n):
                m = diff_sites[i]
                new_diff_sites = diff_sites[0:i] + diff_sites[i+1:n]
            
                n_count, s_count = self.eval_path(new_diff_sites, n_count, s_count)
                # Evaluate reverse substitution order if needed
                if len(new_diff_sites) == 2:
                    new_diff_sites.reverse()
                    n_count, s_count = self.eval_path(new_diff_sites, n_count, s_count)
            
            # Average the changes
            self.n_changes = sum(n_count) / float(len(n_count))
            self.s_changes = sum(s_count) / float(len(s_count))
       
           



    def count_codon_sites(self,codon):
        '''
            Count number of synonymous, nonsynonymous sites in a given codon.
         
            Argument "codon" is the sense codon of interest.
        '''
    
        codon = codon.upper()
        source_aa = self.translation_dict[codon]
    
        # Determine which sense codons are a single nucleotide change from input argument, codon
        # One list for nonsyn and one list for syn
        n_sites = 0.
        s_sites  = 0.
    
        for i in range(3):
            s_numer = 0.
            n_numer = 0.
            denom = 0.
            for n in self.nucleotides:
                if codon[i] != n:
            
                    target = codon[0:i] + n + codon[i+1:3]
                    if target not in self.stop_codons:
                        target_aa = self.translation_dict[target]
                    
                        # Evaluate syn, nonsyn
                        if source_aa == target_aa:
                            s_numer += self.B[codon[i] + n]
                        else:
                            n_numer += self.B[codon[i] + n]
                        denom += self.B[codon[i] + n]        

            # Tally if there were changes
            if s_numer > 0.:
                s_sites += s_numer/denom 
            if n_numer > 0:
                n_sites += n_numer/denom

        return n_sites, s_sites


    def count_sites_over_branch(self):
        '''
            Count number of expected number of nonsyn, syn sites for a given branch. 
        '''

        assert(len(self.path_codon_counts) > 0), "\n\nPath is empty."
        x = float(sum(self.path_codon_counts.values()))
        
        s_sites_raw = 0.
        n_sites_raw = 0.
        for key in self.path_codon_counts:
            ntemp, stemp = self.count_codon_sites(key)
            s_sites_raw += (self.path_codon_counts[key] * stemp)
            n_sites_raw += (self.path_codon_counts[key] * ntemp)
        self.s_sites = (s_sites_raw / x) * self.bl
        self.n_sites = (n_sites_raw / x) * self.bl








class dNdS_Counter(object):
    '''
        Class to count simulated, site-specific dN/dS.
        Required positional arguments:
            alnfile:  alignment file (fasta) produced by pyvolve *with ancestral sequences*
            treefile: file with newick tree 
            
        Optional positional arguments:
            mu_dict: dictionary of nucleotide mutation rates. If missing, assumes equal mutation rates.
    '''
    
    def __init__(self, alnfile, treefile, mu_dict = None):

        # Parse initial tree
        self.tree = newick.read_tree(file = treefile)
        
        # Read alignment, convert to dictionary, and determine nucleotide pi dictionary
        length = 0.
        nuc_counts = np.zeros(4)
        with open(alnfile, "rU") as handle:
            aln = AlignIO.read(handle, "fasta")
            self.alnlen = len(aln[0])/3
            self.alndict = {}
            for rec in aln:
                s = str(rec.seq).upper()
                nuc_counts[0] += s.count("A")
                nuc_counts[1] += s.count("C")
                nuc_counts[2] += s.count("G")
                nuc_counts[3] += s.count("T")
                self.alndict[str(rec.id)] = s
    
        # Build B.
        self.compute_B(mu_dict, nuc_counts)
        
        # How many branch lengths?
        self.num_edges = 0
        self.sum_bl = 0.
        self.tally_bl(self.tree)

        
        # Set up storage for raw site-wise dn, ds quantities.
        self.n_sites   = np.zeros( [self.alnlen, self.num_edges] ) 
        self.s_sites   = np.zeros( [self.alnlen, self.num_edges] ) 
        self.n_changes = np.zeros( [self.alnlen, self.num_edges] ) 
        self.s_changes = np.zeros( [self.alnlen, self.num_edges] ) 
       

    def compute_B(self, mu_dict, nuc_probs):
        '''
            Given nucleotide mutational information, return a dictionary 
            giving the relative probability of subbing nucleotide pairs.
        '''
        nuc_probs /= np.sum(nuc_probs)
        pi_dict = {"A": nuc_probs[0], "C": nuc_probs[1], "G": nuc_probs[2], "T": nuc_probs[3]} 
        
        # Default.
        if mu_dict is None:
            mu_dict = {"AC":1., "AG":1., "AT":1., "CG":1., "CT":1., "GT":1.}
        
        # Fill in missing mutation rates symmetrically
        temp_mu_dict = {}
        for key in mu_dict:
            if key[1] + key[0] not in mu_dict:
                temp_mu_dict[key[1] + key[0]] = mu_dict[key]
        
        # Build unnormalized B
        new = {}
        for key in mu_dict:
            new[key] = mu_dict[key] * pi_dict[key[1]]
        for key in temp_mu_dict:
            new[key] = temp_mu_dict[key] * pi_dict[key[1]]
    
        # Build normalized B
        k = new.keys()
        v = np.array(new.values())
        v /= np.sum(v)
    
        self.B = dict(zip(k,v)) 
        
         
    
    def tally_bl(self, t):
        '''
            Recursive function to retrieve the total number of edges and the sum of all branch lengths.
        ''' 
        if t.branch_length is not None:
            self.sum_bl += t.branch_length
            self.num_edges += 1
        if len(self.tree.children)>0:
            for child in t.children:
                self.tally_bl(child)    
        
    



    def traverse_tree_dnds(self, source_node, target_node, storage_index):
        '''
            Traverse the tree to compute and store dN, dS quantities (sites and changes) at each edge.
        '''
    
        # Compute site-wise dN/dS along this branch, where bl = target_node.branch_length
        try:
            full_target_seq = self.alndict[target_node.name]
        except KeyError:
            print("\n\nTree node names do not have matches in alignment file. Make sure that the provided alignment file *includes ancestral sequences*.")
            sys.exit() 
        full_source_seq = self.alndict[source_node.name]
        bl = target_node.branch_length

        for s in range(0, self.alnlen*3, 3):
            source_seq = full_source_seq[s:s+3]
            target_seq = full_target_seq[s:s+3]
            
            calcer = CalcOverBranch(source_seq, target_seq, bl, self.B)
            n_changes, s_changes, n_sites, s_sites= calcer.compute_branch_dnds()

            # Save quantities
            self.n_changes[s/3][storage_index] = n_changes
            self.s_changes[s/3][storage_index] = s_changes
            self.n_sites[s/3][storage_index] = n_sites
            self.s_sites[s/3][storage_index] = s_sites


        
        storage_index += 1
        # Proceed down the tree only if there are no more children
        #print "children:", target_node.children[0].name, target_node.children[1].name 
        if len(target_node.children) > 0:
            for child in target_node.children:
                storage_index = self.traverse_tree_dnds(target_node, child, storage_index) 
        return storage_index
       

    def calculate_dnds(self, savefile = "counted_dnds.txt"):
        '''
            Function *for users to call* in order to compute, save site-wise dnds quantities over a tree.
        '''
        
        # Obtain all quantities. We have to start with root's children, not root.
        storage_index = 0
        for subtree in self.tree.children:
            storage_index = self.traverse_tree_dnds(self.tree, subtree, storage_index)
            

        # Normalize sites by total branch length
        self.n_sites /= self.sum_bl
        self.s_sites /= self.sum_bl

        # Obtain per-site average of all quantities
        final_nsites   = np.mean(self.n_sites, axis=1)
        final_ssites   = np.mean(self.s_sites, axis=1)
        final_nchanges = np.mean(self.n_changes, axis=1)
        final_schanges = np.mean(self.s_changes, axis=1)
                        
        # Finally, save this csv: site_index, nonsyn_changes, syn_changes, nonsyn_sites, syn_sites
        with open(savefile, "w") as f:
            site = 1
            f.write("site\tns_changes\ts_changes\tns_sites\ts_sites")
            for i in range(len(final_nsites)):
                f.write("\n" + str(site) + "\t" + str(final_nchanges[i]) + "\t" + str(final_schanges[i]) + "\t" + str(final_nsites[i]) + "\t" + str(final_ssites[i]))
                site += 1

