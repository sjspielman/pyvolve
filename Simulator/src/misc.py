# Global variables to be used by all. 
ZERO = 1e-10

class Genetics():
    def __init__(self):
        self.pyrims       = ["C", "T"]
        self.purines      = ["A", "G"]
        self.nucleotides  = ["A", "C", "G", "T"]
        self.amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        self.genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
        self.codon_dict   = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
        self.codons       = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.stop_codons  = ["TAA", "TAG", "TGA"]

class Tree():
    def __init__(self):
        self.name             = None # internal node unique id or leaf name. in future, may change leaves to ints internally.
        self.children         = []   # list of children, each of which is a node
        self.branch_length    = None # Branch length leading up to node
        self.model_flag       = None # Flag for branch heterogeneity
        self.seq              = None # Will be a list of of lists. Outer lists are partitions. Each partition is then a list of Site instances (see below for def).


class Model():
    def __init__(self):
        self.subst_params = {} # parameters pertaining to substitution process
        self.indel_params = {} # parameters pertaining to indel process
        self.Q           = None

class Site():
    def __init__(self):
        self.int_seq    = None # Stores the integer value of the sequence at a site. Note that gaps (of any kind) are -1
        self.letter_seq = None # Stores the real sequence value of a site (letter or gap)
        self.origin    = None # Stores the name of the node of origin for this position's existence. Note that (7/5/14) root is currently the *largest* node value.
        self.state     = None # Stores state of sequence. 0=core, 1=insertion, 2=deleted core, 3=deleted insertion. No need to store a regular deletion, since those are effectively still root.
        
        
        
        
        
        
        
        