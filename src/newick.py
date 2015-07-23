#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Read/parse a newick tree.
'''

import re
import os

class Tree():
    '''
        Defines a Tree() object, which is ultimately comprised of a series of nested Tree() objects.
        In effect, each Tree() instance represents a node within a larger phylogeny.
    '''
    def __init__(self):
            
        self.name           = None # Internal node unique id or leaf name
        self.children       = []   # List of children, each of which is a Tree() object itself. If len(children) == 0, this tree is a tip.
        self.branch_length  = None # Branch length leading up to node
        self.model_flag     = None # Flag indicate that this branch evolves according to a distinct model from parent
        self.seq            = None # Contains sequence (represented by integers) for a given node. Later, this may instead be a list of Site objects.



def read_tree(**kwargs):
    
    '''
        Parse a newick phylogeny, provided either via a file or a string. The tree does not need to be bifurcating, and may be rooted or unrooted.
        Returns a Tree() object along which sequences may be evolved.  
            
        Trees can either read from a file or given directly to ``read_tree`` as a string. One of these two keyword arguments is required.
        
            1. **file**, the name of the file containing a newick tree for parsing. If this argument is provided in addition to tstring, the tree in the file will be used and tstring will be ignored.
            2. **tree**, a newick tree string. If a file is additionally provided, the tstring argument will be ignored.   
        
        To implement branch (temporal) heterogeneity, place "model flags" at particular nodes within the tree. Model flags must be in the format *_flagname_* (i.e. with both a leading and a trailing underscore), and they should be placed *after* the branch lengths.
        Model flags may be repeated throughout the tree, but the model associated with each model flag will always be the same. Note that these model flag names **must** have correspondingly named model objects.


        Examples:
            .. code-block:: python
                
               tree = read_tree(file = "/path/to/tree/file.tre")
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921):0.207);")
               
               # Tree containing model flags named m1 and m2
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921)_m2_:0.207);"
    '''    
    
    filename           = kwargs.get('file')
    tstring            = kwargs.get('tree')
        
    if filename:
        assert (os.path.exists(filename)), "File does not exist. Check path?"
        t = open(filename, 'r')
        tstring = t.read()
        t.close()
    else:
        assert (tstring is not None), "You need to either specify a file with a tree or give your own."
        assert (type(tstring) is str), "Trees provided with the flag `tree` must be in quotes to be considered a string."
    
    # Clean up the string a bit    
    tstring = re.sub(r"\s", "", tstring)
    tstring = re.sub(r":\d+\.*\d*;$", "", tstring) # In case there is a "root bl" at end of string. This mucks up parser.
    tstring = tstring.rstrip(';')

    flags = []
    internal_node_count = 1
    (tree, flags, internal_node_count, index) = _parse_tree(tstring, flags, internal_node_count, 0) 
    nroots = 0
    pf, nroots = _assign_model_flags_to_nodes(nroots, tree)
    assert(nroots == 1), "\n\nYour tree has not been properly specified. Please ensure that all internal nodes and leaves have explicit branch lengths (even if the branch lengths are 0)."
    return tree


def print_tree(tree, level=0):
    '''
        Prints a Tree() object in graphical, nested format. 
        This function takes two arguments:
            
            1. **tree** is a Tree object to print
            2. **level** is used internally for printing. DO NOT PROVIDE THIS ARGUMENT.
        
        Each node in the tree is represented by a string in the format, "name   branch.length   model.flag", and levels are represented by indentation.
        Names for tree tips are taken directly from the provided tree, and internal node names are assigned automatically by the ``read_tree`` function.
        The node with a branch length of None will be the root node where sequence evolution will begin.
        Note that the model.flag field will be None under cases of branch homogeneity.       
        
        For example,
            .. code-block:: python
            
               >>> my_tree = newick.read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921):0.207);")
               >>> print_tree(my_tree)
                    root None None
                        t4 0.785 None
                            internal_node3 0.207 None
                                t3 0.38 None
                                internal_node2 0.921 None
                                    t2 0.806 None
                                    internal_node1 0.762 None
                                        t5 0.612 None
                                        t1 0.66 None
            
               >>> flagged_tree = newick.read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921_m2_):0.207);")
               >>> newick.print_tree(flagged_tree)  
                     root None None
                        t4 0.785 None
                        internal_node3 0.207 None
                            t3 0.38 None
                            internal_node2 0.921 m2
                                t2 0.806 m2
                                internal_node1 0.762 m1
                                    t5 0.612 m1
                                    t1 0.66 m1

                            
    ''' 
    indent=''
    for i in range(level):
        indent+='\t'
    printstring = indent + str(tree.name) + " " + str(tree.branch_length) + " " + str(tree.model_flag)
    print(printstring)
    if len(tree.children) > 0:
        for node in tree.children:
            print_tree(node, level+1)    
    


def _assign_model_flags_to_nodes(nroots, tree, parent_flag = None):
    '''
        Determine the evolutionary model to be used at each node.
        Note that parent_flag = None means root model!!
    '''
    
    # Assign model if there was none in the tree    
    if tree.model_flag is None:
        tree.model_flag = parent_flag
        if tree.name == "root":
            nroots += 1

    if len(tree.children) > 0:
        for node in tree.children:
            parent_flag, nroots = _assign_model_flags_to_nodes(nroots, node, tree.model_flag)
    return parent_flag, nroots
    

def _read_model_flag(tstring, index):
    '''
        Read a model flag id while parsing the tree from the function _parse_tree.
        Model flags are expected to be in the format _flag_, and they must come **after** the branch length associated with that node, before the comma.
    '''
    index += 1 # Skip the leading underscore
    end = index
    while True:
        end+=1
        if tstring[end]=='_':
            break
    model_flag = tstring[index:end]
    return model_flag, end+1
     
     
def _read_node_name(tstring, index):
    '''
        Read a provided internal node name while parsing the tree from the function _parse_tree.
        Importantly, internal node names *MAY NOT* contain colons!!
    '''
    end = index
    while True:
        if end==len(tstring):
            break
        if tstring[end] == ":":
            break
        end += 1
    name = tstring[index:end]
    
    return name, end

     
def _read_branch_length(tstring, index):
    '''
        Read a branch length while parsing the tree from the function _parse_tree.
    '''
    end = index
    while True:
        end += 1
        if end==len(tstring):
            break
        if tstring[end]==',' or tstring[end]==')' or tstring[end] == '_':
            break
    BL = float( tstring[index+1:end] )
    
    return BL, end


def _read_leaf(tstring, index):
    '''
        Read a leaf (taxon name) while parsing the tree from the function _parse_tree.
    '''
    end = index
    node = Tree()
    while True:
        end += 1
        assert( end<len(tstring) ), "\n\nTree parsing error! Please ensure that your tree is a properly specified newick tree with branch lengths for all nodes and tips. Consult the Pyvolve manual for proper internal node name and model flag specification."
        # Leaf has no branch length
        if tstring[end]==',' or tstring[end]==')':
            node.name = tstring[index+1:end]
            node.branch_length = None
            break    
        # Leaf has branch length    
        if tstring[end]==':' :
            node.name = tstring[index:end]
            node.branch_length, end = _read_branch_length(tstring, end)
            break       
    # Does leaf have a model? 
    if tstring[end] == '_':
        node.model_flag, end = _read_model_flag(tstring, end)
    return node, end


def _parse_tree(tstring, flags, internal_node_count, index):
    '''
        Recursively parse a newick tree string and convert to a Tree object. 
        Uses the functions _read_branch_length(), _read_leaf(), _read_model_flag() during the recursion.
    '''
    assert(tstring[index]=='(')
    index += 1
    node = Tree()
    while True:
        
        # New subtree (node) to parse
        if tstring[index]=='(':
            subtree, flags, internal_node_count, index = _parse_tree(tstring, flags, internal_node_count, index)
            node.children.append( subtree )
        
        # March to sister
        elif tstring[index]==',':
            index += 1            
        
        # End of a subtree (node)
        elif tstring[index]==')':
            index+=1
            
            # Now we have either a node name, model flag, BL. Order should be node, BL, model flag (if/when multiple).            
            if index<len(tstring):
                if re.match(r"^[A-Za-z]", tstring[index]):
                    name, index = _read_node_name(tstring, index)
                    node.name = name   
                # Quick warning to prevent users from supply root names
                try:
                    blah = tstring[index]
                except:
                    raise IndexError("\n\nTree parsing error. This error probably occurred because you specified a name for the root node, which you can't do. Pyvolve must assign this node's name to 'root', by default.")
                if tstring[index]==':':
                    BL, index = _read_branch_length(tstring, index)
                    node.branch_length = BL
                if tstring[index]=='_':
                    model_flag, index = _read_model_flag(tstring, index)
                    node.model_flag = model_flag
                    flags.append(model_flag)

            # Assign name to the node, either as internal_code<i> or root (if the branch length is None), if a name was not specified.
            if node.name is None:
                # Root
                if node.branch_length is None:
                    node.name = "root"
                else:
                    # Unnamed internal node
                    node.name = "internal_node" + str(internal_node_count)
                    internal_node_count += 1
            
            # Check that branch lengths and node names were set up
            if node.name != "root":
                assert(node.branch_length is not None), "\nYour tree is missing branch length(s). Please ensure that all nodes and tips have a branch length (even if the branch length is 0!)."
            assert(node.name is not None), "\nInternal node name was neither provided nor assigned, which means your tree has not been properly formatted. Please ensure that you have provided a proper newick tree."
            
            break
            
        # Terminal leaf
        else:
            subtree, index = _read_leaf(tstring, index)
            node.children.append( subtree )
    return node, flags, internal_node_count, index    
    

            
