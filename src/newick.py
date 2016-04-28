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
import warnings 


MODEL_FLAGS = ("_", "#")

class Node():

    '''
        Defines a Node object, a list of nested Node objects. Recursive phylogeny subunit.
    '''
    def __init__(self):
            
        self.name            = None # Internal node unique id or leaf name
        self.children        = []   # List of children, each of which is a Node object itself. If len(children) == 0, this tree is a tip.
        self.branch_length   = None # Branch length leading up to node
        self.model_flag      = None # Flag indicate that this branch evolves according to a distinct model from parent
        self.propagate_model = True # Propagate model flag to the child nodes, default True
        self.seq             = None # Contains sequence (represented by integers) for a given node. Later, this may instead be a list of Site objects.


def read_tree(**kwargs):
    
    '''
        Parse a newick phylogeny, provided either via a file or a string. The tree does not need to be bifurcating, and may be rooted or unrooted.
        Returns a Node object along which sequences may be evolved.  
            
        Trees can either read from a file or given directly to ``read_tree`` as a string. One of these two keyword arguments is required.
        
            1. **file**, the name of the file containing a newick tree for parsing. If this argument is provided in addition to tstring, the tree in the file will be used and tstring will be ignored.
            2. **tree**, a newick tree string. If a file is additionally provided, the tstring argument will be ignored.   
        
        Optional keyword arguments:
            1. **scale_tree** is a float value for scaling all branch lengths by a given multiplier. Default: 1.
        
        To implement branch (temporal) heterogeneity, place "model flags" at particular nodes within the tree. Model flags can be specified with either underscores (_) or hashtags (#), through one of two paradigms:
            + Using trailing and leading symbols, e.g. _flagname_ or #flagname# . Specifying a model flag with this format will cause ALL descendents of that node to also follow this model, unless a new model flag is given downstream.
            + Using *only a leading* symbol, e.g. _flagname or #flagname. Specifying a model flag with this format will cause ONLY that branch/edge to use the provided model. Descendent nodes will NOT inherit this model flag. Useful for changing model along a single branch, or towards a single leaf.
        
        Model flags may be repeated throughout the tree, but the model associated with each model flag will always be the same. Note that these model flag names **must** have correspondingly named model objects.
        
        **IMPORTANT**: Node names must be provided BEFORE a branch length, and model flags be provided AFTER a branch length. For example, this subtree is correct: "...(taxon1:0.5, taxon2:0.2)<NODENAME>:<BL><MODEL FLAG>)...". This subtree is *incorrect* and will raise a cryptic error: "...(taxon1:0.5, taxon2:0.2):<BL><NODENAME><MODEL FLAG>)...". 


        Examples:
            .. code-block:: python
                
               tree = read_tree(file = "/path/to/tree/file.tre")
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921):0.207);")
               
               # Tree containing model flags named m1 and m2, both of which propagate to descendents.
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921)_m2_:0.207);"
               #or
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762#m1#):0.921)#m2#:0.207);"


               # Tree containing model flags named m1 and m2, each of which applies only to that branch.
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660_m1):0.762_m2):0.921):0.207);"
               #or
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660#m1):0.762#m2):0.921):0.207);"


               # Tree with a node demonstrating how to provide both a node name and model flag
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)NODENAME:0.762_m1_):0.921):0.207);" # propagating model flag
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921)NODENAME:0.207#m1);" # non-propagating model flag


               # Tree containing model flags named m1 and m2, where m1 is branch-specific but m2 is propagating.
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660#m1):0.762#m2#):0.921):0.207);" 
               #or
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660_m1):0.762_m2_):0.921):0.207);"
               #or
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660_m1):0.762#m2#):0.921):0.207);"
               #or
               tree = read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660#m1):0.762_m2_):0.921):0.207);"

    '''    
    
    ## Input arguments   
    filename    = kwargs.get('file')
    tstring     = kwargs.get('tree')
    scale_tree  = kwargs.get('scale_tree', 1.)
        
    ## Quick checks on input arguments
    if filename:
        assert (os.path.exists(filename)), "File does not exist. Check path?"
        t = open(filename, 'r')
        tstring = t.read()
        t.close()
    else:
        assert (tstring is not None), "\nYou need to either specify a file with a tree or give your own."
        assert (type(tstring) is str), "\nTrees provided with the flag `tree` must be in quotes to be considered a string."
    try:
        scale_tree = float(scale_tree)
    except:
        raise TypeError("\nThe argument 'scale_tree' must be a number (integer or float).")

    
    # Clean up the string a bit    
    tstring = re.sub(r"\s", "", tstring)
    tstring = re.sub(r":\d+\.*\d*;$", "", tstring) # In case there is a "root bl" at end of string. This mucks up parser.
    tstring = tstring.rstrip(';')

    flags = []
    internal_node_count = 1
    (tree, flags, internal_node_count, index) = _parse_tree(tstring, flags, internal_node_count, scale_tree, 0) 
    nroots = 0
    pf, nroots = _assign_model_flags_to_nodes(nroots, tree)
    assert(nroots == 1), "\n\nYour tree has not been properly specified. Please ensure that all internal nodes and leaves have explicit branch lengths (even if the branch lengths are 0)."
    return tree


def print_tree(tree, level=0):
    '''
        Prints a Node object in graphical, nested format. 
        This function takes two arguments:
            
            1. **tree** is a Node object to print
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
            # By default, model will progagate
            if tree.propagate_model:
                children_model_flag = tree.model_flag 
            else:
                children_model_flag = None
            parent_flag, nroots = _assign_model_flags_to_nodes(nroots, node, children_model_flag)
    return parent_flag, nroots
    
    
    

def _read_model_flag(tstring, index):
    '''
        Read a model flag id while parsing the tree from the function _parse_tree. Flags must come **after** the branch length associated with that node, before the comma.
        Model flags can be indicated with either underscores (_) or hash signs (#). There are two strategies:
            + Leading and trailing, e.g. #flag# or _flag_ . These flags will automatically propagate to all child branches.
            + Trailing only, e.g. #flag or _flag. These flags will be applied *only* to the given branch.
    '''
    
    flag_symbol = tstring[index]
    assert(flag_symbol in MODEL_FLAGS), "\nError: Unknown model flag."

    index += 1 # Skip the leading flag symbol
    end = index
    prop = True # detected flag is propagated, by default
    
    while True:
        end += 1
        if end==len(tstring) or tstring[end] == ":" or tstring[end] == ")" or tstring[end] == ",":
            break
        if tstring[end] == flag_symbol:
            end += 1
            break
    model_flag = tstring[index:end]
    
    # Clean model flag and determine if propagating
    if model_flag.endswith(flag_symbol):
        model_flag = model_flag[:-1]
    else:
        prop = False
    
    # If we had a propagating model, then increment end to remove the trailing symbol
    if tstring[end] == flag_symbol:
        assert(prop is True), "\n\nPyvolve can't tell if your model flag is propagating or not. Please consult docs."
        end += 1
    
    return model_flag, prop, flag_symbol, end

     

     

def _read_node_name(tstring, index):
    '''
        Read a provided internal node name while parsing the tree from the function _parse_tree.
        Importantly, internal node names *MAY NOT* contain colons!!
    '''
    end = index
    while True:
        if end == len(tstring):
            break
        if tstring[end] == ":" or tstring[end] in MODEL_FLAGS:
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
        if tstring[end]==',' or tstring[end]==')' or tstring[end] in MODEL_FLAGS:
            break
    BL = float( tstring[index+1:end] )
    return BL, end


def _read_leaf(tstring, index):
    '''
        Read a leaf (taxon name) while parsing the tree from the function _parse_tree.
    '''
    end = index
    node = Node()

    while True:
        end += 1
        assert( end<len(tstring) ), "\n\nTree parsing error! Please ensure that your tree is a properly specified newick tree with branch lengths for all nodes and tips. Consult the Pyvolve manual for proper internal node name and model flag specification."

        # Leaf has no branch length -> raise error
        if tstring[end]==',' or tstring[end]==')':
            raise AssertionError("\n\nThe leaves on your provided tree do not all have branch lengths. *All* branch lengths must be specified (even if they are 0!).")   
        
        # Leaf has branch length    
        if tstring[end] == ':' :
            node.name = tstring[index:end]
            node.branch_length, end = _read_branch_length(tstring, end)            
            break 
                
    # Does leaf have an associated model? 
    if tstring[end] in MODEL_FLAGS:
        node.model_flag, propagate, flag_symbol, end = _read_model_flag(tstring, end)
        
        # Clean leaf name as needed
        if flag_symbol in node.name:
            node.name, node.model_flag = node.name.split( flag_symbol, 1 )
    return node, end





def _parse_tree(tstring, flags, internal_node_count, scale_tree, index):
    '''
        Recursively parse a newick tree string and convert to a Node object. 
        Uses the functions _read_branch_length(), _read_leaf(), _read_model_flag() during the recursion.
    '''
    assert(tstring[index]=='(')
    index += 1
    node = Node()
    while True:
        
        # New subtree (node) to parse
        if tstring[index]=='(':
            subtree, flags, internal_node_count, index = _parse_tree(tstring, flags, internal_node_count, scale_tree, index)
            node.children.append( subtree )
        
        # March to sister
        elif tstring[index]==',':
            index += 1            
        
        # End of a subtree (node)
        elif tstring[index]==')':
            index += 1
        
            # Now we have either a node name, model flag, BL. Order MUST BE node name, BL, model flag (if/when multiple).            
            if index<len(tstring):
                        
                # Node name
                if re.match(r"^[A-Za-z]", tstring[index]): # Must start w/ letter
                    name, index = _read_node_name(tstring, index)
                    node.name = name
                                               
                # Branch length, with scaling as needed
                if tstring[index]==':':
                    BL, index = _read_branch_length(tstring, index)
                    node.branch_length = BL
                    
                # Model flag
                if tstring[index] in MODEL_FLAGS:
                    node.model_flag, node.propagate_model, flag_symbol, index = _read_model_flag(tstring, index)
                    flags.append( node.model_flag )                                    
            
            # Assign name to the node, either as internal_code<i> or root (if the branch length is None), if a name was not specified.
            if node.name is None:
                # Root
                if node.branch_length is None:
                    node.name = "root"
                else:
                    # Unnamed internal node
                    node.name = "internal_node" + str(internal_node_count)
                    internal_node_count += 1
                    node.branch_length *= scale_tree # scale *internal* branch length
            
            # Check that branch lengths and node names were set up
            if node.name != "root":
                assert(node.branch_length is not None), "\nYour tree is missing branch length(s). Please ensure that all nodes and tips have a branch length (even if the branch length is 0!)."
            assert(node.name is not None), "\nInternal node name was neither provided nor assigned, which means your tree has not been properly formatted. Please ensure that you have provided a proper newick tree."
            
            break
            
        # Terminal leaf
        else:
            subtree, index = _read_leaf(tstring, index)
            subtree.branch_length *= scale_tree # scale *leaf* branch length
            node.children.append( subtree )
    return node, flags, internal_node_count, index    
    

            