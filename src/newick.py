#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Read/parse a newick tree.
'''


from misc import Tree
import re
import os


def read_tree(**kwargs):
    
    '''
    Parses a newick phylogeny, provided either via a file or a string.
    Input tree may contain internal model flags at bifurcations (or trifurcations, polytomies are acceptable!), indicating that the daughter branches should evolve according to a different evolutionary model from parent. 
    Returns a Tree() object along which sequences may be evolved and a list of model flags,  
    
    Example input tree containing model flags, _m1_ and _m2_. Flags *must* be provided in format _flagname_ , i.e. with both a leading a and trailing underscore).
    ((((t1:1.0,t8:1.0):1.0,t7:1.0):1.0_m1_,((t2:1.0,t9:1.0):1.0,t3:1.0):1.0):1.0,(((t6:1.0,t4:1.0):1.0,t5:1.0):1.0_m2_,t10:1.0):1.0);

    ONE OF THESE TWO ARGUMENTS IS REQUIRED:
        1. *file* is the name of the file containing a newick tree for parsing. If this argument is provided in addition to tstring, the tree in the file will be used and tstring will be ignored.
        2. *tstring* is a newick string for a tree. Tree may be rooted or unrooted. If a file is additionally provided, the tstring argument will be ignored.   
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
        assert (type(tstring) is str), "Trees provided with the flag `tstring` must be in quotes to be considered a string."
        
    tstring = re.sub(r"\s", "", tstring)
    tstring = tstring.rstrip(';')
    
    flags = []
    internal_node_count = 1
    (tree, flags, internal_node_count, index) = _parse_tree(tstring, flags, internal_node_count, 0) 
    assert(flags == list(set(flags)) ), "Unique identifiers required for branch model heterogeneity flags."

    return tree, flags



def print_tree(tree, level=0):
    '''
    Prints a Tree() object to stdout in nested format. Mostly used for debugging purposes and/or visualization of tree structure.
    
    Arguments:
    `tree` is a Tree() object. 
    
    `level` is an internal argument used in the recursive printing strategy. Don't provide this argument at all, or the tree will not print properly.
    
    ''' 
    indent=''
    for i in range(level):
        indent+='\t'
    print indent, tree.name, tree.branch_length, tree.model_flag
    if len(tree.children)>0:
        for node in tree.children:
            print_tree(node, level+1)    
    
    
    

def _read_model_flag(tstring, index):
    '''
    Read a model flag id while parsing the tree from the function _parse_tree.
    Model flags are expected to be in the format _flag_, and they must come **after** the branch length associated with that node, before the comma.
    '''
    index +=1 # Skip the leading underscore
    end = index
    while True:
        end+=1
        if tstring[end]=='_':
            break
    model_flag = tstring[index:end]
    return model_flag, end+1
     
     
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
        assert( end<len(tstring) ), "\n\nUh-oh! I seem to have reached the end of the tree, but I'm still trying to parse something. Please check that your tree is in proper newick format."
        # Leaf has a branch length
        if tstring[end]==',' or tstring[end]==')':
            node.name = tstring[index+1:end]
            node.branch_length = None
            break    
        # Leaf has no branch length    
        if tstring[end]==':' :
            node.name = tstring[index:end]
            node.branch_length, end = _read_branch_length(tstring, end)
            break        
    return node, end


def _parse_tree(tstring, flags, internal_node_count, index):
    '''
    Recursively parse a newick tree string and convert to a Tree() object. 
    Uses the functions _read_branch_length(), _read_leaf(), _read_model_flag() during the recursion.
    '''
    assert(tstring[index]=='(')
    index += 1
    node = Tree()
    while True:
        
        # New subtree (node) to parse
        if tstring[index]=='(':
            subtree, flags, internal_node_count, index=_parse_tree(tstring, flags, internal_node_count, index)
            node.children.append( subtree ) 
             
        
        # March to sister
        elif tstring[index]==',':
            index += 1            
        
        # End of a subtree (node)
        elif tstring[index]==')':
            index+=1
            node.name = internal_node_count
            internal_node_count += 1
            # Now we have either a model flag, BL or both. But the BL will be *first*.            
            if index<len(tstring):
                if tstring[index]==':':
                    BL, index = _read_branch_length(tstring, index)
                    node.branch_length = BL
                if tstring[index]=='_':
                    model_flag, index = _read_model_flag(tstring, index)
                    node.model_flag = model_flag
                    flags.append(model_flag)
            break
        # Terminal leaf
        else:
            subtree, index = _read_leaf(tstring, index)
            node.children.append( subtree )
    return node, flags, internal_node_count, index    
    

            
