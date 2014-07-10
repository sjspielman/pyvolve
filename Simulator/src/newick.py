from misc import Tree
import re
import os


def readTree(**kwargs):
    filename    = kwargs.get('file', None)
    tstring     = str(kwargs.get('tree', ''))
    return_model_flags = kwargs.get('flags', False)
        
    if filename:
        assert (os.path.exists(filename)), "File does not exist. Check path?"
        t = open(filename, 'r')
        tstring = t.read()
        t.close()
    else:
        assert (tstring != ''), "You need to either specify a file with a tree or give your own."
        
    tstring = re.sub(r"\s", "", tstring)
    tstring = tstring.rstrip(';')
    
    flags = []
    internal_node_count = 1
    (tree, flags, internal_node_count, index) = parse_tree(tstring, flags, internal_node_count, 0) 
    assert(flags == list(set(flags)) ), "Unique identifiers required for branch model heterogeneity flags."
    if return_model_flags:
        return tree, flags
    else:
        return tree


def read_model_flag(tstring, index):
    ''' Model flags are of the format _flag_ and come **after** the branch length associated with that node, before the comma.'''
    index +=1 # Skip the leading underscore
    end = index
    while True:
        end+=1
        if tstring[end]=='_':
            break
    model_flag = tstring[index:end]
    return model_flag, end+1
     
     
def read_branch_length(tstring, index):
    end = index
    while True:
        end += 1
        if end==len(tstring):
            break
        if tstring[end]==',' or tstring[end]==')' or tstring[end] == '_':
            break
    BL = float( tstring[index+1:end] )
    return BL, end


def read_leaf(tstring, index):
    end = index
    node = Tree()
    while True:
        end += 1
        assert( end<len(tstring) ), "Still trying to parse but have reached the end of your tree. Problematic."
        # Leaf has a branch length
        if tstring[end]==',' or tstring[end]==')':
            node.name = tstring[index+1:end]
            node.branch_length = None
            break    
        # Leaf has no branch length    
        if tstring[end]==':' :
            node.name = tstring[index:end]
            node.branch_length, end = read_branch_length(tstring, end)
            break        
    return node, end


def parse_tree(tstring, flags, internal_node_count, index):
    assert(tstring[index]=='(')
    index += 1
    node = Tree()
    while True:
        
        # New subtree (node) to parse
        if tstring[index]=='(':
            subtree, flags, internal_node_count, index=parse_tree(tstring, flags, internal_node_count, index)
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
                    BL, index = read_branch_length(tstring, index)
                    node.branch_length = BL
                if tstring[index]=='_':
                    model_flag, index = read_model_flag(tstring, index)
                    node.model_flag = model_flag
                    flags.append(model_flag)
            break
        # Terminal leaf
        else:
            subtree, index = read_leaf(tstring, index)
            node.children.append( subtree )
    return node, flags, internal_node_count, index    
    
    
def print_tree(tree, level=0):
    indent=''
    for i in range(level):
        indent+='\t'
    print indent, tree.name, tree.branch_length, tree.model_flag
    if len(tree.children)>0:
        for node in tree.children:
            print_tree(node, level+1)    
    
            
