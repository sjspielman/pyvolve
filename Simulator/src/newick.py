from misc import Tree
import re
import os

def readTree(**kwargs):
    filename = kwargs.get('file', None)
    tstring = str(kwargs.get('tree', ''))
        
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
    (tree, flags, index) = parseTree(tstring,  flags, 0)
    #checkTree(tree)
    
    return tree


def readModelFlag(tstring, index):
    ''' Model flags are of the format _flag_ and come **after** the branch length associated with that node, before the comma.'''
    index +=1 # Skip the leading underscore
    end = index
    while True:
        end+=1
        if tstring[end]=='_':
            break
    model_flag = tstring[index:end]
    return model_flag, end+1
     
     
def readBranchLength(tstring, index):
    end = index
    while True:
        end += 1
        if end==len(tstring):
            break
        if tstring[end]==',' or tstring[end]==')' or tstring[end] == '_':
            break
    BL = float( tstring[index+1:end] )
    return BL, end


def readLeaf(tstring, index):
    end = index
    node = Tree()
    while True:
        end += 1
        assert( end<len(tstring) ), "Still trying to parse but have reached the end of your tree. Problematic."
        # Leaf has a branch length
        if tstring[end]==',' or tstring[end]==')':
            node.name = tstring[index+1:end]
            node.branch = None
            break    
        # Leaf has no branch length    
        if tstring[end]==':' :
            node.name = tstring[index:end]
            node.branch, end = readBranchLength(tstring, end)
            break        
    return node, end
    
def parseTree(tstring, flags, index=0):
    assert(tstring[index]=='(')
    index += 1
    node = Tree()
    while True:    
        if tstring[index]=='(':
            subtree, flags, index=parseTree(tstring, flags, index)
            node.children.append( subtree )        
        elif tstring[index]==',':
            index += 1            
        elif tstring[index]==')':
            index+=1
            # Now we have either a model flag or a BL            
            if index<len(tstring):
                if tstring[index]==':':
                    BL, index = readBranchLength(tstring, index)
                    node.branch = BL
                if tstring[index]=='_':
                    model_flag, index = readModelFlag(tstring, index)
                    node.model_flag = model_flag
                    flags.append(model_flag)
            break
        else:
            subtree, index = readLeaf(tstring, index)
            node.children.append( subtree )
    return node, flags, index    
    
    
def printTree(tree, level=0):
    indent=''
    for i in range(level):
        indent+='\t'
    print indent, tree.name, tree.branch, tree.model_flag, tree.seq
    if len(tree.children)>0:
        for node in tree.children:
            print tree.seq
            printTree(node, level+1)    
    
            
