from misc import Tree
import re
import os

#internalCounter = 0



def readTree(**kwargs):
    filename    = kwargs.get('file', None)
    tstring     = str(kwargs.get('tree', ''))
    returnFlags = kwargs.get('flags', False)
        
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
    internalNodeCount = 0
    (tree, flags, internalNodeCount, index) = parseTree(tstring, flags, internalNodeCount, 0) 
    assert(flags == list(set(flags)) ), "Unique identifiers required for branch model heterogeneity flags."
    if returnFlags:
        return tree, flags
    else:
        return tree


def readModelFlag(tstring, index):
    ''' Model flags are of the format _flag_ and come **after** the branch length associated with that node, before the comma.'''
    index +=1 # Skip the leading underscore
    end = index
    while True:
        end+=1
        if tstring[end]=='_':
            break
    modelFlag = tstring[index:end]
    return modelFlag, end+1
     
     
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


def parseTree(tstring, flags, internalNodeCount, index):
    assert(tstring[index]=='(')
    index += 1
    node = Tree()
    while True:
        
        # New subtree (node) to parse
        if tstring[index]=='(':
            subtree, flags, internalNodeCount, index=parseTree(tstring, flags, internalNodeCount, index)
            node.children.append( subtree ) 
             
        
        # March to sister
        elif tstring[index]==',':
            index += 1            
        
        # End of a subtree (node)
        elif tstring[index]==')':
            internalNodeCount += 1
            index+=1
            # Now we have either a model flag, BL or both. But the BL will be *first*.            
            if index<len(tstring):
                if tstring[index]==':':
                    BL, index = readBranchLength(tstring, index)
                    node.branch = BL
                if tstring[index]=='_':
                    modelFlag, index = readModelFlag(tstring, index)
                    node.modelFlag = modelFlag
                    flags.append(modelFlag)
            node.name = internalNodeCount
            break
        # Terminal leaf
        else:
            subtree, index = readLeaf(tstring, index)
            node.children.append( subtree )
    return node, flags, internalNodeCount, index    
    
    
def printTree(tree, level=0):
    indent=''
    for i in range(level):
        indent+='\t'
    print indent, tree.name, tree.branch, tree.modelFlag, tree.seq
    if len(tree.children)>0:
        for node in tree.children:
            printTree(node, level+1)    
    
            
