import re

class Node:
	def __init__(self):
		self.isBase = False # Assigned to True if base of tree. Note that this should have name=None and BL=None. This will serve as an extra check.
		self.name = None # this can either be None (internal) or a leaf name.
		self.children = [] # list of children, each of which is a node
		self.BL = None # Branch length leading up to node
		self.seq = "" # Sequence can be stored here when simulating


def printTree(tree, level=0):
	indent=''
	for i in range(level):
		indent+='\t'
	print indent, tree.name, tree.BL
	if len(tree.children)>0:
		for node in tree.children:
			printTree(node, level+1)		


def readBranchLength(tstring, index):
	assert(tstring[index]==':')
	
	end = index
	while True:
		end += 1
		if end==len(tstring):
			break
		if tstring[end]==',' or tstring[end]==')':
			break
	BL = float( tstring[index+1:end] )
	return BL, end

def readLeaf(tstring, index):
	end = index
	node = Node()
	while True:
		end += 1
		assert( end<len(tstring) )
		
		# Leaf has a branch length
		if tstring[end]==',' or tstring[end]==')':
			node.name = tstring[index+1:end]
			node.BL = None
			break
		
		# Leaf has no branch length	
		if tstring[end]==':' :
			node.name = tstring[index:end]
			node.BL, end = readBranchLength(tstring, end)
			break
			
	return node, end



def parseTree(tstring, index=0):

	assert(tstring[index]=='(')

	index += 1
	node = Node()
	node.isBase = True
	while True:
	
		if tstring[index]=='(':
			subtree, index=parseTree(tstring, index)
			node.children.append( subtree )
		
		elif tstring[index]==',':
			index += 1
			
		elif tstring[index]==')':
			index+=1
			if index<len(tstring):
				BL, index = readBranchLength(tstring, index)
				node.BL = BL
			break
		
		else:
			subtree, index = readLeaf(tstring, index)
			node.children.append( subtree )
	
	return node, index


def readTree(treefile):

	t = open(treefile, 'r')
	tstring = t.read()
	t.close()
	
	############ Format tree appropriately. ############### 
	tstring = re.sub(r"\s", "", tstring)                  # Remove all whitespace
	tstring = tstring.rstrip(';')                         # Remove trailing semi-colon

	# Go to town
	(tree, index) = parseTree(tstring,  0)
	return tree

	

tree = readTree('small.tre')
printTree(tree)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	