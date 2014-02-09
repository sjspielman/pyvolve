from misc import Tree
import re

def readTree(**kwargs):
	filename = kwargs.get('file', 'tre.tre')
	show = kwargs.get('show', False)
	assert os.path.exists(filename), "Tree file does not exist. Try again."
	
	t = open(filename, 'r')
	tstring = t.read()
	t.close()
	tstring = re.sub(r"\s", "", tstring)
	tstring = tstring.rstrip(';')
	(tree, index) = parseTree(tstring,  0)
	if show:
		printTree(tree)
	return tree
	 
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
	node = Tree()
	while True:
		end += 1
		assert( end<len(tstring) )
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
	
def parseTree(tstring, index=0):
	assert(tstring[index]=='(')
	index += 1
	node = Tree()
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
				node.branch = BL
			break
		else:
			subtree, index = readLeaf(tstring, index)
			node.children.append( subtree )
	return node, index		
	
	
def printTree(tree, level=0):
	indent=''
	for i in range(level):
		indent+='\t'
	print indent, tree.name, tree.branch, tree.seq
	if len(tree.children)>0:
		for node in tree.children:
			print tree.seq
			printTree(node, level+1)	