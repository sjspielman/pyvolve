import re

class Node:
	'''Node class. Each tree itself is a node containing subnodes all the way through.'''
	def __init__(self, name, lchild, rchild, BL):
		self.name = name # this can either be None (internal) or a leaf name.
		self.lchid = lchild # Left-hand child (node). None if leaf
		self.rchild = rchild # right-hand child (node). None if leaf.
		self.BL = BL # Branch length leading up to node


def readTree(treefile):
	''' Takes a file containing a fully bifurcating Newick tree and parses. Branch lengths required at every node.'''
	t = open(treefile, 'r')
	tstring = t.read()
	t.close()
	
	# Remove whitespace, exponent branch lengths, semicolon
	tstring = re.sub(r"\s", "", tstring) 
	tstring = re.sub("\d\.*\d*e-\d+", "0.00001", tstring) # If small enough to be exponent, this is a fine approx
	tstring = tstring.rstrip(';')
	#print tstring
	
	# Go to town
	tree = parseTree(tstring)
	return tree


def parseTree(tstring):
	''' Does the bulk of the tree parsing'''
	
	# Regex pattern subtree:branchlength
	nodePattern = re.compile(r'(?P<subtree>\(.+\))\:(?P<length>\d+\.*\d*)$')
	
	# Search for the pattern (subtree):branchlength
	findNode = nodePattern.search(tstring)
	if findNode:
		if findNode.group("subtree"): tstring = findNode.group("subtree")
		if findNode.group("length"): branchlength = float(findNode.group("length"))
		else: branchlength = 0
	
	# Now we have node's branch length and a tstring. The tstring
	
	
	
	parseSubTree(tstring, length)
	
	
	
	
	
	# Check if we are dealing with a terminal. This is the base condition
	if ',' not in tstring:
		print tstring
		name, branchlength = tstring.split( ":" )
		return( Node(name, None, None, branchlength) )
	

	
	# Grab left and right children
	(ltree, rtree) = parseSubTree(tstring, length)
	
	# Aaaaand keep going. If we're returning down here, it's an internal node so doesn't have a name.
	tree= Node(None, parseTree(left), parseTree(right), branchlength)
	return( tree )

def parseSubTree(subtree):
	'''parses a subtree to get the left and right hand children.'''

	expr = substring[1:-1] #get rid of surrounding parens to isolate our unbounded left,right subtrees
	paren_count = 0
	index = 0
	left='';right=''
	
	# Save left and right based on parentheses. When num open paren=num closed paren, we are done finding left. Remaining string is right hand side
	for i in expr:
		if i=='(':
			paren_count+=1	
		if i==')':
			paren_count-=1
		if paren_count==0 and i==',':
			left = expr[0:index]
			break
		index+=1
	right = expr[(index+1) : ] # use index+1 so we don't include the comma
	
	
	
	return ( Node(None, left, right, length) )
	











def tree2String(tree):

	newString = ""
	# if leaf, just go ahead and print
	if tree.name is not None:
		newString += tree.name
	else: #internal nodes will all have children
		newString += "("

		# descend left
		newString += tree2String(tree.lchild)

		# print comma
		newString += ","

		# descend right
		newString += tree2String(tree.rchild)
		
		# print paren
		newString += ")"

	# regardless of internal/leaf node, return length
	if tree.branchlength:
		newString += ":"
		newString += str(tree.branchlength)
	return newString
	


tree = readTree('small.tre')
#print tree2String(tree)
	
	
	