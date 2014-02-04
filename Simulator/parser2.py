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
	
	############ Format tree appropriately. ############### 
	tstring = re.sub(r"\s", "", tstring)                  # Remove all whitespace
	tstring = re.sub("\d\.*\d*e-\d+", "0.00001", tstring) # If small enough to be exponent, this is a fine approx
	tstring = re.sub(r"\:\.", ":0.", tstring)             # Be sure any branch length that is a decimal has at least one digit before decimal (as in, 0.9 instead of .9)
	tstring = re.sub(r"\)[^:]", "):0", tstring)           # Be sure all nodes have a branch length
	tstring = tstring.rstrip(';')                         # Remove trailing semi-colon
	
	# Go to town
	tree = parseTree(tstring)
	return tree






def parseTree(tstring):
	''' Does the bulk of the tree parsing'''
	
	stopChar = re.compile('[,\)]') # 
	
	## If we have an open parenthesis, there is a subtree to be read
	if tstring[0]=='(':
		(left_string, right_string) = subParse(tstring)
		BL = grabBL(tstring)
		return( Node(None, parseTree(left_string), parseTree(right_string), BL) ) ## I'm missing a step of the logic here, right???
	
	# If no subtree, we've arrived at a leaf
	else:
		name, BL = tstring.split( ":" )
		return( Node(name, None, None, BL) )
	
	
	
	
def subParse(tstring):
	tstring = tstring.lstrip('(') # Remove leading parenthesis
	
	#March down the tstring. When we find (==) character is a comma, we've reached a split
	paren_count=0
	for index in len(tstring):
		if tstring[index]=='(':	paren_count+=1
		if tstring[index]==')':	paren_count-=1
		if parent_count==0 and tstring[index]==',':
			left_string = tstring[ 0 : index ]
			right_string = tstring[ (index+1) : ] # be sure not to include the comma
	return(left_string, right_string)	
	

def grabBL(tstring):
	
	stopChar = re.compile(r'?P<branchlen>(:\d+\.*\d*)[,\)]')
	findLength = stopChar.search(tstring)
	if findLength.group("branchlen"):	BL = float(findLength.group("branchlen"))
	
	return BL
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	