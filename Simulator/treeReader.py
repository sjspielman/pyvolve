from misc import Tree
import re

class Newick():
	def __init__(self, filename):
		self._FILE = filename

		
	def readTree(self, reveal=False):
		t = open(self._FILE, 'r')
		tstring = t.read()
		t.close()
		tstring = re.sub(r"\s", "", tstring)
		tstring = tstring.rstrip(';')
		(tree, index) = self.parseTree(tstring,  0)
		if reveal:
			self.printTree(tree)
		return tree
	 
	def printTree(self, tree, level=0):
		indent=''
		for i in range(level):
			indent+='\t'
		print indent, tree.name, tree.branch, tree.seq
		if len(tree.children)>0:
			for node in tree.children:
				self.printTree(node, level+1)		
	

	def readBranchLength(self, tstring, index):
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
	
	
	
	def readLeaf(self, tstring, index):
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
				node.branch, end = self.readBranchLength(tstring, end)
				break		
		return node, end
		
	
	
	def parseTree(self, tstring, index=0):
		assert(tstring[index]=='(')
		index += 1
		node = Tree()
		while True:	
			if tstring[index]=='(':
				subtree, index=self.parseTree(tstring, index)
				node.children.append( subtree )		
			elif tstring[index]==',':
				index += 1			
			elif tstring[index]==')':
				index+=1
				if index<len(tstring):
					BL, index = self.readBranchLength(tstring, index)
					node.branch = BL
				break
			else:
				subtree, index = self.readLeaf(tstring, index)
				node.children.append( subtree )
		return node, index
		
		