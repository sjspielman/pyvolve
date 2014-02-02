## Figure out how to use dendropy using the rcoal(10) in small.tre. 2/2/14.

from ete2 import Tree

t1 = Tree("small.tre")

for node in t1.traverse("preorder"):
	print node.name