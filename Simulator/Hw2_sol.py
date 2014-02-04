#!/usr/bin/python

import re
import math


def printTree(tree):
#Takes a tree structure, and outputs in Newick format.
#Uses tree2String to convert to Newick string.
        treeStr = tree2String(tree)
        print treeStr + ";" #recall that the + operator concatenates strings


def tree2String(tree):

	newString = ""

	# if leaf, just go ahead and print
	if tree['name'] is not 'internal':
		newString += tree['name']
	else:
		# print paren
		newString += "("

		# descend left
		newString += tree2String(tree["left"])

		# print comma
		newString += ","

		# descend right
		newString += tree2String(tree["right"])
		
		# print paren
		newString += ")"

	# regardless of internal/leaf node, return length
	if tree['length']:
		newString += ":"
		newString += str(tree['length'])
	return newString


def loadTree(filename):
#This function reads a file in Newick format and returns our simple
#dictionary-based data structure for trees.
#Uses parseTree() to interpret input string.
	f = open(filename,'r')
	exp = f.read()
	f.close

	exp = exp.replace(';','') #ignore trailing (or other) semi-colons
	exp = re.sub(r'\s+','',exp) #ignore whitespace
	exp = re.sub(r'\n','',exp)
	exp = re.sub(r'\[.*\]','',exp) #ignore bracketed clauses

	len_present = exp.count(':')
	if len_present == 0:
		return parseNoLengthTree(exp)
	else:
		return parseLengthTree(exp)

def makeLeaf(name,length):
#This function returns a tree structure corresponding to a single leaf
	return { 'left':None, 'right':None, 'name':name, 'length':length }


def parseLengthTree(exp):
#This function takes a string in Newick format and parses it recursively.
#Each clause is expected to be of the general form (a:x,b:y):z
#where a and b may be subtrees in the same format.

	if ',' not in exp: #if this is a leaf
		name, length = exp.split(':')
		length = float(length)
		print "make leaf", exp
		return makeLeaf(name,length)

	#uses the regular expression features of Python
	distPattern = re.compile(r'(?P<tree>\(.+\))\:(?P<length>[e\-\d\.]+)$')
	m = distPattern.search(exp)
	print "full: ",exp
	length = 0
	if m:	
		print "yep"		
		if m.group('length'): length = float( m.group('length') )
		exp = m.group('tree')
		print "tree: ",exp
		print "length: ",length
		print '\n'
	if length == '': length = 0

	#Use the parseExp function to return the left and right hand sides
	#of the expression (e.g., a & b from (a,b))
	lhs, rhs = parseExp(exp)
	print lhs
	print rhs

	#Now package into a tree data structure
	return { "name":"internal",
			 "left":parseLengthTree(lhs), #recursively set up subtrees
			 "right":parseLengthTree(rhs),
			 "length":length }


def parseExp(exp):
	#Parse expression of type "a,b" into a & b where a and b can be
	#Newick formatted strings.
	
	print "exping"
	print exp
	print exp[1:-1]
	chars = list(exp[1:-1]) #get rid of surrounding parens, and split to list
	print chars
	print '\n'
	count = 0
	lhs = True #boolean to distinguish left and right side of the comma
	left = '' #store the left substring
	right = '' #store the right substring

	#a little tricky to deal with nested parens correctly
	for c in chars:
		if c == '(':
			count += 1
		if c == ')':
			count -= 1
		if (c == ',') and (count == 0) and (lhs) :
			lhs = False
			continue

		if lhs: left += c
		else: right += c

	#Now return the left and right substrings
	print "left:",left
	print "right:",right
	return [ left, right ]


tree = loadTree('small.tre')
printTree(tree)
