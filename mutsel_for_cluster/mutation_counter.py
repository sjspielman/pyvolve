#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class MutationCounter:
	
	def evaluateMutationPath( self, path ):
		"""Helper function, called by enumerateMutationPaths().
Should never be called directly.
"""
		
#		print "  Evaluating path", path
		 # we set up temporary variables to store source and target codons for mutations
		self.source = self.codon1
		
		# now loop over all mutations on the path
		s_counts = [ 0, 0, 0 ]
		ns_counts = [ 0, 0, 0 ]
		has_stop = False
		for m in path:
			# generate the target codon from source and mutation
			self.target = self.source[0:m] + self.codon2[m] + self.source[m+1:3]
			# assess whether mutation is synonymous or not
			# these are the only two lines of code in this class that depends on Biopython. You can replace them with your own favorite genetic code table
			aa_source = str(self.source.translate())
			aa_target = str(self.target.translate())
			if aa_source == aa_target:
				s_counts[m] += 1
			else:
				ns_counts[m] += 1
			if aa_target == "*":
				has_stop = True
			
#			print "    Mutation from codon", str( self.source ), "to codon", str( self.target )
#			print "    Counts:", s_counts, ns_counts
		
			# move target codon to source codon
			self.source = self.target
			
		# we disregard all paths with a stop codon
		if not has_stop:
			# add the mutations from the latest path
			self.ns_counts = [sum(pair) for pair in zip(self.ns_counts, ns_counts)] 
			self.s_counts = [sum(pair) for pair in zip(self.s_counts, s_counts)] 
			self.path_count += 1


	def enumerateMutationPaths( self, mutation_list, path=[] ):
		"""Helper function, called by countMutationsInCodons().
Should never be called directly.
"""
		# we take each mutation in the list and then enumerate the remaining
		# mutations recursively
		n = len( mutation_list )
		if n == 0:
			# evaluate for all mutations on the path whether they are synonymous or not
			MutationCounter.evaluateMutationPath( self, path )
		for i in range( n ):
			m = mutation_list[i]
			new_list = mutation_list[0:i] + mutation_list[i+1:n]
			MutationCounter.enumerateMutationPaths( self, new_list, path + [m] )


	def countMutationsInCodons( self, codon1, codon2 ):
		"""Counts the number of synonymous and non-synonymous mutations
going from codon1 to codon2.

The function returns two lists. The first is the number of non-synonymous
mutations at each of the three sites, and the second is the number of synonymous
mutations at each of the three sites.
"""
		assert len( codon1 ) == 3
		assert len( codon2 ) == 3
		self.codon1 = codon1
		self.codon2 = codon2
	
		# find the sites at which there are mutations
		mutation_list = []
		for i in range( 3 ):
			if codon1[i] != codon2[i]:
				mutation_list.append( i )
	
		# now try all possible combinations of mutations
		self.s_counts = [0., 0., 0.]
		self.ns_counts = [0., 0., 0.]
		self.path_count = 0
		MutationCounter.enumerateMutationPaths( self, mutation_list )

		for i in range( 3 ):
			self.s_counts[i] /= self.path_count
			self.ns_counts[i] /= self.path_count
		return ( self.ns_counts, self.s_counts )


	def countMutations( self, seq1, seq2 ):
		"""Counts the number of synonymous and non-synonymous mutations
going from sequence1 to sequence2. The function needs two sequences of the same length and in which all codons are complete.

The function returns two lists. The first is the number of non-synonymous
mutations at each of the three sites, and the second is the number of synonymous
mutations at each of the three sites.
"""
		assert len( seq1 ) == len( seq2 )
		assert len( seq1 ) % 3 == 0
		assert len( seq2 ) % 3 == 0
		
#		print str(seq1)
#		print str(seq2)
		
		ns_counts = []
		s_counts = []
		for i in range( len( seq1 )/ 3 ):
			c1 = seq1[3*i:3*i+3]
			c2 = seq2[3*i:3*i+3]
#			print "Codon", c1, "to", c2
			( ns, s ) = MutationCounter.countMutationsInCodons( self, c1, c2 )
#			print "ns:", ns, "s:", s
			ns_counts += ns
			s_counts += s
		return ( ns_counts, s_counts )

#c1 = Seq("AGTGCCTGAGCAAATTCA", IUPAC.unambiguous_dna)
#c2 = Seq("GCTACCTACCGCAAATTC", IUPAC.unambiguous_dna)

#M = MutationCounter()
#( ns_counts, s_counts ) = M.countMutations( c1, c2 )
#print "    Synonymous mutations:", s_counts, sum( s_counts )
#print "Non-synonymous mutations:", ns_counts, sum( ns_counts )
