#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class SiteCounter:

	def countSitesInCodon( self, codon ):
		"""Counts the number of synonymous and non-synonymous sites in the codon.

The function returns two lists. The first is the number of non-synonymous
sites at each of the three sites, and the second is the number of synonymous
sites at each of the three sites.
"""
		assert len( codon ) == 3
		
		s_counts = [0., 0., 0.]
		ns_counts = [0., 0., 0.]
		
		aa_source = str(codon.translate())
		assert aa_source != '*' # There is no clear way to handle sites that translate to stop codons. Those should be omitted from the sequence
		
		for i in range( 3 ):
			mutations = 0
			for c in ['A', 'C', 'T', 'G']:
				if codon[i].upper() == c:
					continue
				target = codon[0:i] + c + codon[i+1:3]
				
				aa_target = str(target.translate())
				if aa_target == '*': # disregard mutations to stop codons
					continue
				
				if aa_source == aa_target:
					s_counts[i] += 1
				else:
					ns_counts[i] += 1
				mutations += 1
			s_counts[i] /= mutations
			ns_counts[i] /= mutations
			assert mutations <= 3
		return ( ns_counts, s_counts )


	def countSites( self, seq ):
		"""Counts the number of synonymous and non-synonymous sites in a sequence. The function needs a sequence in which all codons are complete.

The function returns two lists. The first is the number of non-synonymous
mutations at each of the three sites, and the second is the number of synonymous
mutations at each of the three sites.
"""
		assert len( seq ) % 3 == 0
		
		ns_counts = []
		s_counts = []
		for i in range( len( seq )/ 3 ):
			c1 = seq[3*i:3*i+3]
#			print "Codon", c1
			( ns, s ) = SiteCounter.countSitesInCodon( self, c1 )
#			print "ns:", ns, "s:", s
			ns_counts += ns
			s_counts += s
		return ( ns_counts, s_counts )

#seq = Seq("AGTGCCTCAGCAAATTCA", IUPAC.unambiguous_dna)

#S = SiteCounter()
#( ns_counts, s_counts ) = S.countSites( seq )
#print "Sequence:", str(seq)
#print "    Synonymous sites:", s_counts
#print "Non-synonymous sites:", ns_counts
