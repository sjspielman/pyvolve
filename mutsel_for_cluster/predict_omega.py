from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from math import *



def calcAAFreqs( seq ):
    aa_seq = seq.translate()
    counts = {}
    for aa in aa_seq:
        if aa in counts:
            counts[aa] += 1
        else:
            counts[aa] = 1
    return counts

def calcCodonFreqs( seq ):
    counts = {}
    s = str(seq)
    for i in xrange( len(s)/3 ):
        codon = s[3*i:3*i+3]
        if codon in counts:
            counts[codon] += 1
        else:
            counts[codon] = 1
    return counts

def calcPfix( f1, f2 ):
    if f1 == f2:
        return 1.
    elif f1 == 0 or f2 == 0:
    	return 0.
    else:
        return (log(f2)-log(f1))/(1.-f1/f2)
        
def calcNonsynPaths( codon, codon_freqs ):
    cod = Seq( codon, IUPAC.unambiguous_dna)
    aa_source = str(cod.translate())
    assert aa_source != '*'
    f1 = codon_freqs[codon]
    rate = 0.
    sites = 0.
    
    for i in range( 3 ):
        for c in ['A', 'C', 'T', 'G']:
            if cod[i].upper() == c:
                continue
            target = cod[0:i] + c + cod[i+1:3]

            aa_target = str(target.translate())
#            print target, aa_target
            if aa_target == '*': # disregard mutations to stop codons
                continue

            if aa_source == aa_target:
                continue
            else:
                sites += 1
                if str(target) in codon_freqs:
                    rate += calcPfix( f1, codon_freqs[str(target)] )
    return ( rate, sites )

def predictdNdS( codon_freqs ):
    num = 0.
    denom = 0.
    rate_sum = 0.
    for codon in codon_freqs.keys():
        #print codon
        (rate, sites) = calcNonsynPaths( codon, codon_freqs )
        #print rate, sites, rate/sites
        num += codon_freqs[codon]*rate
        denom += codon_freqs[codon]*sites
        rate_sum += codon_freqs[codon]*rate/sites
    print "Predicted dN/dS, version 1:", round( num/denom, 5 )
    print "Predicted dN/dS, version 2:", round( rate_sum, 5 )