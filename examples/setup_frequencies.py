#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    This example script shows how to generate customized state frequencies (see documentation for more information: http://sjspielman.org/pyvolve/state_freqs.html )
    Generally, there are five classes for deriving state frequencies:
    1. EqualFrequencies
    2. RandomFrequencies
    3. CustomFrequencies
    4. ReadFrequencies
    5. EmpiricalModelFrequencies

    The first four compute custom frequencies, whereas the final option (EmpiricalModelFrequencies) will merely assign the default frequencies for a given empirical model (JTT, WAG, LG, or ECM).
    The remaining discussion here focuses on the general setup for each class.
    When initializing the frequency calculator, you must provide a single argument ("nuc", "amino", or "codon"). This argument indicates in which alphabet space frequencies will be *computed*.
    For instance, if you say EqualFrequencies("amino"), then a length 20 vector with each entry equal to 0.05 will be computed.
    To return this value, simply call the class instance: EqualFrequencies("amino")().
    Alternatively, you can also return a *converted* version of these frequencies by specifying the argument "type" when calling the instance. 
    In this way, you can convert amino acid frequencies to codon frequencies:  EqualFrequencies("amino")(type = "codon").
'''


from pyvolve import *

# Optional arguments to 