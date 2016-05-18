"""``pyvolve`` package for simulating evolutionary sequences along a phylogeny.
Written by Stephanie J. Spielman.

Python modules
----------------
The package consists of the following Python modules:

* genetics

* model

* parameters_sanity

* newick

* partition

* state_freqs

* matrix_builder

* empirical_matrices

* evolver


"""
__version__ = '0.8.2'
from .model import *
from .newick import *
from .evolver import *
from .genetics import *
from .partition import *
from .state_freqs import *
from .matrix_builder import *
from .parameters_sanity import *
from .empirical_matrices import *


