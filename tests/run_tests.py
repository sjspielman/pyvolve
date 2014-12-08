#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

''' Run all unit tests in coordination with python setup.py test '''

import unittest
from models_test import *
from matrix_builder_test import *
from state_freqs_test import *
from evolver_test import *


if __name__ == '__main__':

    print "Running unit tests for pyvolve"

    print "\n\nRunning tests for state_freqs module"
    run_state_freqs_test()
    print "\n\nRunning tests for matrix_builder module"
    run_matrix_builder_test() 
    print "\n\nRunning tests for models module"
    run_models_test()
    print "\n\nRunning tests for evolver module"
    run_evolver_test()