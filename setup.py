#!/usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
Setup.py script (uses setuptools) for building, testing, and installing pyvolve.

To build and install the package as root (globally), enter (from this directory!) - 
    sudo python setup.py build
    sudo python setup.py test   # OPTIONAL BUT RECOMMENDED. Please contact author with any failed tests! Note that every once in a while the tests for functions which must generate random numbers take excessively long, so just ctrl-C these and run tests again.
    sudo python setup.py install
    

To install for a particular user (locally), enter - 
    python setup.py build
    python setup.py test   # OPTIONAL BUT RECOMMENDED. Please contact author with any failed tests! Note that every once in a while the tests for functions which must generate random numbers take excessively long, so just ctrl-C these and run tests again.
    python setup.py build --user # where user is the computer account to install pyolve in
'''

_VERSION="0.8.2"

from setuptools import setup
setup(name = 'Pyvolve', 
    version = _VERSION, 
    description = 'Sequence simulation along phylogenies according to continuous-time Markov models',
    author = 'Stephanie J. Spielman', 
    author_email = 'stephanie.spielman@gmail.com', 
    url = 'https://github.com/sjspielman/pyvolve',
    download_url = 'https://github.com/sjspielman/pyvolve/tarball/' + _VERSION,
    platforms = 'Tested on Mac OS X.',
    package_dir = {'pyvolve':'src'},
    packages = ['pyvolve', 'tests'],
    package_data = {'tests': ['freqFiles/*', 'evolFiles/*']},
    install_requires=['numpy>=1.7', 'scipy', 'Biopython'],
    test_suite = "tests"
)
