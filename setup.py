#!/usr/bin/env python

from setuptools import setup



# Setup
setup(name = 'pyvolve', 
    version = '0.1', 
    description = 'Evolutionary sequence simulation along a phylogeny',
    author = 'Stephanie J. Spielman', 
    author_email = 'stephanie.spielman@gmail.com', 
    url = 'https://github.com/sjspielman/pyvolve', 
    platforms = 'Tested on Mac OS X.',
    package_dir = {'pyvolve':'src'},
    packages = ['pyvolve', 'tests'],
    install_requires=['numpy', 'scipy', 'Biopython'],
    test_suite = "tests.matrix_builder_test"
)