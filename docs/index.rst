.. pyvolve documentation master file, created by
   sphinx-quickstart on Mon Jan 19 10:26:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyvolve's documentation!
===================================

pyvolve is a python library for simulating sequences along a phylogeny using continuous-time Markov models. pyvolve implements a wide variety of [standard modeling frameworks](http://www.molecularevolution.org/resources/models/):

1. Nucleotide models (GTR and all nested variants)
2. Empirical amino acid models (JTT, WAG, and LG)
3. Codon (dN/dS) models (GY94-style and MG94-style)
4. The ECM empirical codon model
5. Mutation-selection models (Halpern and Bruno, 1998), implemented for both nucleotides and codons

In particular, the pyvolve framework allows you to heavily customize your model, both in terms of stationary frequency parameters and rate matrix elements. 

Download and Installation
-----------

pyvolve is freely available for download and use from https://github.com/sjspielman/pyvolve. Select the latest release (v0.1) from the "releases" tab to download. Once pyvolve has been downloaded, navigate to the directory in the terminal. To install for all users, enter these commands:

.. code-block:: python
   python setup.py build
   sudo python setup.py install

Alternatively, to install locally for a specific user (or if you do not have root privileges), enter these commands:

.. code-block:: python
   python setup.py build
   python setup.py install --user 

Optional tests may be run with the command (which may or may not require ``sudo``, depending on your install choice), 

.. code-block:: python
   python setup.py test


Contents
-----------

.. toctree::
   :maxdepth: 2

   modules

        


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

