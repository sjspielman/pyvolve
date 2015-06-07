.. pyvolve documentation master file, created by
   sphinx-quickstart on Mon Jan 19 10:26:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pyvolve's documentation!
===================================

Pyvolve is a python library for simulating sequences along a phylogeny using continuous-time Markov models. Pyvolve is hosted on github: `https://github.com/sjspielman/pyvolve <https://github.com/sjspielman/pyvolve>`_. Pyvolve implements a variety of `standard modeling frameworks <http://www.molecularevolution.org/resources/models/>`_:


1. Nucleotide models (GTR and all nested variants)
2. Empirical amino acid models (JTT, WAG, and LG)
3. Codon (dN/dS) models (GY94-style and MG94-style)
4. The ECM empirical codon model
5. Mutation-selection models (Halpern and Bruno, 1998), implemented for both nucleotides and codons

In particular, the Pyvolve framework allows you to heavily customize your model, both in terms of stationary frequency parameters and rate matrix elements. Further, you can specify your own rate matrix (and indeed the states to evolve) rather than relying on built-in models.
A detailed user-manual for Pyvolve is available `here <https://github.com/sjspielman/pyvolve/tree/master/user_manual/pyvolve_manual.pdf>`_.

Download and Installation
--------------------------

Pyvolve is freely available under a FreeBSD license. The easiest way to install Pyvolve is using `pip` or `easy_install`:

.. code-block:: bash

   sudo pip install pyvolve
   #or
   sudo easy_install install pyvolve
   
Alternatively, the most recent release of Pyvolve is available for download from `https://github.com/sjspielman/pyvolve/releases <https://github.com/sjspielman/pyvolve/releases>`_. Once Pyvolve has been downloaded, navigate to the directory in the terminal. To install for all users, enter this command:


.. code-block:: bash

   sudo python setup.py install


Alternatively, to install locally for a specific user (or if you do not have root privileges), enter this command:

.. code-block:: python
   
   python setup.py install --user 


Optional tests may be run with the command (which may or may not require ``sudo``, depending on your install choice), 

.. code-block:: python
   
   python setup.py test


Dependencies
-------------

Pyvolve has several dependencies which must be installed: 

1. `NumPy <http://www.numpy.org/>`_
2. `SciPy <http://www.scipy.org/>`_
3. `BioPython <http://biopython.org/wiki/Main_Page>`_

Note that installing Pyvolve with pip and/or easy_install will give you these dependencies if they are missing.

Issues and Questions
---------------------

Please file all bugs, issues, feature requests, and/or questions in the Issues section of Pyvolve's github repository: `https://github.com/sjspielman/pyvolve/issues <https://github.com/sjspielman/pyvolve/issues>`_. Note that you will need a github account to file.


Citation
---------
If you use Pyvolve, please cite it as:

Spielman, SJ and Wilke, CO. 2015. Pyvolve: A flexible Python module for simulating sequence along phylogenies. biorXiv. doi: http://dx.doi.org/10.1101/020214.

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

