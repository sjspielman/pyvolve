.. pyvolve documentation master file, created by
   sphinx-quickstart on Mon Jan 19 10:26:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pyvolve's documentation!
===================================

Pyvolve is a python library for simulating sequences along a phylogeny using continuous-time Markov models. Pyvolve implements a variety of `standard modeling frameworks <http://www.molecularevolution.org/resources/models/>`_:


1. Nucleotide models (GTR and all nested variants)
2. Empirical amino acid models (JTT, WAG, and LG)
3. Codon (dN/dS) models (GY94-style and MG94-style)
4. The ECM empirical codon model
5. Mutation-selection models (Halpern and Bruno, 1998), implemented for both nucleotides and codons

In particular, the Pyvolve framework allows you to heavily customize your model, both in terms of stationary frequency parameters and rate matrix elements. Further, you can specify your own rate matrix rather than relying on built-in models.
A detailed user-manual for Pyvolve is available `here <https://github.com/sjspielman/pyvolve/user_manual/pyvolve_manual.pdf>`_.

Download and Installation
--------------------------

Pyvolve is freely available under a FreeBSD license from `https://github.com/sjspielman/pyvolve/releases <https://github.com/sjspielman/pyvolve/releases>`_. Once pyvolve has been downloaded, navigate to the directory in the terminal. To install for all users, enter these commands:


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


Dependencies
-------------

Pyvolve has several dependencies which must be installed and in your ``PYTHONPATH`` for pyvolve to run properly. 

1. `NumPy <https://www.numpy.org/>`_
2. `SciPy <https://www.scipy.org/>`_
3. `BioPython <http://biopython.org/wiki/Main_Page>`_

Note that both NumPy and SciPy can be installed with ``pip``, but you must download BioPython from source.

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

