This directory contains some possibly useful post-processing modules for computing dN/dS in various ways:

- *compute_dnds_from_mutsel.py* contains a module to compute a dN/dS value from mutation-selection model parameters, as described in Spielman and Wilke 2015 (doi: 10.1093/molbev/msv003). See the module for usage details.

- *count_simulated_dnds.py* contains a module to ascertain the actual simulated site-specific dN/dS values from a pyvolve simulation. This script uses counting approach, similar to the SLAC method [from Kosakovsky Pond and Frost, 2005 (doi: 10.1093/molbev/msi105)], although this module uses the known ancestral sequences as simulated in pyvolve. See the module for usage details.