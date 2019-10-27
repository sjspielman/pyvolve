Pyvolve
============

Pyvolve is an open-source Python module for simulating sequences along a phylogenetic tree according to continuous-time Markov models of sequence evolution. **Please ensure you are using the most up-to-date version of Pyvolve! The current version is 1.0.0.** 

A detailed user manual for Pyvolve is available [here](https://github.com/sjspielman/pyvolve/raw/master/user_manual/pyvolve_manual.pdf), and API documentation for Pyvolve is available at [https://sjspielman.github.io/pyvolve](https://sjspielman.github.io/pyvolve).

Pyvolve has several dependencies:
* [Biopython](http://biopython.org/wiki/Download)
* [Scipy and Numpy](http://www.scipy.org/install.html) (with Numpy >= 1.7)

Currently, `pyvolve` development is occurring in Python3.x with absolutely no plans to return to Python2. The following examples therefore assume Python3 (which you may have on your computer named as `python3` and not `python`!)

You can install Pyvolve directly using `pip` or `easy_install` (note that, if needed, these lines will install any missing dependencies for you!):
```bash
pip install pyvolve
# OR
easy_install install pyvolve
```
Note that these commands might need `sudo` in front, depending on your permissions.

To update your version of Pyvolve (for pip), simply use the `--upgrade` argument (again, possibly w/ sudo):
```bash
pip install --upgrade pyvolve
```

Alternatively, you can download and install from source. Download the most recent version of Pyvolve from the Releases tab ([https://github.com/sjspielman/pyvolve/releases](https://github.com/sjspielman/pyvolve/releases)), uncompress the file, and navigate into the `pyvolve/` directory. From this directory, enter the following commands:
```bash
sudo python setup.py install
sudo python setup.py test  # optional, but recommended. Tests implemented **for Python2 only**
```

If you do not have root privileges, you can install Pyvolve for only you (the user!) with this line instead:
```bash
python setup.py install --user
```

Please file any bugs and/or relay any questions under the Issues tab: [https://github.com/sjspielman/pyvolve/issues](https://github.com/sjspielman/pyvolve/issues).

*If you use Pyvolve in your research, please cite the following:* <br>
Spielman, SJ and Wilke, CO. 2015. [**Pyvolve: A flexible Python module for simulating sequences along phylogenies**](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139047). PLOS ONE. 10(9): e0139047.

```
@article{SpielmanWilke2015,
author = {Spielman, S. J. and Wilke, C. O.},
title = {Pyvolve: A Flexible Python Module for Simulating Sequences along Phylogenies},
journal = {PLOS ONE},
year = {2015},
volume = 10,
pages = {e0139047}
}
```


