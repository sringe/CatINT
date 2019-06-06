Installation
============

.. |rtm| unicode:: U+00AE .. registered trademark sign .. 

CatINT is a program package which combines mass transport simulations with micro-
kinetic modeling. Installing CatINT requires:

- python 2.5 or greater
- numpy
- scipy
- matplotlib
- tables
- logging
.. - (os, sys, re, shutil, math, glob, imp, itertools, collections, subprocess)

After installing all dependencies, you can install CatINT into a folder $CATINT
via:

.. code:: bash

    $ mkdir $CATINT
    $ cd $CATINT
    $ git clone git@github.com:sringe/CatINT.git

As usual update your $PYTHONPATH variable with the $CATINT folder location and
verify that everything worked by importing the central transport module of
CatINT:

.. code::

    $ python >>>from catint import transport

For installing CatMAP, please see https://catmap.readthedocs.io/en/latest/installation.html.

For installing COMSOL |rtm|, please see www.comsol.de. CatiNT writes a java input file
for COMSOL |rtm| which is then compiled and executed. The COMSOL |rtm| executable location
needs to be updated in $CATINT/catint/comsol_wrapper.py so that CatINT finds it.

Documentation (this wiki) is located in the catmap/docs folder, and
is available online at http://catint.readthedocs.org/.
