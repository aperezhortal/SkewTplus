=========================================================
SkewTplus -- Atmospheric Profile Plotting and Diagnostics
=========================================================

SkewTplus provides a few useful tools to help with the plotting and analysis of 
upper atmosphere data. In particular it provides some useful classes to 
handle the awkward skew-x projection.
        
This package is based in the SkewT Python package developed by Thomas Chubb
(https://pypi.python.org/pypi/SkewT/)
        
The main difference with the original *SkewT package* is that the vertical soundings 
plots are handled by an special class (SkewT).
The new *SkewT* class extends the base
`matplolib's Figure <http://matplotlib.org/api/figure_api.html?highlight=figure#module-matplotlib.figure>`_
class with an interface similar to 
`matplolib's pyplot <http://matplotlib.org/api/pyplot_api.html>`_.
It also allows to create Skew-T type plots in a simple way.
This new class allows a complete control over the Figure properties like
multiple plots (normal axis and Skew-T axis).

In addition, the **thermodynamics** module was improved.
All the intensive computations were migrated to Cython and parallelized.
   
The SkewT Python package was a cornerstone of this project.  
We are grateful to all its collaborators.


*Technology builds on technology*


Documentation
=============

The documentations is separated in two big branches. 
The :ref:`user-reference` and the :ref:`developer-reference`.
The user reference provides a quick overview of the most important features of 
the package. For more detailed and a comprehensive understanding of the 
package the reader must consult the :ref:`developer-reference`.


.. toctree::
    :maxdepth: 1
    
    userReference/index
    developerReference/index
    license
    

Dependencies
============

The SkewTplus package need the following dependencies

* matplotlib
* numpy
* cython
* wxPython Phoenix (https://github.com/wxWidgets/Phoenix)
* libgcc

The wxPython Phoenix project is needed for python 3 compatibility.

Installing SkewT
================

If wxPython Phoenix  it's not installed, to install it manually using pip
( could take several minutes)::
    
    pip install --upgrade --pre -f https://wxpython.org/Phoenix/snapshot-builds/ wxPython_Phoenix

If all the dependencies are installed, to install the SkewTplus package
you can download the tarball from the github repository and run::

    python setup.py install

If you want to put it somewhere different than your system files, you can do::
    
    python setup.py install --prefix=/path/to/local/dir


Conda  install: hopefully coming soon
    
    



    




