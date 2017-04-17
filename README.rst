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
All the intensive computations were migrated to Cython and paralellized.
   
The SkewT Python package was a cornerstone of this project.  
We are grateful to all its collaborators.


*Technology builds on technology*


Dependencies
============

The SkewTplus package need the following dependencies

* matplotlib
* numpy
* cython
* libgcc

For running the examples:
* Basemap
* hdf4
* netCDF4


To-Do List
==========
* More column diagnostics.
* Hodographs? Anyone? 


Contributors
============
* Ross Bunn from Monash University is actively developing and finding all my 
  warty bugs.
* Gokhan Sever (North Carolina) is a keen user and has been encouraging me 
  to add more stuff. It's thanks to him that I have finally implemented the 
  CAPE routines.
* Simon Caine.
* Hamish Ramsay (Monash) has promised to at least think about adding some 
  extra diagnostics.
* Holger Wolff as tester


