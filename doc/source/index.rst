.. meta::

   :description  lang=en: The SkewTplus package provides tools to easily read 
    atmospheric sounding data from different formats (University of Wyoming and ARM)
    and create SkewT sounding plots along with parcel diagnostics (CAPE,CIN,etc.).
   
   :description  lang=es: El paquete SkewTplus provee herramientas para la facil manipulación
    y lectura de datos de sondeos verticales en diferentes formatos (University of Wyoming and ARM)
    y crear gráficos tipo SkewT de los sondeos verticales junto con los diagnosticos del sondeo (CAPE,CIN,etc.)
    
   :keywords: python, Skewt plot,skew-t plot, ARM sounding, University of Wyoming sounding


=========================================================
SkewTplus -- Atmospheric Profile Plotting and Diagnostics
=========================================================


The SkewTplus package provides tools to easily read atmospheric sounding data from 
different formats (University of Wyoming and ARM)
and create SkewT sounding plots along with parcel diagnostics (CAPE,CIN,etc.).
        
This package is based on the SkewT Python package developed by Thomas Chubb
(https://pypi.python.org/pypi/SkewT/)
        
The main difference with the original *SkewT package* is that the vertical soundings 
plots are handled by a special class (SkewT).
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

The documentation is separated in two big branches. 
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
* cython (optional)
* netCDF4

For running the examples:

* Basemap



Installing SkewTplus
====================

PIP install
-----------

To install the package using **pip** the numpy package must be already installed.
If is not installed, you can install it by running::
    pip install numpy

After the numpy package was installed, to install the SkewTplus package run::

    pip install SkewTplus


Install from source
-------------------

The latest version can be installed manually by downloading the sources from
https://github.com/aperezhortal/SkewTplus

To install the package manually, the numpy package must be already installed.
If is not installed, you can install it by running::
    pip install numpy
    
Then, you can install the SkewTplus package executing::

    python setup.py install

If you want to put it somewhere different than your system files, you can do::
    
    python setup.py install --prefix=/path/to/local/dir

IMPORTANT: If you install it using this way, all the dependencies need to be already installed! 

Conda install
-------------

If you are using an anaconda environment, to install the package execute::
    
    conda install -c andresperezcba skewtplus
    
    

Contributions
===========

SkewTplus is an open source software project.
Contributions to the package are welcomed from all users.
Feel free to suggest enhancements or report bugs by opening an issue in the github project page: 

https://github.com/aperezhortal/SkewTplus/issues

Thanks for using the SkewTplus package, for any feedback feel free to write to 
andresperezcba AT gmail DOT com


Code
----

The latest source code can be obtained with the command::

    git clone https://github.com/aperezhortal/SkewTplus.git

If you are planning on making changes that you would like included in SkewTplus,
forking the repository is highly recommended.




    




