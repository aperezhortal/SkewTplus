=========================================================
SkewTplus -- Atmospheric Profile Plotting and Diagnostics
=========================================================

SkewTplus provides useful tools for plotting and analyzing atmospheric sounding data.
In particular, it provides useful classes to handle the awkward skew-x projection.
        
This package is based on the SkewT Python package developed by Thomas Chubb
(https://pypi.python.org/pypi/SkewT/)

This package is based on the SkewT Python package developed by Thomas Chubb.
The vertical soundings plots are handled by a new class (SkewT) that provides an
interface similar to the `matplolib's pyplot <http://matplotlib.org/api/pyplot_api.html>`_ api.
This new class allows more control over the Figure properties like multiple plots (normal axis and Skew-T axis).

Also, the **thermodynamics** module was improved.
All the intensive computations were migrated to Cython and parallelized.

*Technology builds on technology*

Documentation
=============

Check the documentation at http://skewtplus.readthedocs.io/en/latest/

Dependencies
============

The SkewTplus package need the following dependencies

* matplotlib
* numpy
* cython (optional)
* netCDF4
* six
* future (python2)
* hdf4
* libgcc >=5
* requests

For running the WRF data example:

* Basemap



Installing SkewTplus
====================

IMPORTANT - OSX installation
----------------------------

Before installing the package, be sure that Numpy is installed.
Then, install the apple's Xcode application by running::
    xcode-select --install

Before running the pip or the setup commands execute::

    export CC=clang ; export CXX=clang

Then you can continue with any of the following installation procedures.
 
Nevertheless **pip** is highly recommended.



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

Conda install - Only available linux users
------------------------------------------

If you are using an anaconda environment, to install the package execute::
    
    conda install -c andresperezcba skewtplus
    

Contributions
===========

SkewTplus is an open source software project.
Contributions to the package are welcomed from all users.
Feel free to suggest enhancements or report bugs by opening an issue in the github project page: 

https://github.com/aperezhortal/SkewTplus/issues


Code
----

The latest source code can be obtained with the command::

    git clone https://github.com/aperezhortal/SkewTplus.git

If you are planning on making changes that you would like included in SkewTplus,
forking the repository is highly recommended.






