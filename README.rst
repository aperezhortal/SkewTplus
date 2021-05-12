=========================================================
SkewTplus -- Atmospheric Profile Plotting and Diagnostics
=========================================================

SkewTplus provides useful tools for plotting and analyzing atmospheric sounding data.
In particular, it provides useful classes to handle the awkward skew-x projection.
        
This package is based on the SkewT Python package developed by Thomas Chubb
(https://pypi.python.org/pypi/SkewT/).
The vertical soundings plots are handled by a new class (SkewT) that provides an
interface similar to the `matplolib's pyplot <http://matplotlib.org/api/pyplot_api.html>`_ api.
This new class allows more control over the Figure properties like multiple plots (normal axis and Skew-T axis).

Also, the **thermodynamics** module was improved.
All the intensive computations were migrated to Cython and parallelized.

*Technology builds on technology*

Documentation
=============

Check the documentation at http://skewtplus.readthedocs.io/en/latest/

Installing SkewTplus
====================

Dependencies
------------

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


OSX users: GNU gcc compiler installation
----------------------------------------

(Adapted from the `Pysteps installation instructions <https://pysteps.readthedocs.io/en/latest/user_guide/install_pysteps.html#osx-users>`_)

SkewTplus uses a Cython extension that need to be compiled with multi-threading
support enabled.
The default OSX Clang compiler does not support OpenMP and the package installation
will fail with an error similar to::

    clang: error: unsupported option '-fopenmp'
    error: command 'gcc' failed with exit status 1

To avoid this error, install the latest GNU gcc using
Homebrew_::

    brew install gcc

.. _Homebrew: https://brew.sh/

This version has multi-threading enabled. Once the gcc compiler is installed, we need
to make sure that this compiler is used during the package installation.
To that end, we have to define the following environment variables::

    # Assuming that gcc version 9 was installed by homebrew
    export CC=gcc-9
    export CXX=g++-9



Under certain circumstances, Homebrew_ does not add the symbolic links for the
gcc executables under `/usr/local/bin`.
If this is the case, you will need to specify the CC and CCX variables using the full
path to the installed gcc. For example::

    export CC=/usr/local/Cellar/gcc/<version>/bin/gcc-9
    export CXX=/usr/local/Cellar/gcc/<version>/bin/g++-9

To find the full gcc path, run `which gcc-9`

Then you can continue with any of the following installation procedures.
 

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


Contributions
=============

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






