"""
Setup 
"""
from __future__ import division, absolute_import, print_function

__author__ = "Andres Perez Hortal"
__copyright__ = "Copyright (c) 2017, Andres A. Perez Hortal"
__license__ = "BSD-3-Clause License, see LICENCE.txt for more details"

from setuptools import setup, find_packages
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
except ImportError:
    raise RuntimeError(
        "Cython required for running the package installation\n"
        + "Try installing it with:\n"
        + "$> pip install cython"
    )

try:
    import numpy
except ImportError:
    raise RuntimeError(
        "Numpy required for running the package installation\n"
        + "Try installing it with:\n"
        + "$> pip install numpy"
    )

_thermodynamics_ExtensionArguments = dict(
    extra_compile_args=["-fopenmp"],
    extra_link_args=["-fopenmp"],
    include_dirs=[numpy.get_include()],
    language="c++",
)

thermodynamicsLibExtension = Extension(
    "SkewTplus._thermodynamics",
    sources=["SkewTplus/_thermodynamics.pyx"],
    **_thermodynamics_ExtensionArguments
)

externalModules = cythonize([thermodynamicsLibExtension])

build_requires = ["matplotlib>=3.1", "cython", "numpy", "netCDF4", "requests", "six"]

setup(
    name="SkewTplus",
    version="1.1.3",
    author="Andres Perez Hortal",
    packages=find_packages(),
    ext_modules=externalModules,
    url="http://pypi.python.org/pypi/SkewTplus/",
    license="LICENSE.txt",
    description="Atmospheric Profile Plotting and Diagnostics",
    long_description=open("README.rst").read(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Cython",
    ],
    setup_requires=build_requires,
    install_requires=build_requires,
)
