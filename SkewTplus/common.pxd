#!python
# cython: boundscheck=False
# cython: cdivision=True
# cython: embedsignature=True
"""
Common definitions and functions  
"""

__author__ = "Andres Perez Hortal"
__copyright__ = "Copyright (c) 2017, Andres A. Perez Hortal, McGill University"
__license__ = "BSD-3-Clause License, see LICENCE.txt for more details"
__email__ = "andresperezcba@gmail.com"


# Common definitions

cimport numpy as numpy
ctypedef numpy.float32_t float32


cdef inline float32 float_abs(float32 a) nogil: return a if a > 0. else -a
""" Return the absolute value of a float """

cdef inline float32 _linearInterpolation(float32 x, float32 x0, float32 x1, float32 y0, float32 y1) nogil:
    """
    Linear interpolation at x.
    y(x) = y0 + (x-x0) * (y1-y0) / (x1-x0)
    """
     
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    


    