



from __future__ import division, absolute_import, print_function

from Cython.Build.Dependencies import cythonize
from numpy.distutils.core import Extension, setup


thermodynamicsLibExtension = Extension( "SkewTplus._thermodynamics",
                                        sources = ['SkewTplus/_thermodynamics.pyx'],
                                        extra_compile_args = ['-fopenmp'],
                                        extra_link_args = ['-fopenmp'] ,
                                        language='c++')

setup(
    name='SkewTplus',
    version='1.1.0',
    author = "Andres Perez Hortal",
    author_email = "andresperezcba@gmail.com",    
    packages=['SkewTplus'],
    ext_modules = cythonize([ thermodynamicsLibExtension]),
    url='http://pypi.python.org/pypi/SkewTplus/',
    license='LICENSE.txt',
    description='Atmospheric Profile Plotting and Diagnostics',
    long_description=open('README.rst').read(),
    classifiers=[
    'Development Status :: 5 - Production/Stable',    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 2'
    'Programming Language :: Python :: 2.7'
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Cython'],
    install_requires=['matplotlib','numpy','cython','libgcc'], 
)


