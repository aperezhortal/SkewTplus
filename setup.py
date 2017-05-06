


from __future__ import division, absolute_import, print_function


from setuptools import setup, find_packages
from setuptools.extension import Extension

try:
    import numpy
except ImportError:
    raise RuntimeError( "Numpy required to pior running the package installation\n" +
                        "Try installing it with:\n" + 
                        "$> pip install numpy" )
    
    
try:
    from Cython.Build.Dependencies import cythonize
    CythonPresent = True
except ImportError:
    CythonPresent = False
    

if CythonPresent:
    
    thermodynamicsLibExtension = Extension( "SkewTplus._thermodynamics",
                                            sources = ['SkewTplus/_thermodynamics.pyx'],
                                            extra_compile_args = ['-fopenmp'],
                                            extra_link_args = ['-fopenmp'] ,
                                            include_dirs=[numpy.get_include()],
                                            language='c++') 
                                           
    externalModules = cythonize([thermodynamicsLibExtension])                                       
else:
    thermodynamicsLibExtension = Extension( "SkewTplus._thermodynamics",
                                            sources = ['SkewTplus/_thermodynamics.cpp'],
                                            extra_compile_args = ['-fopenmp'],
                                            extra_link_args = ['-fopenmp'] ,
                                            include_dirs=[numpy.get_include()],
                                            language='c++') 
    
    externalModules = [thermodynamicsLibExtension]

build_requires=['matplotlib','numpy']

setup(
    name='SkewTplus',
    version='1.1.7',
    author = "Andres Perez Hortal",
    author_email = "andresperezcba@gmail.com",
    packages=find_packages(),
    ext_modules = externalModules,
    url='http://pypi.python.org/pypi/SkewTplus/',
    license='LICENSE.txt',
    description='Atmospheric Profile Plotting and Diagnostics',
    long_description=open('README.rst').read(),
    classifiers=[
    'Development Status :: 5 - Production/Stable',    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Cython'],
    setup_requires=build_requires,
    install_requires=build_requires
    )


