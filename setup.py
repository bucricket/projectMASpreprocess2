#!/usr/bin/env python

from __future__ import print_function
import os


    
try:
    from setuptools import setup
    setup_kwargs = {'entry_points': {'console_scripts':['preprocess=preprocess.prepareData:main']}}
except ImportError:
    from distutils.core import setup
    setup_kwargs = {'scripts': ['bin/preprocess']}
    
from preparepydisalexi import __version__


#=============setup the python scripts============================



setup(
    name="projectmaspreprocess",
    version=__version__,
    description="prepare data for input to pyDisALEXI",
    author="Mitchell Schull",
    author_email="mitch.schull@noaa.gov",
    url="https://github.com/bucricket/projectMASpreprocess2.git",
    packages= ['preprocess'],
    platforms='Posix; MacOS X; Windows',
    license='BSD 3-Clause',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        # Uses dictionary comprehensions ==> 2.7 only
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: GIS',
    ],  
    **setup_kwargs
)

