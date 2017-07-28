#!/usr/bin/python
"""

@author: mcosta

"""
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="spiops",

    # Version number (initial):
    version="0.1.0",

    description="Extension of SPICE functionalities for ESA Planetary Missons",
    long_description=open('README.md').read(),
    license= open(path.join(here,'LICENSE')).read(),
    data_files = [('', ['LICENSE'])],

    # Application author details:
    author="Marc Costa Sitja, ESAC/ESA",
    author_email="marc.costa@esa.int",
    url="https://github.com/marcsit/spiops",


    # Classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science Operations Engineers and Developers, Planetary Scientists',
        'Topic :: Planetary Science :: Geometry Computations',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],

    # Keywords
    keywords=['esa', 'spice', 'naif', 'planetary', 'space', 'geometry'],

    # Packages
    packages=["spiops"],

    # Include additional files into the package
    include_package_data=False,

    # Dependent packages (distributions)
    install_requires=['spiceypy'],
    python_requires='>=3',
)