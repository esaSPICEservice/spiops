#!/usr/bin/python
"""

@author: mcosta

"""
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.txt'), encoding='utf-8') as f:
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
    author_email="marc.costa@.esa.int",
    url="https://spice.esac.esa.int",

    # Packages
    packages=["spiops"],

    # Include additional files into the package
    include_package_data=False,

    # Dependent packages (distributions)
    install_requires=['spiceypy>=0.2.0'],
)