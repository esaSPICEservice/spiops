spiops
======

spiops is a library aimed to help scientists and engineers that deal with 
Solar System Geometry mainly for planetary science. More in particular is 
aimed to assists the users to extend the usage of SPICE.

SPICE is an essential tool for scientists and engineers alike in the 
planetary science field for Solar System Geometry. Please visit the NAIF 
website  for more details about SPICE.
 

Introduction
------------

spiops is a python library thar uses [SpiceyPy](<https://github
.com/AndrewAnnex/SpiceyPy>) to use [SPICE Toolkit](<https://naif.jpl.nasa 
.gov/naif/>) APIS to provide higher-level functions than the ones available 
with SPICE. This functions have been idenfitied from having to implement 
multiple time a series of SPICE APIs to obtain a given derived functionality.

Functionalities vary from the computation of the illumination of a given 
Field-of-View to obtaining the coverage of a given S/C for a particular 
meta-kernel.

The underlying idea of spiops is to be used as a multi-user and 
multi-disciplinary pool of re-usable SPICE based functions to provide cross 
mission and discipline support of SPICE for ESA Planetary and Heliophyiscs 
missions. 

spiops is developed and maintained by the ESA SPICE Service (ESS) 
<https://spice.esac.esa.int>.


Installation
------------

First install the dependencies (spiceypy, pytest) for the project. Then
run ``pip install spiops`` to install from pypi.

If you wish to install spiceypy from source first download or clone the project. Then run ``python setup.py install``.
To uninstall run ``pip uninstall spiceypy``.

Documentation
-------------

The spiops docs are available at:
[spiops.readthedocs.org](<http://spiops.readthedocs.org>).
The documentation for spiops is intentionally abridged so as to utilize the 
excellent [documentation provided by the
NAIF](<http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/index.html>).
Please refer to C and IDL documentation available on the NAIF website
for in-depth explanations. Most of the functions docstrings have a 
link to the
corresponding C function in the NAIF docs that the function is based on.

How to Help
-----------

Feedback and new functionalities are always welcome, if you discover that a 
function is not 
working as expected or if you have a function that you believe can be of 
interest to other people please open an issue or contact [me](<marc.costa@esa
.int>).

Citing spiops
--------------

If the use of spiops is used in a publication, please consider
citing spiops, SpiceyPy and the SPICE toolkit. The citation information
for SPICE can be found on the NAIF website and the citation information for 
SpiceyPy in the GitHub repository.



Known Working Environments:
---------------------------

spiops is compatible with modern Linux and Mac. At this development stage and
with the available resources is hard to guarantee an extension of the 
compatible working environments. If you run into issues with your system 
please submit an issue with details. 

- OS: OS X, Linux
- CPU: 64bit
- Python 3.5

Acknowledgements
----------------

spiops makes an intensive usage of SpiceyPy, a Python wrapper for the 
NAIF C SPICE Toolkit (N66) by [AndrewAnnex](<https://github.com/AndrewAnnex/).

We have to be very grateful to Andrew for having brought SPICE to Python!
