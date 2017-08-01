Welcome to spiops's documentation!
==================================

spiops is a library aimed to help scientists and engineers that deal with
Solar System Geometry mainly for planetary science. More in particular is
aimed to assists the users to extend the usage of SPICE.

Functionalities vary from the computation of the illumination of a given
Field-of-View to obtaining the coverage of a given S/C for a particular
meta-kernel.

The underlying idea of spiops is to be used as a multi-user and
multi-disciplinary pool of re-usable SPICE based functions to provide cross
mission and discipline support of SPICE for ESA Planetary and Heliophyiscs
missions.

Feedback and new functionalities are always welcome, if you discover that a
function is not
working as expected or if you have a function that you believe can be of
interest to other people please open an issue or contact marc.costa@esa.int.

spiops is developed and maintained by the ESA SPICE Service (ESS)
http://spice.esac.esa.int.

SPICE is an essential tool for scientists and engineers alike in the
planetary science field for Solar System Geometry. Please visit the NAIF
website  for more details about SPICE.

Contents
********

.. toctree::
   :maxdepth: 2

   documentation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Installation
************

First install the dependencies (spiceypy, pytest) for the project. Then run

.. code-block:: python

   pip3 install spiops

to install from pypi. If you wish to install spiops from source first
download or clone the project from: https://mcosta@repos.cosmos.esa.int/socci/scm/spice/spiops.git

Then run

.. code-block:: python

   python3 setup.py install.



Citations
*********

If the use of spiops is used in a publication, please consider
citing spiops, SpiceyPy and the SPICE toolkit. The citation information
for SPICE can be found on the NAIF website and the citation information for
SpiceyPy in the GitHub repository.


To cite SpiceyPy use the following DOI: |Citation|

.. |Citation| image:: https://zenodo.org/badge/16987/AndrewAnnex/SpiceyPy.svg
   :target: https://zenodo.org/badge/latestdoi/16987/AndrewAnnex/SpiceyPy