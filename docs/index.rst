PyTantan |Stars|
================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pytantan.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pytantan/stargazers
   :class: dark-light

`Cython <https://cython.org/>`_ *bindings and Python interface to* `Tantan <https://gitlab.com/mcfrith/tantan>`_, *a fast method for identifying repeats in DNA and protein sequences.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pytantan/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pytantan/actions
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pytantan?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pytantan/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/pytantan.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pytantan
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pytantan?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pytantan
   :class: dark-light

.. |AUR| image:: https://img.shields.io/aur/version/python-pytantan?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pytantan
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pytantan?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pytantan/#files
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pytantan.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pytantan/#files
   :class: dark-light

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pytantan.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pytantan/#files
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pytantan/
   :class: dark-light

.. |Mirror| image:: https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=3600&style=flat-square
   :target:https://git.lumc.nl/mflarralde/pytantan/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pytantan.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pytantan/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/pytantan?style=flat-square&maxAge=3600
   :target: http://pytantan.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pytantan/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://img.shields.io/pypi/dm/pytantan?style=flat-square&color=303f9f&maxAge=3600&label=downloads
   :target: https://pepy.tech/project/pytantan
   :class: dark-light


Overview
--------

`Tantan <https://gitlab.com/mcfrith/tantan>`_ is a fast method developed
by Martin Frith to identify simple repeats in DNA or protein sequences. It can 
be used to mask repeat regions in reference sequences, and avoid false homology 
predictions between repeated regions.

PyTantan is a Python module that provides bindings to 
`Tantan <https://gitlab.com/mcfrith/tantan>`_ using `Cython <https://cython.org/>`_. 
It implements a user-friendly, Pythonic interface to mask a sequence with 
various parameters. It interacts with the Tantan interface rather than with 
the CLI, which has the following advantages:

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pytantan`` as a ``pip`` or ``conda`` dependency, no need
      for the ``tantan`` binary or any external dependency.

   .. grid-item-card:: :fas:`gears` Flexible

      Pass any `str` or `bytes`-like object containing the raw 
      sequence as an input, and get a `str` as the output.

   .. grid-item-card:: :fas:`screwdriver-wrench` Configurable

      Use any scoring matrix from the `scoring-matrices <https://scoring-matrices.readthedocs.io/>`_ 
      package or build your own.

   .. grid-item-card:: :fas:`server` Parallel

      Easily run computations in parallel querying thread-safe 
      `~pytantan.RepeatFinder` with several sequences in parallel.

   .. grid-item-card:: :fas:`code-compare` Consistent

      Get the same results as Tantan! You are using the same code under 
      the hood.

   .. grid-item-card:: :fas:`microchip` Portable

      Get SIMD-acceleration on any supported platform without having
      to build the package from scratch.


Setup
-----

PyTantan is available for all modern Python versions (3.7+).

Run ``pip install pytantan`` in a shell to download the latest release 
from PyPi, or have a look at the :doc:`Installation page <guide/install>` to find 
other ways to install ``pytantan``.

Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


License
-------

This library is provided under the `GNU General Public License v3.0 or later <https://choosealicense.com/licenses/gpl-3.0/>`_.
Tantan is developed by `Martin Frith <https://sites.google.com/site/mcfrith/martin-frith>`_ and is 
distributed under the terms of the GPLv3 or later as well. See 
the :doc:`Copyright page <guide/copyright>` for more information.

*This project was developed by* `Martin Larralde <https://github.com/althonos/>`_
*during his PhD project at the* `Leiden University Medical Center <https://www.lumc.nl/en/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
