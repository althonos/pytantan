PyTantan |Stars|
================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pytantan.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pytantan/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `Tantan <https://gitlab.com/mcfrith/tantan>`_, *a fast method for identifying repeats in DNA and protein sequences.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pytantan/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pytantan/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pytantan?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pytantan/

.. |PyPI| image:: https://img.shields.io/pypi/v/pytantan.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pytantan

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pytantan?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pytantan

.. |AUR| image:: https://img.shields.io/aur/version/python-pytantan?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pytantan

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pytantan?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pytantan/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pytantan.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pytantan/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pytantan.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pytantan/#files

.. |License| image:: https://img.shields.io/badge/license-GPLv3+-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/gpl-3.0/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pytantan/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pytantan.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pytantan/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pytantan?style=flat-square&maxAge=3600
   :target: http://pytantan.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pytantan/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/pytantan?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pytantan


.. currentmodule:: scoring_matrices


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

- **no binary dependency**: PyTantan is distributed as a Python package, so
  you can add it as a dependency to your project, and stop worrying about the
  `tantan` binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the Tantan CLI using a sub-process
  and temporary files.
- **better portability**: Tantan uses SIMD to accelerate alignment scoring, 
  but doesn't support dynamic dispatch, so it has to be compiled on the local
  machine to be able to use the full capabilities of the local CPU. PyTantan
  ships several versions of Tantan instead, each compiled with different 
  target features, and selects the best one for the local platform at runtime.


Setup
-----

Run ``pip install pytantan`` in a shell to download the latest release 
from PyPi, or have a look at the :doc:`Installation page <install>` to find 
other ways to install ``pytantan``.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Contributing <contributing>
   API Reference <api/index>
   Changelog <changes>


License
-------

This library is provided under the `GNU General Public License v3.0 or later <https://choosealicense.com/licenses/gpl-3.0/>`_.
Tantan is developed by `Martin Frith <https://sites.google.com/site/mcfrith/martin-frith>`_ and is 
distributed under the terms of the GPLv3 or later as well. See 
`vendor/tantan/COPYING.txt` in the source distribution for more information.

*This project was developed by* `Martin Larralde <https://github.com/althonos/>`_
*during his PhD project at the* `Leiden University Medical Center <https://www.lumc.nl/en/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
