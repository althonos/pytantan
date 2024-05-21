# 🐍🔁 PyTantan [![Stars](https://img.shields.io/github/stars/althonos/pytantan.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pytantan/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [Tantan](https://gitlab.com/mcfrith/tantan), a fast method for identifying repeats in DNA and protein sequences.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pytantan/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pytantan/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pytantan?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pytantan/)
[![License](https://img.shields.io/badge/license-GPLv3+-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pytantan.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pytantan)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pytantan?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pytantan)
[![AUR](https://img.shields.io/aur/version/python-pytantan?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pytantan)
[![Wheel](https://img.shields.io/pypi/wheel/pytantan.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pytantan/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pytantan.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pytantan/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pytantan.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pytantan/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pytantan/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003eaa?style=flat-square&maxAge=2678400)](https://git.lumc.nl/mflarralde/pytantan/)
[![Issues](https://img.shields.io/github/issues/althonos/pytantan.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pytantan/issues)
[![Docs](https://img.shields.io/readthedocs/pytantan/latest?style=flat-square&maxAge=600)](https://pytantan.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pytantan/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pytantan?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pytantan)


## 🗺️ Overview

[Tantan](https://gitlab.com/mcfrith/tantan) is a fast method developed
by Martin Frith[\[1\]](#ref1) to identify simple repeats in DNA or protein 
sequences. It can be used to mask repeat regions in reference sequences, and 
avoid false homology predictions between repeated regions.

PyTantan is a Python module that provides bindings to [Tantan](https://gitlab.com/mcfrith/tantan)
using [Cython](https://cython.org/). It implements a user-friendly, Pythonic
interface to mask a sequence with various parameters. It interacts with the 
Tantan interface rather than with the CLI, which has the following advantages:

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


## 🔧 Installing

PyTantan is available for all modern versions (3.6+), depending only on the
[`scoring-matrices`](https://pypi.org/project/scoring-matrices) package, and
optionally on the lightweight [`archspec`](https://pypi.org/project/archspec)
package for runtime CPU feature detection.

It can be installed directly from [PyPI](https://pypi.org/project/pytantan/),
which hosts some pre-built wheels for Linux and MacOS, as well as the code 
required to compile from source with Cython:
```console
$ pip install pytantan
```

<!-- Otherwise, PyTantan is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pytantan
``` -->

Check the [*install* page](https://pytantan.readthedocs.io/en/stable/install.html)
of the documentation for other ways to install PyTantan on your machine.

## 💡 Example

The top-level function `pytantan.mask_repeats` can be used to mask a sequence
without having to manage intermediate objects:

```python
import pytantan
masked = pytantan.mask_repeats("ATTATTATTATTATT")
print(masked)                 # ATTattattattatt
```

The mask symbol (and other parameters) can be given as keyword arguments:

```python
import pytantan
masked = pytantan.mask_repeats("ATTATTATTATTATT", mask='N')
print(masked)                 # ATTNNNNNNNNNNNN
```

To mask several sequences iteratively with the same parameters, consider 
creating a `RepeatFinder` once and calling the `mask_repeats` method for 
each sequence to avoid resource re-initialization.

<!-- See the [API documentation](https://pytantan.readthedocs.io/en/stable/api/index.html) 
for more examples, including how to use the internal API, and detailed 
reference of the parameters and result types. -->

<!-- ## 🧶 Thread-safety -->

<!-- ## ⏱️ Benchmarks -->


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/pytantan/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pytantan/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pytantan/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
Tantan is developed by [Martin Frith](https://sites.google.com/site/mcfrith/martin-frith) and is distributed under the
terms of the GPLv3 or later as well. See `vendor/tantan/COPYING.txt` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [Tantan authors](https://github.com/Martinsos). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [Leiden University Medical Center](https://www.lumc.nl/en/) in
the [Zeller team](https://github.com/zellerlab).*


## 📚 References

- <a id="ref1">\[1\]</a> Frith, Martin C. “A new repeat-masking method enables specific detection of homologous sequences.” Nucleic acids research vol. 39,4 (2011): e23. [doi:10.1093/nar/gkq1212](https://doi.org/10.1093/nar/gkq1212)

