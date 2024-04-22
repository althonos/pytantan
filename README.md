# 🐍🌈🪨 PyOpal [![Stars](https://img.shields.io/github/stars/althonos/pytantan.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pytantan/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [Tantan](https://gitlab.com/mcfrith/tantan), a fast method for identifying repeats in DNA and protein sequences.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pytantan/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pytantan/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pytantan?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pytantan/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![PyPI](https://img.shields.io/pypi/v/pytantan.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pytantan)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pytantan?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pytantan)
[![AUR](https://img.shields.io/aur/version/python-pytantan?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pytantan)
[![Wheel](https://img.shields.io/pypi/wheel/pytantan.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pytantan/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pytantan.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pytantan/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pytantan.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pytantan/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pytantan/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pytantan/)
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

PyTantan is available for all modern versions (3.6+), optionally depending on
the lightweight Python package [`archspec`](https://pypi.org/project/archspec)
for runtime CPU feature detection.

<!-- It can be installed directly from [PyPI](https://pypi.org/project/pytantan/),
which hosts some pre-built x86-64 wheels for Linux, MacOS, and Windows,
Aarch64 wheels for Linux and MacOS, as well as the code required to 
compile from source with Cython:
```console
$ pip install pytantan
``` -->

<!-- Otherwise, PyTantan is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pytantan
``` -->

<!-- Check the [*install* page](https://pytantan.readthedocs.io/en/stable/install.html)
of the documentation for other ways to install PyTantan on your machine. -->

## 💡 Example

All classes are imported in the main namespace `pytantan`:
```python
import pytantan
```

Tantan needs a scoring matrix to score matches and mismatches between
sequence characters. Use the `ScoreMatrix` object to load a default one:
```python
matrix = pytantan.ScoreMatrix.dna()
```

A `RepeatFinder` can be created from a score matrix, and from additional
arguments to control the `tantan` parameters:
```python
tantan = pytantan.RepeatFinder(matrix)
```

The repeat finder can then be used to mask repeats inside sequences.
```
masked = tantan.mask("TGCAAGCTATTAGGCTTAGGTCAGTGCTTAGGCTTAGGTCAGTGCAACATA")
print(masked)       # TGCAAGCTATTAGGCTTAGGTCAGTGCttaagcttaggtcagtgcAACATA
```

<!-- The top-level function `pytantan.mask` can be used to mask a sequence
without having to manage intermediate objects. -->

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
at the [Leiden University Medical Center](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*


## 📚 References

- <a id="ref1">\[1\]</a> Korpar Matija, Martin Šošić, Dino Blažeka, Mile Šikić. SW#db: ‘GPU-Accelerated Exact Sequence Similarity Database Search’. PLoS One. 2015 Dec 31;10(12):e0145857. [doi:10.1371/journal.pone.0145857](https://doi.org/10.1371/journal.pone.0145857). [PMID:26719890](https://pubmed.ncbi.nlm.nih.gov/26719890). [PMC4699916](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4699916/).
- <a id="ref2">\[2\]</a> Rognes Torbjørn. Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation. BMC Bioinformatics. 2011 Jun 1;12:221. [doi:10.1186/1471-2105-12-221](https://doi.org/10.1186/1471-2105-12-221). [PMID:21631914](https://pubmed.ncbi.nlm.nih.gov/21631914/).[PMC3120707](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc3120707/).
