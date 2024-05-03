# coding: utf-8
"""Test doctest contained tests in every file of the module.
"""

import configparser
import doctest
import importlib
import os
import pkgutil
import re
import shutil
import sys
import types
import warnings
from unittest import mock

from scoring_matrices import ScoringMatrix

import pytantan


def _load_tests_from_module(tests, module, globs, setUp=None, tearDown=None):
    """Load tests from module, iterating through submodules."""
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith("_")):
        if isinstance(attr, types.ModuleType):
            suite = doctest.DocTestSuite(
                attr,
                globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
            tests.addTests(suite)
    return tests


def load_tests(loader, tests, ignore):
    """`load_test` function used by unittest to find the doctests."""
    _current_cwd = os.getcwd()
    # demonstrate how to use Biopython substitution matrices without
    # actually requiring Biopython
    Bio = mock.Mock()
    Bio.Align = mock.Mock()
    Bio.Align.substitution_matrices = mock.Mock()
    Bio.Align.substitution_matrices.load = mock.Mock()
    Bio.Align.substitution_matrices.load.return_value = feng = mock.Mock()

    data = [ [-1 for _ in range(20)] for _ in range(20) ]
    for i in range(20):
        data[i][i] = 1

    feng.alphabet = "ARNDCQEGHILKMFPSTWYV"
    feng.__len__ = mock.Mock(return_value=20)
    feng.__iter__ = mock.Mock(wraps=data.__iter__)

    def setUp(self):
        warnings.simplefilter("ignore")
        # os.chdir(os.path.realpath(os.path.join(__file__, os.path.pardir, "data")))

    def tearDown(self):
        # os.chdir(_current_cwd)
        warnings.simplefilter(warnings.defaultaction)

    # doctests are not compatible with `green`, so we may want to bail out
    # early if `green` is running the tests
    if sys.argv[0].endswith("green"):
        return tests

    # recursively traverse all library submodules and load tests from them
    packages = [None, pytantan.lib]

    for pkg in iter(packages.pop, None):
        globs = dict(pyopal=pytantan, ScoringMatrix=ScoringMatrix, Bio=Bio, **pkg.__dict__)
        tests.addTests(
            doctest.DocTestSuite(
                pkg,
                globs=globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
        )
        # for (_, subpkgname, subispkg) in pkgutil.walk_packages(pkg.__path__):
        #     # do not import __main__ module to avoid side effects!
        #     if subpkgname == "__main__" or subpkgname.startswith("tests"):
        #         continue
        #     # import the submodule and add it to the tests
        #     module = importlib.import_module(".".join([pkg.__name__, subpkgname]))
        #     if hasattr(module, "__test__"):
        #         del module.__test__
        #     globs = dict(pyopal=pytantan, Bio=Bio, **module.__dict__)
        #     tests.addTests(
        #         doctest.DocTestSuite(
        #             module,
        #             globs=globs,
        #             setUp=setUp,
        #             tearDown=tearDown,
        #             optionflags=+doctest.ELLIPSIS,
        #         )
        #     )
        #     # if the submodule is a package, we need to process its submodules
        #     # as well, so we add it to the package queue, unless it is the
        #     # `pyopal.tests` package (we don't need to check for doctests there)
        #     if subispkg and subpkgname != "tests":
        #         packages.append(module)

    return tests
