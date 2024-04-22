# noqa: D104

from . import (
    test_alphabet,
    test_doctest,
    test_score_matrix,
)


def load_tests(loader, suite, pattern):
    test_doctest.load_tests(loader, suite, pattern)
    suite.addTests(loader.loadTestsFromModule(test_alphabet))
    suite.addTests(loader.loadTestsFromModule(test_score_matrix))
    return suite