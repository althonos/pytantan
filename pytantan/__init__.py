# noqa: D104
from ._version import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3+"
__all__ = [
    "Alphabet",
    "RepeatFinder",
    "LikelihoodMatrix",
    "mask_repeats",
    "default_scoring_matrix",
]

from . import lib
from ._mask import mask_repeats
from .lib import Alphabet, RepeatFinder, LikelihoodMatrix, default_scoring_matrix

# Use the library documentation
__doc__ = lib.__doc__

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pytantan.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
