# noqa: D104
from ._version import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3+"
__all__ = [
    "Alphabet",
    "RepeatFinder",
    "mask_repeats",
]

from . import lib
from ._mask import mask_repeats
from .lib import Alphabet, RepeatFinder
