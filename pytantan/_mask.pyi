import typing
from typing import Union, Optional

from scoring_matrices import ScoringMatrix

def mask_repeats(
    sequence: Union[str, bytes, bytearray, memoryview],
    *,
    protein: bool = False,
    scoring_matrix: Optional[ScoringMatrix] = None,
    match_score: Optional[float] = None,
    mismatch_cost: Optional[float] = None,
    repeat_start: float = 0.005,
    repeat_end: float = 0.05,
    decay: float = 0.9,
    mask: Optional[str] = None,
    threshold: float = 0.5,
) -> str: ...
