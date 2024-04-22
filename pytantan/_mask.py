import collections

from .lib import Alphabet, ScoreMatrix, RepeatFinder

_DNA = Alphabet.dna()
_PROTEIN = Alphabet.protein()


def mask_repeats(
    sequence,
    *,
    protein=False,
    # preserve_case=False,
    score_matrix=None,
    match_score=None,
    mismatch_cost=None,
    repeat_start=0.005,
    repeat_end=0.05,
    repeat_period=None,
    decay=0.9,
    # gap_open=None,
    # gap_extend=None,
    mask=None,
    threshold=0.5,
):
    """Mask regions predicted as repeats in the given sequence.

    Arguments:
        sequence (`str` or byte-like object): The sequence containing
            the repeats to mask.
        protein (`bool`): Set to `True` to treat the input sequence as
            a protein sequence.
        score_matrix (`str` or `~pytantan.ScoreMatrix`): A score matrix
            to use for scoring character matches and mismatches. Either
            pass a matrix name (such as ``BLOSUM62``) to load built-in
            matrix, or a pre-initialized `ScoreMatrix` object.
        match_score (`int`): The score for scoring character matches.
            Must be set along `mismatch_cost`. Incompatible with the
            `score_matrix` option.
        match_score (`int`): The penalty for scoring character mismatches.
            Must be set along `match_score`. Incompatible with the
            `score_matrix` option.
        repeat_start (`float`): The probability of a repeat starting
            per position.
        repeat_end (`float`): The probability of a repeat ending per
            position.
        decay (`float`): The probability decay per period.
        threshold (`float`): The probability threshold above which to
            mask sequence characters.
        mask (`str` or `None`): A single mask character to use for 
            masking positions. If `None` given, masking uses the
            lowercase letters of the original sequence.

    """
    if match_score is not None and mismatch_cost is None:
        raise ValueError("Cannot set `match_score` without setting `mismatch_cost`")
    elif match_score is None and mismatch_cost is not None:
        raise ValueError("Cannot set `mismatch_cost` without setting `match_score`")
    if score_matrix is not None and match_score is not None:
        raise ValueError("Cannot set both `score_matrix` and `match_score`")

    if isinstance(score_matrix, ScoreMatrix):
        matrix = score_matrix
    elif isinstance(score_matrix, str) and protein:

        matrix = ScoreMatrix(score_matrix)
    elif match_score is not None:
        matrix = ScoreMatrix.match_mismatch(
            _PROTEIN if protein else _DNA, 
            match_score=match_score, 
            mismatch_cost=mismatch_cost
        )
    else:
        matrix = ScoreMatrix.match_mismatch(
            _PROTEIN if protein else _DNA
        )

    repeat_finder = RepeatFinder(
        matrix,
        repeat_start=repeat_start,
        repeat_end=repeat_end,
        decay=decay
    )        
    return repeat_finder.mask_repeats(
        sequence,
        threshold=threshold,
        mask=mask,
    )