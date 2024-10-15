import collections

from .lib import Alphabet, RepeatFinder, default_scoring_matrix

from scoring_matrices import ScoringMatrix


_DNA = Alphabet.dna()
_PROTEIN = Alphabet.protein()


def mask_repeats(
    sequence,
    *,
    protein=False,
    # preserve_case=False,
    scoring_matrix=None,
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
        scoring_matrix (`str` or `~scoring_matrices.ScoringMatrix`): A
            scoring matrix to use for scoring character matches and
            mismatches. Either pass a matrix name (such as ``BLOSUM62``)
            to load a built-in matrix, or a pre-initialized `ScoringMatrix`
            object.
        match_score (`int`): The score for scoring character matches.
            Must be set along `mismatch_cost`. Incompatible with the
            `scoring_matrix` option.
        match_score (`int`): The penalty for scoring character mismatches.
            Must be set along `match_score`. Incompatible with the
            `scoring_matrix` option.
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
    if scoring_matrix is not None and match_score is not None:
        raise ValueError("Cannot set both `scoring_matrix` and `match_score`")

    if isinstance(scoring_matrix, ScoringMatrix):
        matrix = scoring_matrix
    elif isinstance(scoring_matrix, str):
        matrix = ScoringMatrix.from_name(scoring_matrix)
    elif scoring_matrix is None:
        matrix = default_scoring_matrix(protein, match_score, mismatch_cost)
    else:
        ty = type(scoring_matrix).__name__
        raise TypeError(f"expected ScoringMatrix, str or None, got {ty}")

    repeat_finder = RepeatFinder(
        matrix,
        repeat_start=repeat_start,
        repeat_end=repeat_end,
        decay=decay,
        protein=protein,
    )
    return repeat_finder.mask_repeats(
        sequence,
        threshold=threshold,
        mask=mask,
    )
