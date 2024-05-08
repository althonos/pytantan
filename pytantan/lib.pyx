# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True
"""Bindings to tantan, a method for finding repeats in biological sequences.

References:
    - Frith, Martin C. *A new repeat-masking method enables specific detection 
      of homologous sequences.* Nucleic acids research vol. 39,4 (2011): e23. 
      :doi:`10.1093/nar/gkq1212`.

"""

# --- C imports ----------------------------------------------------------------

cimport cython
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBUF_WRITE

from libc.math cimport exp
from libc.limits cimport INT_MAX, UCHAR_MAX
from libc.stdlib cimport calloc, free
from libc.string cimport memcmp
from libcpp.vector cimport vector

from .tantan cimport SCORE_MATRIX_SIZE
from .tantan.options cimport TantanOptions, OutputType
from .tantan.alphabet cimport Alphabet as _Alphabet, DNA, PROTEIN
from .tantan.score_matrix cimport ScoreMatrix as _ScoreMatrix, BLOSUM62
from .tantan.utils cimport unstringify
from .tantan.lc cimport LambdaCalculator

from .platform cimport _mask_fn_t, _probas_fn_t
from .platform.generic cimport maskSequencesNone, getProbabilitiesNone
if NEON_BUILD_SUPPORT:
    from .platform.neon cimport maskSequencesNEON, getProbabilitiesNEON
# if SSE2_BUILD_SUPPORT:
#     from .platform.sse2 cimport maskSequencesSSE2, getProbabilitiesSSE2
if SSE4_BUILD_SUPPORT:
    from .platform.sse4 cimport maskSequencesSSE4, getProbabilitiesSSE4
if AVX2_BUILD_SUPPORT:
    from .platform.avx2 cimport maskSequencesAVX2, getProbabilitiesAVX2

from scoring_matrices cimport ScoringMatrix

cdef extern from "<cctype>" namespace "std" nogil:
    cdef bint isalpha(int ch)
    cdef int toupper(int ch)
    cdef int tolower(int ch)

# --- Python imports -----------------------------------------------------------

import array
import itertools

# --- Runtime CPU detection ----------------------------------------------------

try:
    import archspec.cpu
    _HOST_CPU             = archspec.cpu.host()
    _HOST_FEATURES        = _HOST_CPU.features
except ImportError:
    _HOST_CPU             = None
    _HOST_FEATURES        = set()

_SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
_SSE4_BUILD_SUPPORT   = SSE4_BUILD_SUPPORT
_AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
_NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
_SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT and "sse2" in _HOST_FEATURES
_SSE4_RUNTIME_SUPPORT = SSE4_BUILD_SUPPORT and "sse4_1" in _HOST_FEATURES
_AVX2_RUNTIME_SUPPORT = AVX2_BUILD_SUPPORT and "avx2" in _HOST_FEATURES
_NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT and "neon" in _HOST_FEATURES

# NOTE(@althonos): NEON is always supported on Aarch64 so we should only check
#                  that the extension was built with NEON support.
if TARGET_CPU == "aarch64":
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT

# NOTE(@althonos): SSE2 is always supported on x86-64 so we should only check
#                  that the extension was built with SSE2 support.
if TARGET_CPU == "x86_64":
    _SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT


# --- Utils --------------------------------------------------------------------

cpdef ScoringMatrix default_scoring_matrix(
    bint protein = False,
    object match_score = None,
    object mismatch_cost = None
):
    """Get the default Tantan scoring matrix for the given parameters.
    """
    if (match_score is None) != (mismatch_cost is None):
        raise ValueError("Cannot set either `match_score` or `mismatch_cost` without setting the two")
    if match_score is not None:
        alphabet = PROTEIN.decode('ascii') if protein else DNA.decode('ascii')
        return ScoringMatrix.from_match_mismatch(
            alphabet=alphabet,
            match_score=match_score,
            mismatch_score=-mismatch_cost,
        )
    elif protein:
        return ScoringMatrix.from_name("BLOSUM62").shuffle(PROTEIN.decode('ascii'))
    else:
        return ScoringMatrix.from_match_mismatch(1.0, -1.0, alphabet=DNA.decode('ascii'))


# --- Parameters ---------------------------------------------------------------

cdef class Alphabet:
    """An alphabet used for encoding sequences with ordinal encoding.
    """
    cdef readonly str       letters
    cdef          bint      _protein
    cdef          _Alphabet _abc

    @classmethod
    def dna(cls):
        return cls(DNA.decode('ascii'), protein=False)

    @classmethod
    def protein(cls):
        return cls(PROTEIN.decode('ascii'), protein=True)

    def __init__(self, str letters not None, bint protein = False):
        if len(letters) != len(set(letters)):
            raise ValueError(f"Duplicate letters found in alphabet: {letters!r}")
        elif not letters.isupper():
            raise ValueError("Alphabet must only contain uppercase characters")
        elif not letters.isalpha():
            raise ValueError("Alphabet must only contain alphabetic characters")
        self.letters = letters
        self._protein = protein
        self._abc.fromString(self.letters.encode('ascii'))

    def __len__(self):
        return self._abc.size

    def __contains__(self, object item):
        return item in self.letters

    def __getitem__(self, ssize_t index):
        cdef ssize_t index_ = index
        if index_ < 0:
            index_ += self._abc.size
        if index_ < 0 or index_ >= self._abc.size:
            raise IndexError(index)
        return chr(self._abc.numbersToLetters[index_])

    def __reduce__(self):
        return type(self), (self.letters,)

    def __repr__(self):
        return f"{type(self).__name__}({self.letters!r})"

    def __str__(self):
        return self.letters

    def __eq__(self, object item):
        if isinstance(item, str):
            return self.letters == item
        elif isinstance(item, Alphabet):
            return self.letters == item.letters
        else:
            return False

    cpdef void encode_into(self, const unsigned char[:] sequence, unsigned char[:] encoded):
        r"""Encode a sequence to ordinal-encoding into the given buffer.
        """
        cdef ssize_t       i
        cdef unsigned char code
        cdef unsigned char letter

        if sequence.shape[0] != encoded.shape[0]:
            raise ValueError("Buffers do not have the same dimensions")

        with nogil:
            for i in range(sequence.shape[0]):
                letter = sequence[i]
                if not isalpha(letter):
                    raise ValueError(f"Character outside ASCII range: {chr(letter)!r}")
                code = self._abc.lettersToNumbers[<unsigned char> letter]
                if code == 255:
                    raise ValueError(f"Non-alphabet character in sequence: {chr(letter)!r}")
                encoded[i] = code

    cpdef void decode_into(self, const unsigned char[:] encoded, unsigned char[:] sequence):
        r"""Decode a sequence from ordinal-encoding into the given buffer.
        """
        cdef unsigned char code
        cdef size_t        length = len(self.letters)

        if sequence.shape[0] != encoded.shape[0]:
            raise ValueError("Buffers do not have the same dimensions")

        with nogil:
            for i in range(encoded.shape[0]):
                code = encoded[i]
                sequence[i] = self._abc.numbersToLetters[code]

    cpdef bytes encode(self, object sequence):
        r"""Encode a sequence to an ordinal-encoded sequence using the alphabet.

        Arguments:
            sequence (`str` or byte-like object): The sequence to encode.

        Returns:
            `bytes`: The encoded sequence.

        Raises:
            `ValueError`: When the sequence contains invalid characters, or
            unknown sequence characters while the alphabet contains no
            wildcard character.

        Example:
            >>> alphabet = Alphabet("ACGT")
            >>> alphabet.encode("GATACA")
            b'\x02\x00\x03\x00\x01\x00'

        """
        cdef bytearray encoded = bytearray(len(sequence))
        if isinstance(sequence, str):
            sequence = sequence.encode('ascii')
        self.encode_into(sequence, encoded)
        return bytes(encoded)

    cpdef str decode(self, object encoded):
        r"""Decode an ordinal-encoded sequence using the alphabet.

        Arguments:
            sequence (byte-like object): The sequence to decode.

        Returns:
            `str`: The decoded sequence.

        Raises:
            `ValueError`: When the sequence contains invalid indices.

        Example:
            >>> alphabet = Alphabet("ACGT")
            >>> alphabet.decode(bytearray([2, 0, 3, 0, 1, 0]))
            'GATACA'

        """
        cdef bytearray decoded = bytearray(len(encoded))
        self.decode_into(encoded, decoded)
        return decoded.decode('ascii')


cdef class LikelihoodMatrix:
    r"""A likelihood ratio matrix derived from a scoring matrix.

    This matrix is in the form:

    .. math::

        \mathcal{L}_{x,y} = \frac{Q_{x,y}}{P_x P_y}

    where :math:`Q_{x,y}` is the probability of seeing letters :math:`x` and
    :math:`y` in equivalent positions of a tandem repeat, and :math:`P_x` is
    the background probability of letter :math:`x`.

    This matrix is related to a scoring matrix :math:`\mathcal{M}` with an
    implicit scale factor :math:`\lambda` by the following equation:

    .. math::

        \mathcal{L}_{x,y} = e^{ \lambda \mathcal{M}_{x,y} }

    Attributes:
        alphabet (`Alphabet`): The alphabet used to index the rows and
            columns of the matrix.
        scoring_matrix (`~scoring_matrices.ScoringMatrix` or `None`): The
            scoring matrix this likelihood ratio matrix was derived from,
            or `None` if it was created with explicit values.

    """
    cdef Py_ssize_t             _shape[2]

    cdef readonly Alphabet      alphabet
    cdef readonly ScoringMatrix scoring_matrix

    cdef int**    fastMatrixPointers
    cdef double** probMatrixPointers

    cdef int _allocate_matrix(self) except 1:
        cdef size_t i

        # allocate data and pointer arrays
        self.probMatrixPointers = <double**> calloc(sizeof(double*), 64)
        if self.probMatrixPointers is NULL:
            raise MemoryError
        self.probMatrixPointers[0] = <double*> calloc(sizeof(double), 64*64)
        if self.probMatrixPointers[0] is NULL:
            raise MemoryError
        self.fastMatrixPointers = <int**> calloc(sizeof(int*), 64)
        if self.fastMatrixPointers is NULL:
            raise MemoryError
        self.fastMatrixPointers[0] = <int*> calloc(sizeof(int), 64*64)
        if self.fastMatrixPointers[0] is NULL:
            raise MemoryError

        # prepare pointers to data
        for i in range(64):
            self.fastMatrixPointers[i] = &self.fastMatrixPointers[0][i * 64]
            self.probMatrixPointers[i] = &self.probMatrixPointers[0][i * 64]

    cdef int _make_fast_matrix(self, vector[vector[int]]& scores, vector[char]& letters) except 1:
        cdef size_t i
        cdef size_t j
        cdef int    a
        cdef int    b
        cdef int    c
        cdef int    d
        cdef int    x
        cdef int    y
        cdef int    default = INT_MAX

        # find default score
        for i in range(scores.size()):
            for j in range(scores.size()):
                if j >= scores[i].size():
                    raise IndexError(j)
                if scores[i][j] < default:
                    default = scores[i][j]

        # fill matrix with default score (FIXME)
        for i in range(SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE):
            self.fastMatrixPointers[0][i] = default

        # update matrix with scores
        for i, x in enumerate(letters):
            a = self.alphabet._abc.lettersToNumbers[toupper(x)]
            b = self.alphabet._abc.lettersToNumbers[tolower(x)]

            for j, y in enumerate(letters):
                c = self.alphabet._abc.lettersToNumbers[toupper(y)]
                d = self.alphabet._abc.lettersToNumbers[tolower(y)]

                self.fastMatrixPointers[a][c] = scores[i][j]
                self.fastMatrixPointers[a][d] = scores[i][j]
                self.fastMatrixPointers[b][c] = scores[i][j]
                self.fastMatrixPointers[b][d] = scores[i][j]

    cdef int _make_likelihood_matrix(self):
        cdef size_t           i
        cdef LambdaCalculator lc = LambdaCalculator()

        # compute lambda factor
        lc.calculate(self.fastMatrixPointers, self.alphabet._abc.size)
        if lc.isBad():
            raise RuntimeError("cannot calculate probabilities for score matrix")

        # compute likelihood matrix
        lambda_ = lc.lambda_()
        for i in range(SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE):
            x = lambda_ * self.fastMatrixPointers[0][i]
            self.probMatrixPointers[0][i] = exp(x)

    def __init__(
        self,
        Alphabet alphabet not None,
        ScoringMatrix scoring_matrix not None,
    ):
        """Create a new score matrix from the given alphabet and scores.

        Arguments:
            alphabet (`str` or `~pytantan.Alphabet`): The alphabet of the
                similarity matrix.
            scoring_matrix (`~scoring_matrices.ScoringMatrix`): The scoring
                matrix, indexed by the alphabet characters.

        """
        cdef size_t              i
        cdef size_t              colsize
        cdef double              lambda_
        cdef const float**       matrix
        cdef size_t              size
        cdef vector[char]        cols    = vector[char]()
        cdef vector[vector[int]] scores  = vector[vector[int]]()

        # record alphabet
        self.alphabet = alphabet
        self._shape[0] = self._shape[1] = SCORE_MATRIX_SIZE

        # handle scoring matrix
        if not scoring_matrix.is_integer():
            raise ValueError("Expected integer scoring matrix")
        self.scoring_matrix = scoring_matrix
        for l in self.scoring_matrix.alphabet:
            cols.push_back(ord(l))

        # extract scores
        with nogil:
            matrix = self.scoring_matrix.matrix_ptr()
            size = self.scoring_matrix.size()
            for i in range(size):
                scores.push_back(vector[int]())
                for j in range(size):
                    scores[i].push_back(<int> matrix[i][j])

        # prepare fast matrix
        self._allocate_matrix()
        self._make_fast_matrix(scores, cols)
        self._make_likelihood_matrix()

    def __dealloc__(self):
        if self.probMatrixPointers is not NULL:
            free(self.probMatrixPointers[0])
        free(self.probMatrixPointers)
        if self.fastMatrixPointers is not NULL:
            free(self.fastMatrixPointers[0])
        free(self.fastMatrixPointers)


# --- RepeatFinder -------------------------------------------------------------

cdef class RepeatFinder:
    """A repeat finder using the Tantan method.

    Attributes:
        alphabet (`Alphabet`): The alphabet to use for encoding and decoding
            sequences.
        scoring_matrix (`~scoring_matrices.ScoringMatrix` or `None`): The
            scoring matrix used to derive the likelihood ratio matrix,
            if any.
        likelihood_matrix (`LikelihoodMatrix`): The likelihood ratio matrix
            used for scoring letter pairs in tandem repeats.

    """
    cdef          TantanOptions _options
    cdef          _mask_fn_t    _mask_sequences
    cdef          _probas_fn_t  _get_probabilities

    cdef readonly Alphabet         alphabet
    cdef readonly ScoringMatrix    scoring_matrix
    cdef readonly LikelihoodMatrix likelihood_matrix

    def __cinit__(self):
        self._mask_sequences = NULL
        self._get_probabilities = NULL

    def __init__(
        self,
        ScoringMatrix scoring_matrix not None,
        *,
        double repeat_start = 0.005,
        double repeat_end = 0.05,
        double decay = 0.9,
        bint protein = False,
    ):
        """Create a new repeat finder.

        Arguments:
            scoring_matrix (`~scoring_matrices.ScoringMatrix`): The scoring 
                matrix to use for scoring sequence alignments.
            repeat_start (`float`): The probability of a repeat starting
                per position.
            repeat_end (`float`): The probability of a repeat ending per
                position.
            decay (`float`): The probability decay per period.
            protein (`bool`): Set to `True` to treat the input sequence as
                a protein sequence.

        """
        # store score matrix and alphabet
        self.alphabet = Alphabet(scoring_matrix.alphabet)
        self.scoring_matrix = scoring_matrix
        self.likelihood_matrix = LikelihoodMatrix(self.alphabet, self.scoring_matrix)
        # initialize options
        self._options.isProtein = protein
        if self._options.maxCycleLength < 0:
            self._options.maxCycleLength = 50 if protein else 100
        if self._options.mismatchCost == 0:
            self._options.mismatchCost = INT_MAX
        self._options.repeatProb = repeat_start
        self._options.repeatEndProb = repeat_end
        self._options.repeatOffsetProbDecay = decay
        # Select the best available SIMD backend
        # if SSE2_BUILD_SUPPORT and _SSE2_RUNTIME_SUPPORT:
        #     self._get_probabilities = getProbabilitiesSSE2
        #     self._mask_sequences = maskSequencesSSE2
        if SSE4_BUILD_SUPPORT and _SSE4_RUNTIME_SUPPORT:
            self._get_probabilities = getProbabilitiesSSE4
            self._mask_sequences = maskSequencesSSE4
        if AVX2_BUILD_SUPPORT and _AVX2_RUNTIME_SUPPORT:
            self._get_probabilities = getProbabilitiesAVX2
            self._mask_sequences = maskSequencesAVX2
        if NEON_BUILD_SUPPORT and _NEON_RUNTIME_SUPPORT:
            self._get_probabilities = getProbabilitiesNEON
            self._mask_sequences = maskSequencesNEON
        if self._get_probabilities is NULL:
            self._get_probabilities = getProbabilitiesNone
        if self._mask_sequences is NULL:
            self._mask_sequences = maskSequencesNone

    @cython.boundscheck(False)
    cpdef object get_probabilities(self, object sequence):
        """Get the probabilities of being a repeat for each sequence position.

        Arguments:
            sequence (`str` or byte-like object): The sequence to compute
                probabilities for.

        Returns:
            `~array.array`: A float array containing the repeat probability
            per sequence position.

        """
        cdef unsigned char[::1] seq
        # extract sequence (FIXME)
        if isinstance(sequence, str):
            seq = bytearray(sequence, 'ascii')
        elif not isinstance(sequence, bytearray):
            seq = bytearray(sequence)
        else:
            seq = sequence.copy()
        # prepare probability array
        cdef object     probas = array.array("f", itertools.repeat(0.0, seq.shape[0]))
        cdef float[::1] p      = probas
        # compute probabilities
        self.alphabet.encode_into(seq, seq)
        with nogil:
            getProbabilitiesNone(
                &seq[0],
                &seq[seq.shape[0]],
                self._options.maxCycleLength,
                <const double**> self.likelihood_matrix.probMatrixPointers,
                self._options.repeatProb,
                self._options.repeatEndProb,
                self._options.repeatOffsetProbDecay,
                0.0, #firstGapProb,
                0.0, #otherGapProb,
                &p[0],
            )
        return probas

    @cython.boundscheck(False)
    cpdef str mask_repeats(
        self,
        object sequence,
        double threshold = 0.5,
        object mask = None,
    ):
        """Mask regions predicted as repeats in the given sequence.

        Arguments:
            sequence (`str` or byte-like object): The sequence containing
                the repeats to mask.
            threshold (`float`): The probability threshold above which to
                mask sequence characters.
            mask (`str` or `None`): A single mask character to use for
                masking positions. If `None` given, masking uses the
                lowercase letters of the original sequence.

        Returns:
            `str`: The input sequence with repeat regions masked.

        Example:
            >>> matrix = ScoringMatrix(
            ...     [[ 1, -1, -1, -1],
            ...      [-1,  1, -1, -1],
            ...      [-1, -1,  1, -1],
            ...      [-1, -1, -1,  1]],
            ...     alphabet="ACGT",
            ... )
            >>> tantan = RepeatFinder(matrix)
            >>> tantan.mask_repeats("ATTATTATTATTATT")
            'ATTattattattatt'
            >>> tantan.mask_repeats("ATTATTATTATTATT", mask='N')
            'ATTNNNNNNNNNNNN'

        """
        cdef unsigned char[::1]    seq
        cdef vector[unsigned char] mask_data
        cdef unsigned char*        mask_ptr
        # extract sequence
        if isinstance(sequence, str):
            seq = bytearray(sequence, 'ascii')
        elif not isinstance(sequence, bytearray):
            seq = bytearray(sequence)
        else:
            seq = sequence.copy()
        # build mask
        if mask is not None:
            if ord(mask) >= UCHAR_MAX:
                raise ValueError(f"Invalid mask symbol: {mask!r}")
            mask_char = self.alphabet._abc.lettersToNumbers[ord(mask)]
            if mask_char == UCHAR_MAX:
                raise ValueError(f"Invalid mask symbol: {mask!r}")
            mask_data.resize(UCHAR_MAX, mask_char)
            mask_ptr = mask_data.data()
        else:
            mask_ptr = self.alphabet._abc.numbersToLowercase
        # mask sequence
        self.alphabet.encode_into(seq, seq)
        with nogil:
            maskSequencesNone(
                &seq[0],
                &seq[seq.shape[0]],
                self._options.maxCycleLength,
                <const double**> self.likelihood_matrix.probMatrixPointers,
                self._options.repeatProb,
                self._options.repeatEndProb,
                self._options.repeatOffsetProbDecay,
                0.0, #firstGapProb,
                0.0, #otherGapProb,
                threshold,
                mask_ptr,
            )
        return self.alphabet.decode(seq)