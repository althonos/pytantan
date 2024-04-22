# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

# --- C imports ----------------------------------------------------------------

from libc.math cimport exp
from libc.limits cimport INT_MAX, UCHAR_MAX
from libc.stdlib cimport calloc, free
from libc.string cimport memcmp
from libcpp.vector cimport vector

from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBUF_WRITE

from .tantan cimport SCORE_MATRIX_SIZE
from .tantan.options cimport TantanOptions, OutputType
from .tantan.alphabet cimport Alphabet as _Alphabet, DNA, PROTEIN
from .tantan.score_matrix cimport ScoreMatrix as _ScoreMatrix, BLOSUM62
from .tantan.utils cimport unstringify
from .tantan.lc cimport LambdaCalculator

from pytantan.platform.generic cimport maskSequencesNone, getProbabilitiesNone
# if NEON_BUILD_SUPPORT:
#     from pytantan.platform.neon cimport maskSequencesNEON
# if SSE2_BUILD_SUPPORT:
    # from .platform.sse2 cimport maskSequencesSSE2
# if SSE4_BUILD_SUPPORT:
#     from pytantan.platform.sse4 cimport maskSequencesSSE4
# if AVX2_BUILD_SUPPORT:
#     from pytantan.platform.avx2 cimport maskSequencesAVX2

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

# --- Parameters ---------------------------------------------------------------

cdef class Alphabet:
    cdef readonly str       letters
    cdef          bint      _protein
    cdef          _Alphabet _abc

    @classmethod
    def dna(cls):
        return cls(DNA.decode('ascii'), protein=False)

    @classmethod
    def amino(cls):
        return cls(PROTEIN.decode('ascii'), protein=True)

    def __init__(self, str letters not None, bint protein = False):
        self.letters = letters
        self._protein = protein
        self._abc.fromString(self.letters.encode('ascii'))

    def __len__(self):
        return self.length

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


cdef class ScoreMatrix:
    cdef readonly Alphabet alphabet
    cdef Py_ssize_t        _shape[2]

    cdef double** probMatrixPointers
    cdef int**    fastMatrixPointers

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

    cdef int _make_fast_matrix(self, vector[vector[int]]& scores) except 1:
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
        for i in range(self.alphabet._abc.size):
            for j in range(self.alphabet._abc.size):
                if scores[i][j] < default:
                    default = scores[i][j]

        # fill matrix with default score (FIXME)
        for i in range(SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE):
            self.fastMatrixPointers[0][i] = default

        # update matrix with scores
        for i in range(self.alphabet._abc.size):
            x = self.alphabet._abc.numbersToLetters[i]
            a = self.alphabet._abc.lettersToNumbers[toupper(x)]
            b = self.alphabet._abc.lettersToNumbers[tolower(x)]

            for j in range(self.alphabet._abc.size):
                y = self.alphabet._abc.numbersToLetters[j]
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

    def __init__(self, Alphabet alphabet not None, object matrix not None):
        cdef double              lambda_
        cdef vector[vector[int]] scores  = vector[vector[int]]()

        # record alphabet
        self.alphabet = alphabet
        self._shape[0] = self._shape[1] = self.alphabet._abc.size

        # extract scores
        if len(matrix) != self.alphabet._abc.size:
            raise ValueError("Matrix length should be equal to alphabet length")
        scores.resize(alphabet._abc.size)
        for i, row in enumerate(matrix):
            scores[i].resize(alphabet._abc.size)
            if len(row) != alphabet._abc.size:
                raise ValueError("Matrix width should be equal to alphabet length")
            for j, x in enumerate(row):
                scores[i][j] = x

        # prepare fast matrix
        self._allocate_matrix()
        self._make_fast_matrix(scores)
        self._make_likelihood_matrix()

    def __dealloc__(self):
        if self.probMatrixPointers is not NULL:
            free(self.probMatrixPointers[0])
        free(self.probMatrixPointers)
        if self.fastMatrixPointers is not NULL:
            free(self.fastMatrixPointers[0])
        free(self.fastMatrixPointers)

    def __eq__(self, object other):
        cdef ScoreMatrix matrix
        cdef int         result

        if not isinstance(other, ScoreMatrix):
            return NotImplemented
        matrix = other
        if self.alphabet != other.alphabet:
            return False
        with nogil:
            result = memcmp(
                self.fastMatrixPointers[0],
                matrix.fastMatrixPointers[0],
                SCORE_MATRIX_SIZE*SCORE_MATRIX_SIZE*sizeof(int)
            )
        return result == 0

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.alphabet!r}, {self.matrix!r})"

    def __reduce__(self):
        return (type(self), (self.alphabet, self.matrix))

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        if flags & PyBUF_FORMAT:
            buffer.format = b"i"
        else:
            buffer.format = NULL
        buffer.buf = self.fastMatrixPointers[0]
        buffer.internal = NULL
        buffer.itemsize = sizeof(int)
        buffer.len = self._shape[0] * self._shape[1] * sizeof(int)
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = <Py_ssize_t*> &self._shape
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def matrix(self):
        """`list` of `list` of `int`: The score matrix.
        """
        cdef int        i
        cdef int        j
        cdef int        length = self.alphabet._abc.size

        return [
            [ self.fastMatrixPointers[i][j] for j in range(length) ]
            for i in range(length)
        ]


# --- RepeatFinder -------------------------------------------------------------

cdef class RepeatFinder:
    cdef          TantanOptions _options
    cdef readonly Alphabet      alphabet
    cdef readonly ScoreMatrix   score_matrix

    def __init__(
        self,
        ScoreMatrix score_matrix not None,
        *,
        double repeat_start = 0.005,
        double repeat_end = 0.05,
        double decay = 0.9,
    ):
        # store score matrix and alphabet
        self.score_matrix = score_matrix
        self.alphabet = score_matrix.alphabet
        # initialize options
        self._options.isProtein = self.alphabet._protein
        if self._options.maxCycleLength < 0:
            self._options.maxCycleLength = 50 if self.alphabet._protein else 100
        if self._options.mismatchCost == 0:
            self._options.mismatchCost = INT_MAX
        self._options.repeatProb = repeat_start
        self._options.repeatEndProb = repeat_end
        self._options.repeatOffsetProbDecay = decay

    cpdef object get_probabilities(self, object sequence):
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
                <const double**> self.score_matrix.probMatrixPointers,
                self._options.repeatProb,
                self._options.repeatEndProb,
                self._options.repeatOffsetProbDecay,
                0.0, #firstGapProb,
                0.0, #otherGapProb,
                &p[0],
            )
        return probas

    cpdef str mask_repeats(
        self,
        object sequence,
        double threshold = 0.5,
        object mask = None,
    ):
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
                <const double**> self.score_matrix.probMatrixPointers,
                self._options.repeatProb,
                self._options.repeatEndProb,
                self._options.repeatOffsetProbDecay,
                0.0, #firstGapProb,
                0.0, #otherGapProb,
                threshold,
                mask_ptr,
            )
        return self.alphabet.decode(seq)