# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

from libc.math cimport exp
from libc.limits cimport INT_MAX
from libc.stdlib cimport calloc, free
from libcpp.vector cimport vector

from .tantan cimport SCORE_MATRIX_SIZE
from .tantan.options cimport TantanOptions, OutputType
from .tantan.alphabet cimport Alphabet as _Alphabet, DNA, PROTEIN
from .tantan.score_matrix cimport ScoreMatrix as _ScoreMatrix, BLOSUM62
from .tantan.utils cimport unstringify
from .tantan.lc cimport LambdaCalculator

from pytantan.platform.generic cimport maskSequencesNone
# if NEON_BUILD_SUPPORT:
#     from pytantan.platform.neon cimport maskSequencesNEON
# if SSE2_BUILD_SUPPORT:
    # from .platform.sse2 cimport maskSequencesSSE2
# if SSE4_BUILD_SUPPORT:
#     from pytantan.platform.sse4 cimport maskSequencesSSE4
# if AVX2_BUILD_SUPPORT:
#     from pytantan.platform.avx2 cimport maskSequencesAVX2

cdef extern from "<cctype>" namespace "std" nogil:
    cdef int toupper(int ch)
    cdef int tolower(int ch)


cdef class Alphabet:
    cdef _Alphabet _abc
    cdef bint      _protein

    @classmethod
    def dna(cls):
        return cls(DNA, protein=False)

    @classmethod
    def amino(cls):
        return cls(PROTEIN, protein=True)

    def __init__(self, object letters not None, bint protein = False):
        if isinstance(letters, str):
            letters = letters.encode('ascii')
        self._protein = protein
        self._abc.fromString(letters)

    @property
    def letters(self):
        cdef size_t i
        return ''.join(
            chr(self._abc.numbersToLetters[i])
            for i in range(self._abc.size)
        )


cdef class ScoreMatrix:
    cdef readonly Alphabet     alphabet

    cdef double*  probMatrix
    cdef double** probMatrixPointers
    cdef int*     fastMatrix
    cdef int**    fastMatrixPointers

    cdef int _allocate_matrix(self) except 1:
        # allocate data and pointer arrays
        self.probMatrix = <double*> calloc(sizeof(double), 64*64)
        if self.probMatrix is NULL:
            raise MemoryError
        self.probMatrixPointers = <double**> calloc(sizeof(double*), 64)
        if self.probMatrixPointers is NULL:
            raise MemoryError
        self.fastMatrix = <int*> calloc(sizeof(int), 64*64)
        if self.fastMatrix is NULL:
            raise MemoryError
        self.fastMatrixPointers = <int**> calloc(sizeof(int*), 64)
        if self.fastMatrixPointers is NULL:
            raise MemoryError
        # prepare pointers to data
        for i in range(64):
            self.fastMatrixPointers[i] = &self.fastMatrix[i * 64]
            self.probMatrixPointers[i] = &self.probMatrix[i * 64]
        return 0

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
        for i in range(SCORE_MATRIX_SIZE):
            for j in range(SCORE_MATRIX_SIZE):
                self.fastMatrixPointers[i][j] = default

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

    def __init__(self, Alphabet alphabet not None, object matrix not None):
        cdef double              lambda_  
        cdef vector[vector[int]] scores  = vector[vector[int]]()
        cdef LambdaCalculator    lc      = LambdaCalculator()

        # record alphabet
        self.alphabet = alphabet

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

        # compute lambda factor
        lc.calculate(self.fastMatrixPointers, self.alphabet._abc.size)
        if lc.isBad():
            raise RuntimeError("cannot calculate probabilities for score matrix")

        # compute likelihood matrix
        lambda_ = lc.lambda_()
        for i in range(64):
            for j in range(64):
                x = lambda_ * self.fastMatrixPointers[i][j]
                self.probMatrixPointers[i][j] = exp(x)


cdef class Tantan:
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

    def mask(self, object sequence):
        # extract sequence (FIXME)
        if isinstance(sequence, str):
            sequence = bytearray(sequence, 'ascii')
        elif not isinstance(sequence, bytearray):
            sequence = bytearray(sequence)

        cdef unsigned char[::1] seq = sequence

        with nogil:
            self.alphabet._abc.encodeInPlace(&seq[0], &seq[seq.shape[0]])
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
                self._options.minMaskProb,
                self.alphabet._abc.numbersToLowercase,
            )
            self.alphabet._abc.decodeInPlace(&seq[0], &seq[seq.shape[0]])

        return sequence.decode()