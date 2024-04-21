# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

from libc.math cimport exp
from libc.limits cimport INT_MAX
from libc.stdlib cimport calloc, free

from .tantan.options cimport TantanOptions, OutputType
from .tantan.alphabet cimport Alphabet as _Alphabet, DNA, PROTEIN
from .tantan.score_matrix cimport ScoreMatrix as _ScoreMatrix, BLOSUM62
from .tantan.utils cimport unstringify
from .tantan.lc cimport LambdaCalculator

# if NEON_BUILD_SUPPORT:
#     from pytantan.platform.neon cimport maskSequencesNEON
# if SSE2_BUILD_SUPPORT:
    # from .platform.sse2 cimport maskSequencesSSE2
if SSE4_BUILD_SUPPORT:
    from pytantan.platform.sse4 cimport maskSequencesSSE4
# if AVX2_BUILD_SUPPORT:
#     from pytantan.platform.avx2 cimport maskSequencesAVX2


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
    cdef          _ScoreMatrix scoring

    cdef double*  probMatrix
    cdef double** probMatrixPointers
    cdef int*     fastMatrix
    cdef int**    fastMatrixPointers

    def __init__(self, object alphabet not None):
        #
        if not isinstance(alphabet, Alphabet):
            alphabet = Alphabet(alphabet)
        self.alphabet = alphabet

        # make scoring matrix
        if alphabet.letters == PROTEIN.decode():
            unstringify(self.scoring, BLOSUM62)
        else:
            self.scoring.initMatchMismatch(1, 1, alphabet.letters.encode('ascii'))

        # allocate matrices
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

        # prepare pointers
        for i in range(64):
            self.fastMatrixPointers[i] = &self.fastMatrix[i * 64]
            self.probMatrixPointers[i] = &self.probMatrix[i * 64]
        for i in range(64*64):
            self.fastMatrix[i] = 0
            self.probMatrix[i] = 0.0

        # compute fast matrix
        self.scoring.makeFastMatrix(
            self.fastMatrixPointers,
            64,
            self.alphabet._abc.lettersToNumbers,
            self.scoring.minScore(),
            False
        )

        # compute lambda factor
        cdef LambdaCalculator lc = LambdaCalculator()
        lc.calculate(self.fastMatrixPointers, self.alphabet._abc.size)
        if lc.isBad():
            raise RuntimeError("cannot calculate probabilities for score matrix")

        # compute likelihood matrix
        cdef double matrixLambda = lc.lambda_()
        for i in range(64):
            for j in range(64):
                x = matrixLambda * self.fastMatrixPointers[i][j]
                self.probMatrixPointers[i][j] = exp(x)



cdef class Tantan:
    cdef          TantanOptions _options
    cdef readonly Alphabet      alphabet
    cdef readonly ScoreMatrix   score_matrix

    def __init__(
        self,
        ScoreMatrix score_matrix not None,
        *,
        bint protein = False,
    ):
        # store score matrix and alphabet
        self.score_matrix = score_matrix
        self.alphabet = score_matrix.alphabet
        # initialize options
        self._options.isProtein = protein
        if self._options.maxCycleLength < 0:
            self._options.maxCycleLength = 50 if self.alphabet._protein else 100
        if self._options.mismatchCost == 0:
            self._options.mismatchCost = INT_MAX

    def mask(self, object sequence):
        # extract sequence (FIXME)
        if isinstance(sequence, str):
            sequence = bytearray(sequence, 'ascii')
        elif not isinstance(sequence, bytearray):
            sequence = bytearray(sequence)

        cdef unsigned char[::1] seq = sequence

        self.alphabet._abc.encodeInPlace(&seq[0], &seq[seq.shape[0]])
        maskSequencesSSE4(
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