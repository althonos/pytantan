# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

from libc.math cimport exp
from libc.limits cimport INT_MAX
from libc.stdlib cimport calloc, free

from .tantan.options cimport TantanOptions, OutputType
from .tantan.alphabet cimport Alphabet as _Alphabet, DNA, PROTEIN
from .tantan.score_matrix cimport ScoreMatrix as _ScoreMatrix
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

    @classmethod
    def dna(cls):
        cdef Alphabet abc = cls.__new__(cls)
        abc._abc.fromString(DNA)
        return abc

    @classmethod
    def protein(cls):
        cdef Alphabet abc = cls.__new__(cls)
        abc._abc.fromString(PROTEIN)
        return abc


cdef class ScoreMatrix:
    cdef readonly Alphabet     alphabet
    cdef          _ScoreMatrix scoring

    cdef double* probMatrix
    cdef double** probMatrixPointers
    cdef int* fastMatrix
    cdef int** fastMatrixPointers

    def __init__(self, Alphabet alphabet not None):
        # 
        self.alphabet = alphabet      

        # make scoring matrix
        self.scoring.initMatchMismatch(1, 1, b"ACGTU")

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

        for i in range(64):
            self.fastMatrixPointers[i] = &self.fastMatrix[i * 64]
            self.probMatrixPointers[i] = &self.probMatrix[i * 64]
        for i in range(64*64):
            self.fastMatrix[i] = 0
            self.probMatrix[i] = 0.0

        self.scoring.makeFastMatrix(
            self.fastMatrixPointers, 
            64, 
            self.alphabet._abc.lettersToNumbers, 
            self.scoring.minScore(), 
            False
        )

        cdef LambdaCalculator lc = LambdaCalculator()
        lc.calculate(self.fastMatrixPointers, self.alphabet._abc.size)
        if lc.isBad():
            raise RuntimeError("cannot calculate probabilities for score matrix")
        
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
            self._options.maxCycleLength = 50 if protein else 100
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
            100,
            <const double**> self.score_matrix.probMatrixPointers,
            0.005,
            0.05,
            0.9,
            0.0, #firstGapProb,
            0.0, #otherGapProb,
            0.5,
            self.alphabet._abc.numbersToLowercase,
        )

        self.alphabet._abc.decodeInPlace(&seq[0], &seq[seq.shape[0]])
        return sequence.decode()