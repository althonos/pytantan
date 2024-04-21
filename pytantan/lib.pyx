# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True

from libc.math cimport exp
from libc.stdlib cimport calloc, free

from .tantan.options cimport TantanOptions, OutputType
from .tantan.alphabet cimport Alphabet, DNA
from .tantan.score_matrix cimport ScoreMatrix, BLOSUM62
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


cdef class Tantan:

    def mask(self, object sequence):
        # get default options (FIXME)
        cdef TantanOptions options = TantanOptions()

        # make alphabet (FIXME)
        cdef Alphabet alphabet = Alphabet()
        alphabet.fromString(DNA)
        
        # extract sequence (FIXME)
        if isinstance(sequence, str):
            sequence = bytearray(sequence, 'ascii')
        elif not isinstance(sequence, bytearray):
            sequence = bytearray(sequence)
        cdef unsigned char[::1] seq = sequence

        # make scoring matrix
        cdef ScoreMatrix   scoring = ScoreMatrix()
        scoring.initMatchMismatch(1, 1, b"ACGTU")

        cdef double* probMatrix = <double*> calloc(sizeof(double), 64*64)
        if probMatrix is NULL:
            raise MemoryError
        cdef double** probMatrixPointers = <double**> calloc(sizeof(double*), 64)
        if probMatrixPointers is NULL:
            raise MemoryError      
        cdef int* fastMatrix = <int*> calloc(sizeof(int), 64*64)
        if fastMatrix is NULL:
            raise MemoryError
        cdef int** fastMatrixPointers = <int**> calloc(sizeof(int*), 64)
        if fastMatrixPointers is NULL:
            raise MemoryError

        for i in range(64):
            fastMatrixPointers[i] = &fastMatrix[i * 64]
            probMatrixPointers[i] = &probMatrix[i * 64]
        for i in range(64*64):
            fastMatrix[i] = 0
            probMatrix[i] = 0.0

        scoring.makeFastMatrix(
            fastMatrixPointers, 
            64, 
            alphabet.lettersToNumbers, 
            scoring.minScore(), 
            False
        )

        for i in range(alphabet.size):
            for j in range(alphabet.size):
                print(fastMatrixPointers[i][j], end=" ")
            print()

        cdef LambdaCalculator lc = LambdaCalculator()
        lc.calculate(fastMatrixPointers, alphabet.size)
        if lc.isBad():
            raise RuntimeError("cannot calculate probabilities for score matrix")
        
        cdef double matrixLambda = lc.lambda_()
        for i in range(64):
            for j in range(64):
                x = matrixLambda * fastMatrixPointers[i][j]
                if options.outputType != OutputType.repOut:
                    x = exp(x)
                probMatrixPointers[i][j] = x

        alphabet.encodeInPlace(&seq[0], &seq[seq.shape[0]])
        maskSequencesSSE4(
            &seq[0],
            &seq[seq.shape[0]],
            100,
            <const double**> probMatrixPointers,
            0.005,
            0.05,
            0.9,
            0.0, #firstGapProb,
            0.0, #otherGapProb,
            0.5,
            alphabet.numbersToLowercase,
        )

        alphabet.decodeInPlace(&seq[0], &seq[seq.shape[0]])
        return sequence.decode()