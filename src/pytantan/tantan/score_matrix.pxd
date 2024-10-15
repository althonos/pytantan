from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "mcf_score_matrix.hh" namespace "mcf" nogil:

    cppclass ScoreMatrix:
        ScoreMatrix()

        void initMatchMismatch(int matchScore, int mismatchCost, const string& letters) except +
        int minScore() except +
        void makeFastMatrix(int** fastMatrix, int fastMatrixSize, const unsigned char* letterToIndex, int defaultScore, bool isCaseSensitive) except +


cdef extern from "mcf_score_matrix.hh" namespace "mcf::ScoreMatrix" nogil:
    const char* BLOSUM62 "mcf::ScoreMatrix::blosum62"