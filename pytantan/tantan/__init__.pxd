cdef extern from "tantan.hh" namespace "tantan" nogil:

    void maskSequences(
        unsigned char* seqBeg,
        unsigned char* seqEnd,
        int maxRepeatOffset,
        const double** likelihoodRationMatrix,
        double repeatProb,
        double repeatEndProb,
        double repeatOffsetProbDecay,
        double firstGapProb,
        double otherGapProb,
        double minMaskProb,
        const unsigned char* maskTable
    ) except +


cdef extern from * nogil:
    """
    const size_t SCORE_MATRIX_SIZE = 64;
    """
    const size_t SCORE_MATRIX_SIZE

