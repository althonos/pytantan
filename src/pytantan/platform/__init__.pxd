ctypedef int (*_mask_fn_t)(
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
) except 1 nogil

ctypedef int (*_probas_fn_t)(
    const unsigned char* seqBeg,
    const unsigned char* seqEnd,
    int maxRepeatOffset,
    const double** likelihoodRationMatrix,
    double repeatProb,
    double repeatEndProb,
    double repeatOffsetProbDecay,
    double firstGapProb,
    double otherGapProb,
    float* probabilities
) except 1 nogil