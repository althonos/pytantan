from libcpp cimport bool

cdef extern from "mcf_tantan_options.hh" namespace "mcf::TantanOptions" nogil:

    enum OutputType:
        maskOut
        probOut
        countOut
        bedOut
        repOut

cdef extern from "mcf_tantan_options.hh" namespace "mcf" nogil:

    cppclass TantanOptions:
        bool isProtein
        char maskSymbol
        bool isPreserveLowercase
        const char *scoreMatrixFileName
        double repeatProb
        double repeatEndProb
        int maxCycleLength
        double repeatOffsetProbDecay
        int matchScore
        int mismatchCost
        int gapExistenceCost
        int gapExtensionCost
        double minMaskProb
        double minCopyNumber
        OutputType outputType

        TantanOptions()


