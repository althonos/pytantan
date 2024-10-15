from libcpp cimport bool

cdef extern from "LambdaCalculator.hh" namespace "cbrc" nogil:

    cppclass LambdaCalculator:
        LambdaCalculator()
        
        void calculate(const int** matrix, int alphSize)
        void setBad()
        bool isBad()
        double lambda_ "lambda"()