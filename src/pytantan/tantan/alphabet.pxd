from libcpp.string cimport string


cdef extern from "mcf_alphabet.hh" namespace "mcf::Alphabet" nogil:
    const unsigned CAPACITY "mcf::Alphabet::capacity"
    const char* DNA "mcf::Alphabet::dna"
    const char* PROTEIN "mcf::Alphabet::protein"

cdef extern from "mcf_alphabet.hh" namespace "mcf" nogil:

    cppclass Alphabet:
        Alphabet()
        
        void fromString(const string& normalLetters)
        void makeCaseInsensitive()
        void encodeInPlace(unsigned char* begin, unsigned char* end)
        void decodeInPlace(unsigned char* begin, unsigned char* end)

        unsigned char size
        unsigned char[CAPACITY] lettersToNumbers
        unsigned char[CAPACITY] numbersToLetters
        unsigned char[CAPACITY] numbersToUppercase
        unsigned char[CAPACITY] numbersToLowercase
