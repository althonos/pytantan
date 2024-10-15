from libcpp.string cimport string

cdef extern from "mcf_util.hh" namespace "mcf" nogil:
    void unstringify[T](T& x, string& s) except +