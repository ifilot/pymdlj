# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr

cdef extern from "cell.cpp":
    pass

cdef extern from "parameters.cpp":
    pass

cdef extern from "parameters.h":
    cdef cppclass Parameters:
        Parameters(string input) except +
        double get_param[double](string varname) except +
        uint get_param[uint](string varname) except +

cdef extern from "cell.h":
    cdef cppclass Cell:
        Cell(shared_ptr[Parameters]) except +
        void integrate(double, double) except +
