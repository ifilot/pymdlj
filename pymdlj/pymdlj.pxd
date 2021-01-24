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
        double integrate(double, double) except +
        double get_etot() except +
        double get_ekin() except +
        double get_epot() except +
        vector[double] get_positions() except +
        vector[double] get_velocities() except +
        vector[double] get_initial_positions() except +
        vector[double] get_initial_velocities() except +
