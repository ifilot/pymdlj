# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from .pymdlj cimport Cell, Parameters
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
import numpy as np
import sys
import os

cdef class PyMDLJ:
    def __cinit__(self):
        pass

    def simulate(self, filename):
        cdef shared_ptr[Parameters] parameters
        cdef Cell* cell

        # check if filename exists, else throw error
        if not os.path.exists(filename):
            raise Exception('Cannot open file: ' + filename)

        # build parameter set
        parameters = shared_ptr[Parameters](new Parameters(filename.encode('utf-8')))

        # get some parameters
        stepsize = parameters.get().get_param[double]("stepsize".encode('utf-8'))
        nrsteps = parameters.get().get_param[uint]("nrsteps".encode('utf-8'))

        # build simulation object
        cell = new Cell(parameters)
        
        for i in range(0, nrsteps):
            cell.integrate(i, stepsize);