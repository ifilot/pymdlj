# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

from .pymdlj cimport Cell, Parameters
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr,unique_ptr
import numpy as np
import sys
import os

cdef class PyMDLJ:
    def __cinit__(self):
        pass

    def simulate(self, filename):
        cdef shared_ptr[Parameters] parameters
        cdef unique_ptr[Cell] cell

        # check if filename exists, else throw error
        if not os.path.exists(filename):
            raise Exception('Cannot open file: ' + filename)

        # build parameter set
        parameters = shared_ptr[Parameters](new Parameters(filename.encode('utf-8')))

        # get some parameters
        stepsize = parameters.get().get_param[double]("stepsize".encode('utf-8'))
        nrsteps = parameters.get().get_param[uint]("nrsteps".encode('utf-8'))
        nrparticles = parameters.get().get_param[uint]("nr_particles".encode('utf-8'))

        # build simulation object
        cell = unique_ptr[Cell](new Cell(parameters))

        # store intermediate results
        ekin = []
        epot = []
        etot = []

        # perform integration
        for i in range(0, nrsteps):
            cell.get().integrate(i, stepsize)

            # store intermediate energies
            ekin.append(cell.get().get_ekin())
            epot.append(cell.get().get_epot())
            etot.append(cell.get().get_etot())

        results = {
            'ekin' : ekin,
            'epot' : epot,
            'etot' : etot,
            'nrsteps': nrsteps,
            'nrparticles': nrparticles,
            'positions': np.array(cell.get().get_positions()).reshape(nrparticles,3),
            'velocities': np.array(cell.get().get_velocities()).reshape(nrparticles,3),
            'positions_initial': np.array(cell.get().get_initial_positions()).reshape(nrparticles,3),
            'velocities_initial': np.array(cell.get().get_initial_velocities()).reshape(nrparticles,3)
        }

        return results
