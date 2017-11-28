# -*- coding: utf-8 -*-

# Python extension through Cython following:
# http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id10
# Nice example:
# https://gist.github.com/phaustin/4973792

from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)


# Imports Cython declarations for numpy.
# cimport is used to import C data types, functions, etc defined in 
# other Cython file. Details: 
# http://cython.readthedocs.io/en/latest/src/userguide/sharing_declarations.html#the-cimport-statement
# In this case we are not importing Python numpy. We are loading the 
# Cython code that allows it to interact with numpy.
cimport numpy as cnp

# Enable Numpy-C-API access.
# Interesting:
# https://github.com/cython/cython/wiki/tutorials-numpy
# Numpy-C-API:
# https://docs.scipy.org/doc/numpy/reference/c-api.html
#np.import_array()
import numpy as np

# This tells Cython that there is a c_evolveBD function defined elsewhere
# whose header is in "c_evolveBD.h".
# cdef is used to define c functions.
# http://notes-on-cython.readthedocs.io/en/latest/function_declarations.html
# cdef extern especifies that the function is defined elsewhere. This 
# can be used without the from "file.hpp" part, but I don't understand how.
# http://cython-docs2.readthedocs.io/en/latest/src/userguide/external_C_code.html
cdef extern from "c_evolve.h":
    void c_evolveBD(long int* in_heights, long int* out_heights, int length, 
                    long int nsteps, long int* pbc)
    void c_evolveRD(long int* in_heights, long int* out_heights, int length, 
                    long int nsteps, long int* pbc)
    void c_depositBD(int j_latt, long int *heights, long int* pbc);
    void c_depositRD(int j_latt, long int *heights, long int* pbc);

cdef extern from "dranxor2/dranxor2C.h":
    void dranini_(int*)


# Define a new c type for a function taking the arguments of an
# of an evolution function. This is needed for passing functions
# as arguments in Cython functions.
# Do the same with deposit functions
ctypedef void (*ev_func)(long int*, long int*, int, long int, long int*)
ctypedef void (*dep_func)(int, long int*, long int*)


# Define a wrapper function that will act as a bridge between Python and 
# the C function <--- no se hasta que punto es esto totalmente cierto
# not None: by default, Cython allows arguments that meet the specified
# data type or that are None. In order to prevent the last behaviour, we 
# must add not None after the parameter name.
# http://docs.cython.org/en/latest/src/userguide/extension_types.html#extension-types-and-none

# I can't find where the sintax np.ndarray[...] is explained. However, from 
# this example we can see that the first argument is the datatype, ndim 
# refers to the dimension of the array and I think mode determines how 
# the array is stored in the memory. In this case we would use the c-way
# (whatever that is). What I have found about the matter:
# Here they say it is realated to "efficient buffer access". Other modes are
# also presented.
# https://github.com/cython/cython/wiki/enhancements-buffer
# The <...> before argument are type casts:
# http://cython.readthedocs.io/en/latest/src/reference/language_basics.html#type-casting
# If we use a Python variable as an argument of a Cython function with a
# specified type, automatic conversion will be attempted.
# http://cython.readthedocs.io/en/latest/src/userguide/language_basics.html

cdef evolve_wrapper(cnp.ndarray[long int, ndim=1, mode="c"] in_heights,
                    int n,
                    cnp.ndarray[long int, ndim=1, mode="c"] pbc,
                    ev_func evolvefunc):
    """Evolves in_heights lattice heights throwing nsteps particles.
    
    The evolution algorithm used is defined in evolvefunc.
    
    Parameters
    ----------
        in_heights : 1-d int array
            Array with the initial height of each lattice point.
        n : int
            Number of particles to be thrown over the lattice.

    Returns:
        out_heights : 1-d int array
            Array with the final height of each lattice point after
            the deposition of the n particles.

    """
    length = in_heights.size
    cdef cnp.ndarray[long int, ndim=1, mode="c"] out_heights = \
                                                    np.zeros(length, dtype=int)

    evolvefunc(
               <long int*> cnp.PyArray_DATA(in_heights),    
               <long int*> cnp.PyArray_DATA(out_heights), 
               length, n,
               <long int*> cnp.PyArray_DATA(pbc))

    return out_heights


cdef deposit_wrapper(int j_latt,
                     cnp.ndarray[long int, ndim=1, mode="c"] in_heights,
                     cnp.ndarray[long int, ndim=1, mode="c"] pbc,
                     dep_func depositfunc):
    """Deposit a particle at lattice point 'j_lat'.
    
    Takes the surface given by in_heights and returns the resulting 
    lattice after depositing the particle in the given point.

    The sticking algorithm used is defined in "depositfunc".
    
    Parameters
    ----------
        j_lat : int
            Lattice point where the particle will be depositated.
        in_heights : 1-d int array
            Array with the initial height of each lattice point.

    Returns:
        out_heights : 1-d int array
            Array with the final height of each lattice point after
            the deposition.

    """
    length = in_heights.size

    if j_latt >= length or j_latt < 0:
        raise IndexError("j_latt is out of range")
    
    cdef cnp.ndarray[long int, ndim=1, mode="c"] out_heights = \
                                                    np.zeros(length, dtype=int)

    out_heights = np.copy(in_heights, order="c")
    
    depositfunc(<int> j_latt,
                <long int*> cnp.PyArray_DATA(out_heights),
                <long int*> cnp.PyArray_DATA(pbc))
    
    return out_heights


def evolveBD(cnp.ndarray[long int, ndim=1, mode="c"] in_heights not None,
              int n, cnp.ndarray[long int, ndim=1, mode="c"] pbc not None):
    """Evolve the heights according to the ballistic deposition model.

    """
    out_heights = evolve_wrapper(in_heights, n, pbc, c_evolveBD)
    return out_heights


def evolveRD(cnp.ndarray[long int, ndim=1, mode="c"] in_heights not None,
              int n, cnp.ndarray[long int, ndim=1, mode="c"] pbc not None):
    """Evolve the heights according to the random deposition model.

    """
    out_heights = evolve_wrapper(in_heights, n, pbc, c_evolveRD)
    return out_heights


def depositBD(int j_latt,
              cnp.ndarray[long int, ndim=1, mode="c"] in_heights not None,
              cnp.ndarray[long int, ndim=1, mode="c"] pbc not None):
    """Deposit a particle using the ballistic deposition model.

    """
    out_heights = deposit_wrapper(j_latt, in_heights, pbc, c_depositBD)
    return out_heights


def depositRD(int j_latt,
              cnp.ndarray[long int, ndim=1, mode="c"] in_heights not None,
              cnp.ndarray[long int, ndim=1, mode="c"] pbc not None):
    """Deposit a particle using the random deposition model.
 
    """
    out_heights = deposit_wrapper(j_latt, in_heights, pbc, c_depositRD)
    return out_heights
    

def seed(int iseed): 
    """Initialize the random number generator in the c extension.
    
    Parameters
    ----------
        iseed : int
            Seed for the random number generator.

    """
    dranini_(&iseed)
    return    
