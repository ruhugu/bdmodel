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

# This tells Cython that there is a c_evolve function defined elsewhere
# whose header is in "c_evolve.h".
# cdef is used to define c functions.
# http://notes-on-cython.readthedocs.io/en/latest/function_declarations.html
# cdef extern especifies that the function is defined elsewhere. This 
# can be used without the from "file.hpp" part, but I don't understand how.
# http://cython-docs2.readthedocs.io/en/latest/src/userguide/external_C_code.html
cdef extern from "c_evolve.h":
    void c_evolve(long int* in_heights, long int* out_heights, int length, 
                  long int nsteps)
cdef extern from "dranxor2/dranxor2C.h":
    void dranini_(int*)

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
def evolve(cnp.ndarray[long int, ndim=1, mode="c"] in_heights not None,
           n):
    """Evolves in_heights lattice heights throwing nsteps particles.
    
    This function uses the nearest neighbour sticking rule.
    
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

    c_evolve(<long int*> cnp.PyArray_DATA(in_heights),    
             <long int*> cnp.PyArray_DATA(out_heights), 
             length, n)

    return out_heights

def seed(int iseed): 
    dranini_(&iseed)
    return    
    
