from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

import numpy as np
import random
from . import c_evolve

class Lattice:

    """Lattice object.
    
    Parameters
    ----------
        length : int
            Number of lattice points.

    Attributes
    ----------
        

    """

    def __init__(self, length, seed=None):
        self._length = length 

        # Initialize the random number generator
        if seed == None:
            seed = random.SystemRandom().randint(0, 32767)
        c_evolve.seed(seed)

        # Create vector for the implementation of the periodic
        # boundary conditions.
        self._pbc = np.zeros(length + 2, dtype=int)
        self._pbc[0] = length - 1;
        self._pbc[length + 1] = 0;
        for i in range(1, length + 1):
            self._pbc[i] = i - 1;
        
        # Initialize the height of the lattice points
        self._heights = np.zeros(length, dtype=int)

    @property
    def length(self):
        """Get length of lattice."""
        return self._length

    @property 
    def heights(self):
        """Get heights of the lattice points."""
        return self._heights

    def evolve(self, n):
        """Evolve the system throwing n particles.
            
        Parameters
        ----------
            n : int
                Number of particles to throw.

        """
        self._heights = c_evolve.evolve(self._heights, n, self._pbc)
        return 

    def deposit(self, j):
        """Deposit a particle at lattice point j.
        
        """
        self._heights = c_evolve.deposit(j, self._heights, self._pbc)
        return
        
        
        
