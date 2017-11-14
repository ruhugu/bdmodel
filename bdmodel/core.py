from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

from matplotlib import pyplot as plt
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
        self.resetheights()

    @property
    def length(self):
        """Get length of lattice."""
        return self._length

    @property 
    def heights(self):
        """Get heights of the lattice points."""
        return self._heights

    def resetheights(self):
        """Reset lattice heights."""
        self._heights = np.zeros(self._length, dtype=int)
        return

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

    def meanheight(self):
        """Calculate mean height of lattice."""
        return np.mean(self._heights)
    
    def width(self):
        """Calculate the interface width of the lattice."""
        # TODO: This could be optimized in order to avoid
        # calculating the mean twice.
        sumsqr = np.power((self._heights - self.meanheight()), 2)
        return np.sqrt(np.sum(sumsqr))/self._length

    def _measure_run(self, prtcl_measure, nmeasures):
        """Measure mean height and interface width over one run.
        
        Note: this resets the lattice heights.
    
        Parameters
        ----------
            prtcl_measure : int
                Number of particles deposited between measures.
            nmeasures : int
                Number of measures.
        """
        # Reset the lattice
        self.resetheights() 

        # Store measure time and number of measures
        # TODO: add to docs
        self._nsteps = nsteps
        self._measuretime = measuretime
        
        # Create vector to store the mean heights results
        # TODO: add to docs
        self._mh_vec = np.zeros(nsteps)
        # Vector storing the width results
        # TODO: add to docs
        self._w_vec = np.zeros(nsteps)
        
        for j_m in range(nsteps - 1):
            self._mh_vec[j_m] = self.meanheight()
            self._w_vec[j_m] = self.width()
            self.evolve(prtcl_measure)

        self._mh_vec[nsteps - 1] = self.meanheight()
        self._w_vec[nsteps - 1] = self.width()

        return
        

    def measure(self, measuretime, nmeasures, nruns=1):
        """Measure interface width and mean height over time.

        Parameters
        ----------
            measuretime : int
                Duration of the run in Monte Carlo steps (MCS) (in our
                system 1 MSC = self.length particles thrown).
            nmeasures : int
                Number of measures.
            nruns : number of simulation runs to average over.

        """
        # TODO: DOCS: where are the results stored?
        # Store the number of runs.
        self._nruns = nruns

        # Create vectors where the accumulated sum over runs will
        # be saved.
        mh_vec_sum = np.zeros(nmeasures)
        w_vec_sum = np.zeros(nmeasures)

        # Number of steps between measures in MCS.
        steps = measuretime/nmeasures
        
        # Number of particle deposited between measures
        prtcl_measure = int(steps*self.length)

        # Loop over runs
        for j_run in range(nruns):
            # Measure the run
            self._measure_run(prtcl_measure, nmeasures)

            # Add the values to the average
            mh_vec_sum += self._mh_vec 
            w_vec_sum += self._w_vec 

        # Store the average in the object
        self._mh_vec_ravg = mh_vec_sum/nruns
        self._w_vec_ravg = w_vec_sum/nruns
        
        return

    # Plot methods
    def _plot_measures(self, y_vec, log=False):
        """Auxilar plot function.
            
        Plots the y_vec versus the time of each measure.
        """
        # TODO: improve this doc
        if (y_vec.size != self._nsteps):
            raise ValueError("y_vec size is different from the number of"
                                                                "measures.")
        t_vec = self._measuretime*np.arange(self._nsteps)        
        
        if log==False:
            plt.plot(t_vec, y_vec)
        else:
            plt.loglog(t_vec, y_vec)
        plt.show()
        
        return 

    def plot_mheight(self, log=False):
        self._plot_measures(self._mh_vec, log=log)
        return

    def plot_width(self, log=True):
        self._plot_measures(self._w_vec, log=log)
        return

    def plot_mheight_ravg(self, log=False):
        self._plot_measures(self._mh_vec_ravg, log=log)
        return
    
    def plot_width_ravg(self, log=True):
        self._plot_measures(self._w_vec_ravg, log=log)
        return

    
