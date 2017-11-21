from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

from matplotlib import pyplot as plt
import numpy as np
import random
from . import c_evolve

class Lattice:
    """Lattice object.
    
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

    def _measure_run(self, ts_prtcl):
        """Measure mean height and interface width over one run.
        
        Note: this resets the lattice heights.
    
        Parameters
        ----------
            ts_prtcl : int array
                Times when measures are taken in units of deposited
                particles.
        """
        self.resetheights() 
        
        nmeasures = ts_prtcl.size

        # Create vector to store the mean heights and interface width
        # TODO: add to docs
        self._meanheights = np.zeros(nmeasures)
        self._widths = np.zeros(nmeasures)
        
        last_t = 0
        for j_t, t in enumerate(ts_prtcl):
            # Calculate the number of particles to be deposited
            # before the next measurement    
            n_prtcls = t - last_t
            self.evolve(n_prtcls)

            # Measure
            self._meanheights[j_t] = self.meanheight()
            self._widths[j_t] = self.width()

            last_t = t

        return
        

    def measure(self, measuretime, nmeasures, nruns=1, log=True):
        """Measure interface width and mean height over time.

        Parameters
        ----------
            measuretime : int
                Duration of the run in Monte Carlo steps (MCS) (in our
                system 1 MSC = self.length particles thrown).
            nmeasures : int
                Number of measures.
            nruns : int
                Number of simulation runs to average over.
            log : bool
                If true, the measures are taken in logarithmically 
                spaced intervals. Else, they are taken in linear
                intervals.

        """
        # TODO: DOCS: where are the results stored?
        # TODO: document the arguments and write methods if required
        # Store the arguments
        self._nruns = nruns
        self._measuretime = measuretime
        self._nmeasures = nmeasures
        self._logt = log

        # Generate a linearly or logarithmically spaced time sequence
        # according to the value of log argument.
        if log:
            # A logartihmically spaced sequence of times {t_i} can be
            # written like t_i = factor^{i}*t_0.
            # This factor verifies factor = (t_nmeasures/t_0)^(1/nmeasures)
            # with t_nmeasures = measuretime.
            t_ini = 10
            # Calculate the spacing factor
            spacefactor = np.power(self._measuretime/t_ini, 1./self._nmeasures)
            self._t_MCS_vec = t_ini*np.logspace(0, self._nmeasures - 1, 
                                                num=self._nmeasures, 
                                                base=spacefactor, dtype=int)
            # Remove possible repeated values
            self._t_MCS_vec = np.unique(self._t_MCS_vec)

            # Update the number of measures
            self._nmeasures = self._t_MCS_vec.size
            
        else:
            # Check the intervals are big enought so that we do not get
            # repeated time values when we truncate to integers.
            if nmeasures > measuretime:            
                raise ValueError("The number of measures is bigger than"
                                                    "the number timesteps.")
            # Number of steps between measures in MCS.
            steps = measuretime/nmeasures
            self._t_MCS_vec = np.linspace(1, self._measuretime, 
                                          num=self._nmeasures, dtype=int)
        
        # Calculate the number of particles deposited before each measure
        self._ts_prtcl = self._length*self._t_MCS_vec

        # Create vectors where the accumulated sum over runs will
        # be saved.
        meanheights_sum = np.zeros(self._nmeasures)
        widths_sum = np.zeros(self._nmeasures)

        # Loop over runs
        for j_run in range(nruns):
            # Measure the run
            self._measure_run(self._ts_prtcl)

            # Add the values to the average
            meanheights_sum += self._meanheights 
            widths_sum += self._widths 

        # Store the average in the object
        self._meanheights_ravg = meanheights_sum/nruns
        self._widths_ravg = widths_sum/nruns
        
        return

    # Plot methods
    def _plot_measures(self, y_vec, log=False, **plt_args):
        """Auxilar plot function.
            
        Plots the y_vec versus the time of each measure.
        """
        # TODO: improve this doc
        if (y_vec.size != self._nmeasures):
            raise ValueError("y_vec size is different from the number of"
                                                                "measures.")
        if "linestyle" in plt_args:
            plt_args["linestyle"] = ""
        
        size = 4
        plt.figure(figsize=(1.62*size, size))

        if log==False:
            plt.plot(self._t_MCS_vec, y_vec, **plt_args)
        else:
            plt.loglog(self._t_MCS_vec, y_vec, **plt_args)

        plt.xlabel(r"$t (MCS)$")
        plt.tight_layout()
        
        return 

    def plot_mheight(self, log=False, **plt_args):
        self._plot_measures(self._meanheights, log=log, **plt_args)
        plt.show()
        return

    def plot_width(self, log=True, **plt_args):
        self._plot_measures(self._widths, log=log, **plt_args)
        plt.show()
        return

    def plot_mheight_ravg(self, log=False, **plt_args):
        self._plot_measures(self._meanheights_ravg, log=log, **plt_args)
        plt.show()
        return
    
    def plot_width_ravg(self, log=True, **plt_args):
        """Plot the run averaged interface width against time.

        Parameters
        ----------
            log : bool
                If true (default), loglog axis will be used.

        """ 
        self._plot_measures(self._widths_ravg, log=log, **plt_args)
        plt.ylabel(r"$w$")
        plt.show()

        return
