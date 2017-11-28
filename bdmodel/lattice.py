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
        """Reset surface heights."""
        self._heights = np.zeros(self._length, dtype=int)
        return

    def evolve(self, nprtcls):
        """Evolve the system by depositing nprtcls particles.

        This method is to be overridden when writing child 
        classes with specific implementation of the models.
            
        Parameters
        ----------
            nprtcls : int
                Number of particles to deposit.

        """
        pass
        return 

    def deposit(self, j_latt):
        """Deposit a particle at lattice point j_latt.

        This method is to be overridden when writing child 
        classes with specific implementation of the models.
        
        """
        pass
        return

    def meanheight(self):
        """Calculate mean height of lattice."""
        return np.mean(self._heights)
    
    def width(self):
        """Calculate the interface width of the lattice."""
        mean_sqrheight = np.mean(self._heights*self._heights)
        return np.sqrt(mean_sqrheight - np.power(self.meanheight(), 2))

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
                system 1 MSC = self.length deposited particles).
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
            self._ts_MCS = t_ini*np.logspace(0, self._nmeasures - 1, 
                                                num=self._nmeasures, 
                                                base=spacefactor, dtype=int)
            # Remove possible repeated values
            self._ts_MCS = np.unique(self._ts_MCS)

            # Update the number of measures if it has changed
            if self._nmeasures != (self._ts_MCS).size:
                # TODO: Specify the Lattice instance name
                self._nmeasures = (self._ts_MCS).size
                print("The number of measures has been "
                                        "reduced to {}".format(self._nmeasures))
            
        else:
            # Check the intervals are big enought so that we do not get
            # repeated time values when we truncate to integers.
            if nmeasures > measuretime:            
                raise ValueError("The number of measures is bigger than"
                                                    "the number timesteps.")
            # Number of steps between measures in MCS.
            steps = measuretime/nmeasures
            self._ts_MCS = np.linspace(1, self._measuretime, 
                                          num=self._nmeasures, dtype=int)
        
        # Calculate the number of particles deposited before each measure
        self._ts_prtcl = self._length*self._ts_MCS

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
    def _plot_measures(self, ys, log=False, **plt_args):
        """Auxilar plot function.
            
        Plots the ys versus the time of each measure.
        """
        # TODO: improve this doc
        if (ys.size != self._nmeasures):
            raise ValueError("ys size is different from the number of"
                                                                "measures.")
        # Use makers only by default
        if all(not k in plt_args for k in ("linestyle", "marker")):
            plt_args["linestyle"] = ""
            plt_args["marker"] = "o"
            plt_args["markersize"] = np.clip(300./self._nmeasures, 2, 8)
        
        size = 4
        plt.figure(figsize=(1.62*size, size))

        if log==False:
            plt.plot(self._ts_MCS, ys, **plt_args)
        else:
            plt.loglog(self._ts_MCS, ys, **plt_args)

        plt.xlabel(r"$t (MCS)$")
        plt.tight_layout()
        
        return 

    def plot_meanheight(self, log=False, **plt_args):
        self._plot_measures(self._meanheights, log=log, **plt_args)
        plt.show()
        return

    def plot_width(self, log=True, **plt_args):
        self._plot_measures(self._widths, log=log, **plt_args)
        
        # Draw line showing expected slope
        # TODO: add offset to line and label with the value of beta
        def logline(t):
            return self._widths[0]*np.power(t/self._ts_MCS[0], self.beta)

        plt.plot(self._ts_MCS, logline(self._ts_MCS))
        plt.show()
        return

    # TODO: merge wiwth non averaged function
    def plot_meanheight_ravg(self, log=False, **plt_args):
        self._plot_measures(self._meanheights_ravg, log=log, **plt_args)
        plt.show()
        return
    
    # TODO: merge with non averaged function
    def plot_width_ravg(self, log=True, **plt_args):
        """Plot the run averaged interface width against time.

        Parameters
        ----------
            log : bool
                If true (default), loglog axis will be used.
            **plt_args : 
                Keyword arguments to be passed to 
                :func:`matplotlib.pyplot.plot` function.

        """ 
        self._plot_measures(self._widths_ravg, log=log, **plt_args)

        # Draw line showing expected slope
        # TODO: add offset to line and label with the value of beta
        def logline(t):
            return self._widths_ravg[0]*np.power(t/self._ts_MCS[0], self.beta)

        plt.plot(self._ts_MCS, logline(self._ts_MCS))
        plt.ylabel(r"$w$")
        plt.show()

        return


class BDLattice(Lattice):
    """Lattice with evolution following the ballistic deposition model.

    """

    # Kardar-Parisi-Zhang scaling exponent
    _z_KPZ = 3./2.
    _alpha_KPZ = 0.5
    _beta_KPZ = 1./3.
    
    def __init__(self, length, seed=None):
        Lattice.__init__(self, length, seed=seed)

        # Initialize the scaling exponents of the lattice to those
        # of the KPZ universality class.
        self.z = BDLattice._z_KPZ
        self.alpha = BDLattice._alpha_KPZ
        self.beta = BDLattice._beta_KPZ


    def evolve(self, nprtcls):
        """Evolve the system by depositing nprtcls particles.
    
        The depositions are made using the random deposition model.
            
        Parameters
        ----------
            nprtcls : int
                Number of particles to deposit.

        """
        self._heights = c_evolve.evolveBD(self._heights, nprtcls, self._pbc)
        return 

    def deposit(self, j_latt):
        """Deposit a particle at lattice point j_latt.

        Ballistic deposition model is used.
        
        """
        self._heights = c_evolve.depositBD(j_latt, self._heights, self._pbc)
        return
    
class RDLattice(Lattice):
    """Lattice with evolution following the random deposition model.

    """

    # Random deposition scaling exponents
    _beta_RD = 1./2.
    
    def __init__(self, length, seed=None):
        Lattice.__init__(self, length, seed=seed)

        # Initialize the scaling exponents of the lattice to those
        # of the KPZ universality class.
        self.beta = RDLattice._beta_RD


    def evolve(self, nprtcls):
        """Evolve the system by depositing nprtcls particles.
    
        The depositions are made using the random deposition model.
            
        Parameters
        ----------
            nprtcls : int
                Number of particles to deposit.

        """
        self._heights = c_evolve.evolveRD(self._heights, nprtcls, self._pbc)
        return 

    def deposit(self, j_latt):
        """Deposit a particle at lattice point j_latt.

        Ballistic deposition model is used.
        
        """
        self._heights = c_evolve.depositRD(j_latt, self._heights, self._pbc)
        return
     


def plot_scaledwidth(latt_list, alpha=None, z=None, log=True):
    """Scale different (interface width)-time plots to show data collapse.

    Parameters
    ----------
        latt_list : :class:`Lattice` list 
            List with simulation on different lattices.
        alpha : float
            Scaling exponent :math:`\\alpha`. If not provided uses
            the attribute :any:`Lattice.alpha` in the :class:`Lattice`
            instances.
        z : float
            Scaling exponent :math:`z`. If not provided uses
            the attribute :any:`Lattice.z` in the :class:`Lattice`
            instances.

    """
    size = 4
    plt.figure(figsize=(1.62*size, size))
    
    alpha_is_given = not (alpha == None)
    z_is_given = not (z == None)

    for latt in latt_list:
        # Scaling factors
        if not alpha_is_given:
            alpha = latt.alpha
        if not z_is_given:
            z = latt.z

        w_scalefactor = np.power(latt.length, - float(alpha))
        t_scalefactor = np.power(latt.length, - float(z))

        plt.loglog(t_scalefactor*latt._ts_MCS, w_scalefactor*latt._widths_ravg, 
                   label="L={}".format(latt.length))

    plt.xlabel(r"$t/L^z (MCS)$")
    plt.ylabel(r"$w/L^\alpha$")
    plt.legend()
    
    plt.show()

    return
