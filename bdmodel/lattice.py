from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

from matplotlib import pyplot as plt
import numpy as np
import random
from . import c_evolve

class Lattice(object):
    """Lattice object.
    
    Attributes
    ----------
        
    """

    List of attributes to be saved/loaded from file
    self._attrList = [  "length",
                        "heights",
                        "_nruns ",
                        "_measuretime",
                        "_nmeasures",
                        "_logt",
                        "_ts_MCS",
                        "_ts_prtcl"
                        "_meanheights",
                        "_widths",
                        "_meanheights_ravg",
                        "_widths_ravg",
                     ]
    
    def __init__(self, length, seed=None, heighttype=int):
        self._length = length 
        self._heighttype = heighttype

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
        self._heights = np.zeros(self._length, dtype=self._heighttype)
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

    # File I/O methods

    def savetxt(self, fname):
        with open(fname, 'w') as f:
            classname = type(self).__name__
            f.write('{}\n'.format(classname)) 
        
            for attr in self.attrlist:
                if type(attr) is np.ndarray:
                    s_attr = attr.tostring()
                else: 
                    s_attr = str(attr)
                f.write('{}\n'.format(s_attr))

    def loadtxt(self, fname):
        with open(fname, 'r') as f:
            line = f.readline()
            classname = line.strip()

            for
            
            classname = type(self).__name__
            f.write('{}\n'.format(classname)) 
        
            for attr in self.attrlist:
                if type(attr) is np.ndarray:
                    s_attr = attr.tostring()
                else: 
                    s_attr = str(attr)
                f.write('{}\n'.format(s_attr))


    # Plot methods
    def _plot_measures(self, ys, relinterval=[0,1], log=False,
                       reg=True, **plt_args):
        # TODO update docs
        """Auxilar plot function.
            
        Plots the results of the measures "ys" versus the time 
        (in :term:`MCS`) of each measure (self._ts_MCS).

        Parameters
        ----------
            ys : numpy array 
                Measure values to be plotted in the y-axis. 
            log : bool
                If True, plots in loglog scale. Defaults to False.
            **plt_args : keyword argument dict
                Dictionary with the keyword arguments to be passed
                to the :py:func:~matplotlib.pyplot.plot() function.
        """
        # TODO: improve this doc
        if (ys.size != self._nmeasures):
            raise ValueError("ys size is different from the number of"
                                                            "measures.")
        # Use only markers by default
        if all(not k in plt_args for k in ("linestyle", "marker")):
            plt_args["linestyle"] = ""
            plt_args["marker"] = "o"
            plt_args["markersize"] = np.clip(300./self._nmeasures, 2, 8)

        # Use the data in the given interval
        relinterval = np.array(relinterval)
        [ini, end] = (self._nmeasures*relinterval).astype(int)
        end -= 1
        
        # Define the x vector
        xs = self._ts_MCS[ini:end]
        ys = ys[ini:end]
        
        # Calculate the fitting line according to the plot type
        if reg:
            if log:
                reg_coef = np.polyfit(np.log10(xs), np.log10(ys), 1)
                def reg_func(t): return (np.power(10, reg_coef[1])
                                                        *np.power(t, reg_coef[0]))
            else:
                reg_coef = np.polyfit(xs, ys, 1)
                reg_func = np.poly1d(reg_coef) 
        
        size = 4
        plt.figure(figsize=(1.62*size, size))

        if log:
            plt.loglog(xs, ys, **plt_args)
        else:
            plt.plot(xs, ys, **plt_args)
        
        # Plot regression line
        if reg:
            plt.plot(xs, reg_func(xs)) 
            
            text_relx = 0.8
            text_rely = 0.05
            #relxpos = 1./4
            #relyoffset = 0.1
            #if log:
            #    text_x = np.power(xs[0], 1. - relxpos)*np.power(xs[-1], relxpos)
            #    text_y = reg_func(text_x)*np.power(ys[-1]/ys[0], relyoffset)
            #else:
            #    text_x = (1 - relposx)*xs[0] - relposx*xs[-1]
            #    text_y = reg_func(text_x) + relyoffset*(ys[-1] - ys[0])
            
            plt.text(text_relx, text_rely, "m = {:.5f}".format(reg_coef[0]), 
                     transform=(plt.gca()).transAxes)
            
        plt.xlabel(r"$t\,(MCS)$")
        plt.tight_layout()
        
        if reg:
            return reg_coef
        else:
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
    def plot_meanheight_ravg(self, log=False, reg=True, **plt_args):
        self._plot_measures(self._meanheights_ravg, log=log, reg=True,
                            **plt_args)
        plt.show()
        return
    
    # TODO: merge with non averaged function
    def plot_width_ravg(self, relinterval=[0,1], log=True, reg=False, 
                        **plt_args):
        # TODO update docs
        """Plot the run averaged interface width against time.

        Parameters
        ----------
            log : bool
                If true (default), loglog axis will be used.
            **plt_args : 
                Keyword arguments to be passed to 
                :func:`matplotlib.pyplot.plot` function.

        """ 
        reg_coef = self._plot_measures(self._widths_ravg, relinterval=relinterval,
                                       log=log, reg=reg, **plt_args)
        if reg:
            self.reg_beta = reg_coef[0]

        # Draw line showing expected slope
        # TODO: add offset to line and label with the value of beta
        #def logline(t):
        #    return self._widths_ravg[0]*np.power(t/self._ts_MCS[0], self.beta)

        #plt.plot(self._ts_MCS, logline(self._ts_MCS))
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

    def plot_width_ravg(self, log=True, reg=True, **plt_args):
        Lattice.plot_width_ravg(self, log=log, reg=reg, **plt_args)
        return 
     

class RDdiffLattice(Lattice):
    """Lattice with evolution following the random deposition model.

    """

    # Random deposition scaling exponents
    _beta_RD = 1./2.
    
    def __init__(self, length, ht=0.1, seed=None):
        # ht viene dado en unidades de MCS
        Lattice.__init__(self, length, seed=seed, heighttype=float)

        self._ht = ht

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
        # Translate prtcls into MCS and then to timesteps
        nsteps = int((float(nprtcls)/self.length)/self._ht)
        
        self._heights = c_evolve.evolveRDdiff(self._heights, self._ht, 
                                              nsteps, self._pbc)
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
