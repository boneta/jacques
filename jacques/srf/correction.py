"""
=======================================================================
  Correction
=======================================================================

Functions
---------

    corr_1D
    corr_2D

"""

import numpy as np

from .interpolation import array2array, grid2grid


def corr_1D(grid, z, potential, correction=None, imethod='lowess', **kwargs):
    """
        Correction 1D

        Correct free energy values with differences
        between potential and high-level correction

        Parameters
        ----------
        grid : ndarray(n)
            array of grid coordinates
        z : ndarray(n)
            array of integrated free energy values
        potential : ndarray(m,2)
            array of coordinates and potential energies
        correction : ndarray(m), optional
            array of high-level corrected energies corresponding
            to the potential coordinates and energies
            if None, the values provided as potential energies
            are taken as differences directly
        imethod : {lowess}, optional
            interpolation method for the points (def: lowess w/ span=0.2)

        Returns
        -------
        z_corr : ndarray(n)
            array of corrected free energies
    """

    # default kwargs options
    kwargs_def = {'span': 0.2}
    kwargs = {**kwargs_def, **kwargs}

    # energy differences  ->  high_level - low_level
    diff = np.copy(potential)
    if correction is not None:
        diff[:,1] = correction[:] - potential[:,1]
    diff[:,1] -= np.min(diff[:,1])

    # interpolate smoothly grid coord to differences and sum to free energy
    z_corr = z + array2array(diff[:,0], diff[:,1], grid, imethod=imethod, **kwargs)
    z_corr -= np.min(z_corr)

    return z_corr

def corr_2D(grid, z, potential, correction=None, imethod='lowess', **kwargs):
    '''
        Correction 2D

        Correct free energy values with differences
        between potential and high-level correction

        Parameters
        ----------
        grid : ndarray(n,2)
            array of grid coordinates
        G : ndarray(n)
            array of integrated free energy values
        potential : ndarray(m,3)
            array of coordinates and potential energies
        correction : ndarray(m), optional
            array of high-level corrected energies corresponding
            to the potential coordinates and energies
            if None, the values provided as potential energies
            are taken as differences directly
        imethod : {lowess, scipy}, optional
            interpolation method for the points (def: lowess w/ span=0.2)

        Returns
        -------
        z_corr : ndarray(n)
            array of corrected free energies
    '''

    # default kwargs options
    kwargs_def = {'span': 0.2}
    kwargs = {**kwargs_def, **kwargs}

    # energy differences  ->  high_level - low_level
    diff = np.copy(potential)
    if correction is not None:
        diff[:,2] = correction[:] - potential[:,2]
    diff[:,2] -= np.min(diff[:,2])

    # interpolate smoothly grid coord to differences and sum to free energy
    z_corr = z + grid2grid(diff[:,0:2], diff[:,2], grid, imethod=imethod, fortran=True, **kwargs)[1]
    z_corr -= np.min(z_corr)

    return z_corr
