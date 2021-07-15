"""
=======================================================================
  Interpolation Operations
=======================================================================

Functions
---------

    point_in
    mv_inside
    smooth_array
    array2point
    array2array
    grid2point
    grid2grid

"""

import sys

import numpy as np
from scipy import interpolate as sp_interpolate

# fortranized functions
try:
    from . import interpolation_fortran
    _fortran_local = True
except:
    sys.stderr.write(
        "WARNING: Interpolation fortran subroutines could not be imported\n")
    _fortran_local = False


def point_in(coord, grid, thr, fortran=True):
    """
        Evaluate if some coordinates are inside a grid (close to any grid point)

        Parameters
        ----------
        coord : ndarray(n)
            array of coordinates of the point
        grid : ndarray(m,n)
            array of grid coordinates
        thr : float
            maximum value to consider close
        fortran: bool, optional
            use faster function written in Fortran (def: True)

        Returns
        -------
        bool
            True if inside, False otherwise
    """

    coord = np.asarray(coord)
    if fortran and _fortran_local:
        return bool(interpolation_fortran.point_in(coord, grid, thr))
    else:
        # absolute distance to input coordinates for all points
        dist = grid[:] - coord
        # calculate norm of distance vectors if not 1D
        if coord.size != 1:
            dist = np.linalg.norm(dist, axis=1)
            # diff = np.sqrt((diff*diff).sum(axis=1))
        return np.min(dist) <= thr

def mv_inside(coord, grid, thr, force=False):
    """
        Move coordinates to the closest point inside a grid if are out

        Parameters
        ----------
        coord : ndarray(n)
            array of coordinates of the point
        grid : ndarray(m,n)
            array of grid coordinates
        thr : float
            maximum value to consider close
        force : bool, optional
            do not check if the point is already inside the grid (def: False)

        Returns
        -------
        ndarray(n)
            array of coordinates of the point inside the grid
    """

    coord = np.asarray(coord)
    dist = grid[:] - coord
    if coord.size != 1:
        dist = np.linalg.norm(dist, axis=1)
    # closest coordinates or original
    return grid[np.argmin(np.abs(dist))] if force or np.min(np.abs(dist)) >= thr else coord

def smooth_array(z, method='sma', fortran=True, **kwargs):
    """
        Smooth the values of an array

        Parameters
        ----------
        z : ndarray(n)
            array of values
        method : str, optional
            algorithm for smoothing (def: sma)
                sma : simple moving average (central)
                    neighbours : int, optional
                        number of points around to average (def: 5)
        fortran : bool, optional
            use faster function written in Fortran (def: True)

        Returns
        -------
        z_smooth : ndarray(n)
            array of smoothed values
    """

    # default kwargs options
    kwargs_def = {'neighbours': 5}
    kwargs = {**kwargs_def, **kwargs}

    n = z.shape[0]
    z_smooth = np.zeros_like(z)

    # simple moving average (central)
    if method.lower() == 'sma':
        if fortran and _fortran_local:
            z_smooth = interpolation_fortran.smooth_sma(z, kwargs['neighbours'])
        else:
            # main averages along array
            for i in range(kwargs['neighbours'], n-kwargs['neighbours']):
                z_smooth[i] = np.mean(z[i-kwargs['neighbours']:i+kwargs['neighbours']+1])
            # tails
            for n in range(1, kwargs['neighbours'])[::-1]:
                z_smooth[n] = np.mean(z[:n*2+1])
                z_smooth[-n-1] = np.mean(z[-n*2:])
            # extremes
            z_smooth[0], z_smooth[-1] = z[0], z[-1]

    else:
        raise ValueError('Unkown smoothing array method')

    return z_smooth

def array2point(coord, grid, z, imethod='gauss', **kwargs):
    """
        Interpolate values from an array to a single point

        Parameters
        ----------
        coord : float
            point to interpolate
        grid : ndarray(n)
            array of grid coordinates
        z : ndarray(n)
            array of corresponding values
        imethod : str, optional
            interpolation method (def: gauss)
                nearest : nearest grid point
                gauss : gaussian smoothing
                    gauss_factor : int, optional
                        gaussian smoothing factor, bigger smoother (def: 0.01)
        Returns
        -------
        point : float
            interpolated value on the point
    """

    # default kwargs options
    kwargs_def = {'gauss_factor': 0.01}
    kwargs = {**kwargs_def, **kwargs}

    # value taken from the nearest grid point
    if imethod.lower() == 'nearest':
        point = z[np.argmin(np.abs(grid-coord))]

    # gaussian smoothing (based on grids:regular)
    elif imethod.lower() == 'gauss':
        gauss_inv = 1./kwargs['gauss_factor']  # gaussian parameter (?)
        dat_x = np.abs((grid-coord)*gauss_inv)
        w = np.exp(-np.power(dat_x,2))
        point = np.sum(z*w) / np.sum(w)

    else:
        raise ValueError('Unkown smoothing 1D method')

    return point

def array2array(grid, z, grid2=None, imethod='gauss', **kwargs):
    """
        Interpolate values from an array base to another array

        Parameters
        ----------
        grid : ndarray(n)
            array of coordinates
        z : ndarray(n)
            array of corresponding values
        grid2 : ndarray(m), optional
            final array of coordinates
            if None, the original array is taken
        imethod : str, optional
            algorithm for smoothing (def: gauss)
                gauss : gaussian smoothing (array2point)
                    gauss_factor : int, optional
                        gaussian smoothing factor, bigger smoother (def: 0.4)
        Returns
        -------
        grid2 : ndarray(m)
            final array of coordinates
        z_smooth : ndarray(n)
            array of smoothed values
    """

    # default kwargs options
    kwargs_def = {'gauss_factor': 0.4}
    kwargs = {**kwargs_def, **kwargs}

    grid2 = grid if grid2 is None else grid2

    # interpolate all the points
    z_smooth = np.asarray([array2point(i, grid, z, imethod, **kwargs) for i in grid2])

    return grid2, z_smooth

def grid2point(coord, grid, z, imethod='legacy', fortran=True, **kwargs):
    """
        Interpolate values from grid to a single point

        Parameters
        ----------
        coord : array(2)
            array of coordinates of the point
        grid : ndarray(n,2)
            array of grid coordinates
        z : ndarray(n)
            array of grid values
        imethod : str, optional
            interpolation method (def: legacy)
                nearest : nearest grid point
                legacy : legacy method
                    factor : float, optional
                        gaussian factor (def: 0.01)
                scipy : scipy's interp2d interpolation
                    kind : {linear, cubic, quintic}, optional
                        kind of spline interpolation to use (def: linear)
        fortran : bool, optional
            use function writen in Fortran (def: True)

        Returns
        -------
        point : float
            interpolated value on the point
    """

    # default kwargs options
    kwargs_def = {'factor': 0.01,
                  'kind': 'linear'}
    kwargs = {**kwargs_def, **kwargs}

    # value taken from the nearest grid point
    if imethod.lower() == 'nearest':
        # calculate distance of all grid points to the point
        dist = np.linalg.norm(grid[:]-coord, axis=1)
        # value corresponding to the point with minimum distance
        point = z[np.argmin(dist)]

    # legacy method, probably by Jean-Pierre Moreau
    # http://jean-pierre.moreau.pagesperso-orange.fr/f_function.html
    elif imethod.lower() == 'legacy':
        if fortran and _fortran_local:
            point = interpolation_fortran.grid2point_legacy(coord, grid, z, kwargs['factor'])
        else:
            gauss = 1./kwargs['factor']
            dat_x = abs((grid[:, 0]-coord[0])*gauss)
            dat_y = abs((grid[:, 1]-coord[1])*gauss)
            w = np.exp(-np.power([x * np.sqrt(1.+y*y/(x*x))
                                  if x > y else 0.0
                                  if y == 0.0 else y * np.sqrt(1.+x*x/(y*y))
                                  for x, y in zip(dat_x, dat_y)], 2))
            point = np.sum(z*w) / np.sum(w)

    # scipy's interp2d interpolation
    elif imethod.lower() == 'scipy':
        # build interpolation function
        inter = sp_interpolate.interp2d(x = grid[:, 0],
                                        y = grid[:, 1],
                                        z = z,
                                        kind = kwargs['kind'],
                                        bounds_error = False,
                                        fill_value = None)
        # call function
        point = inter(coord[0], coord[1], assume_sorted=False)

    else:
        raise NameError('Unkown interpolation method for grid2point')

    return point

def grid2grid(grid, z, grid2=None, imethod='legacy', spline=False, fortran=True, **kwargs):
    """
        Interpolate values from grid to grid

        Parameters
        ----------
        grid : ndarray(n,2)
            original array of grid coordinates
        z : ndarray(n)
            original array of grid values
        grid2 : ndarray(m,2), optional
            final array of grid coordinates
            if None, the original grid is taken
        imethod : str, optional
            interpolation method (def: legacy)
                legacy : legacy method
                    factor : float, optional
                        gaussian factor (def: 0.15)
                scipy : scipy's griddata interpolation
                    kind : {linear, nearest, cubic}, optional
                        method of interpolation (def: linear)
                    fill_value : bool, optional
                        fill the NaN values not covered by standard
                        interpolation with nearest aproximation,
                        otherwise delete them from grid (def: True)
        spline : bool, optional
            apply spline to correction values,
            recomendedwith scipy interpolation (def: False)
            smooth : float, optional
                smoothing factor for spline, needs to be
                large to avoid trashy results (def: 100000)
        fortran : bool, optional
            use function writen in Fortran (def: True)

        Returns
        -------
        grid2 : ndarray(m,2)
            final array of grid coordinates
        Zf : ndarray(m)
            final array of interpolated grid values
    """

    # default kwargs options
    kwargs_def = {'factor' : 0.15,
                  'kind' : 'cubic',
                  'fill_value' : True,
                  'smooth' : 100000}
    kwargs = {**kwargs_def, **kwargs}

    grid2 = grid if grid2 is None else grid2

    # legacy method, probably by Jean-Pierre Moreau
    # http://jean-pierre.moreau.pagesperso-orange.fr/f_function.html
    if imethod.lower() == 'legacy':
        if fortran and _fortran_local:
            zf = interpolation_fortran.grid2grid_legacy(grid, z, grid2, kwargs['factor'])
        else:
            gauss = 1./kwargs['factor']
            zf = np.zeros((grid2.shape[0]))
            for i in range(zf.shape[0]):
                dat_x = abs((grid[:, 0]-grid2[i, 0])*gauss)
                dat_y = abs((grid[:, 1]-grid2[i, 1])*gauss)
                w = np.exp(-np.power([x * np.sqrt(1.+y*y/(x*x))
                                      if x > y else 0.0
                                      if y == 0.0 else y * np.sqrt(1.+x*x/(y*y))
                                      for x, y in zip(dat_x, dat_y)], 2))
                zf[i] = np.sum(z*w) / np.sum(w)

    # scipy's griddata interpolation
    elif imethod.lower() == 'scipy':
        zf = sp_interpolate.griddata(points=grid, values=z, xi=grid2, method=kwargs['kind'])
        # fill NaN values with nearest interpolation method
        if kwargs['fill_value']:
            Zf_near = sp_interpolate.griddata(points=grid, values=z, xi=grid2, method='nearest')
            zf = np.asarray([near if np.isnan(grid) else grid for grid, near in zip(zf, Zf_near)])
        # discart NaN values
        else:
            i = 0
            while i < zf.size:
                if np.isnan(zf[i]):
                    grid2 = np.delete(grid2, i, 0)
                    zf = np.delete(zf, i, 0)
                else:
                    i += 1

    else:
        raise ValueError('Unkown interpolation method for grid2grid')

    # SPLINEs
    if spline:
        inter = sp_interpolate.SmoothBivariateSpline(grid2[:, 0],
                                                     grid2[:, 1],
                                                     zf,
                                                     kx = 5,
                                                     ky = 5,
                                                     s = kwargs['smooth'])      # big number to work well, otherwise sh*t
        zf = inter.ev(grid2[:, 0], grid2[:, 1])

    return grid2, zf
