"""
=======================================================================
  Interpolation Operations
=======================================================================

Functions
---------

    point_in
    mv_inside
    smooth_array
    grid2point
    grid2grid

"""

import sys

import numpy as np
from scipy import interpolate as sp_interpolate

# fortranized functions
try:
    from jacques.utils import interpolation_fortran
    _fortran_local = True
except:
    sys.stderr.write(
        "WARNING: Interpolation fortran subroutines could not be imported\n")
    _fortran_local = False

def point_in(coord, grid, thr, fortran=True):
    '''
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
    '''

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
    '''
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
    '''

    coord = np.asarray(coord)
    dist = grid[:] - coord
    if coord.size != 1:
        dist = np.linalg.norm(dist, axis=1)
    # closest coordinates or original
    return grid[np.argmin(dist)] if force or np.min(dist) >= thr else coord

def smooth_array(grid, neighbours=5, method='sma', fortran=True):
    '''
        Smooth the values of an array

        Parameters
        ----------
        grid : ndarray(n)
            array of values
        neighbours : int, optional
            number of close points to consider on each side (def: 5)
        method : str, optional
            algorithm for smoothing (def: sma)
                sma: simple moving average (central)
        fortran : bool, optional
            use faster function written in Fortran (def: True)

        Returns
        -------
        grid_smooth : ndarray(n)
            array of smoothed values
    '''

    n = grid.shape[0]
    grid_smooth = np.zeros_like(grid)

    # simple moving average (central)
    if method.lower() == 'sma':
        if fortran and _fortran_local:
            grid_smooth = interpolation_fortran.smooth_sma(grid, neighbours)
        else:
            # main averages along array
            for i in range(neighbours, n-neighbours):
                grid_smooth[i] = np.mean(grid[i-neighbours:i+neighbours+1])
            # tails
            for n in range(1, neighbours)[::-1]:
                grid_smooth[n] = np.mean(grid[:n*2+1])
                grid_smooth[-n-1] = np.mean(grid[-n*2:])
            # extremes
            grid_smooth[0], grid_smooth[-1] = grid[0], grid[-1]

    else:
        raise ValueError('Unkown smoothing array method')

    return grid_smooth

def grid2point(coord, grid, Z, method='legacy', fortran=True, **kwargs):
    '''
        Interpolate values from grid to a single point

        Parameters
        ----------
        coord : array(2)
            array of coordinates of the point
        grid : ndarray(n,2)
            array of grid coordinates
        Z : ndarray(n)
            array of grid values
        method : str, optional
            interpolation algorithm (def: legacy)
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
    '''

    # default kwargs options
    kwargs_def = {'factor' : 0.01,
                  'kind' : 'linear'}
    kwargs = {**kwargs_def, **kwargs}

    # value taken from the nearest grid point
    if method.lower() == 'nearest':
        # calculate distance of all grid points to the point
        dist = np.linalg.norm(grid[:]-coord, axis=1)
        # value corresponding to the point with minimum distance
        point = Z[np.argmin(dist)]

    # legacy method, probably by Jean-Pierre Moreau
    # http://jean-pierre.moreau.pagesperso-orange.fr/f_function.html
    elif method.lower() == 'legacy':
        if fortran and _fortran_local:
            point = interpolation_fortran.grid2point_legacy(coord, grid, Z, kwargs['factor'])
        else:
            gauss = 1./kwargs['factor']
            dat_x = abs((grid[:, 0]-coord[0])*gauss)
            dat_y = abs((grid[:, 1]-coord[1])*gauss)
            w = np.exp(-np.power([x * np.sqrt(1.+y*y/(x*x))
                                  if x > y else 0.0
                                  if y == 0.0 else y * np.sqrt(1.+x*x/(y*y))
                                  for x, y in zip(dat_x, dat_y)], 2))
            point = np.sum(Z*w) / np.sum(w)

    # scipy's interp2d interpolation
    elif method.lower() == 'scipy':
        # build interpolation function
        inter = sp_interpolate.interp2d(x = grid[:, 0],
                                        y = grid[:, 1],
                                        z = Z,
                                        kind = kwargs['kind'],
                                        bounds_error = False,
                                        fill_value = None)
        # call function
        point = inter(coord[0], coord[1], assume_sorted=False)

    else:
        raise NameError('Unkown interpolation method for grid2point')

    return point

def grid2grid(grid, Z, grid2, method='legacy', spline=False, fortran=True, **kwargs):
    '''
        Interpolate values from grid to grid

        Parameters
        ----------
        grid : ndarray(n,2)
            original array of grid coordinates
        Z : ndarray(n)
            original array of grid values
        grid2 : ndarray(m,2)
            final array of grid coordinates
        method : str, optional
            interpolation algorithm (def: legacy)
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
    '''

    # default kwargs options
    kwargs_def = {'factor' : 0.15,
                  'kind' : 'cubic',
                  'fill_value' : True,
                  'smooth' : 100000}
    kwargs = {**kwargs_def, **kwargs}

    # legacy method, probably by Jean-Pierre Moreau
    # http://jean-pierre.moreau.pagesperso-orange.fr/f_function.html
    if method.lower() == 'legacy':
        if fortran and _fortran_local:
            Zf = interpolation_fortran.grid2grid_legacy(grid, Z, grid2, kwargs['factor'])
        else:
            gauss = 1./kwargs['factor']
            Zf = np.zeros((grid2.shape[0]))
            for i in range(Zf.shape[0]):
                dat_x = abs((grid[:, 0]-grid2[i, 0])*gauss)
                dat_y = abs((grid[:, 1]-grid2[i, 1])*gauss)
                w = np.exp(-np.power([x * np.sqrt(1.+y*y/(x*x))
                                      if x > y else 0.0
                                      if y == 0.0 else y * np.sqrt(1.+x*x/(y*y))
                                      for x, y in zip(dat_x, dat_y)], 2))
                Zf[i] = np.sum(Z*w) / np.sum(w)

    # scipy's griddata interpolation
    elif method.lower() == 'scipy':
        Zf = sp_interpolate.griddata(points=grid, values=Z, xi=grid2, method=kwargs['kind'])
        # fill NaN values with nearest interpolation method
        if kwargs['fill_value']:
            Zf_near = sp_interpolate.griddata(points=grid, values=Z, xi=grid2, method='nearest')
            Zf = np.asarray([near if np.isnan(grid) else grid for grid, near in zip(Zf, Zf_near)])
        # discart NaN values
        else:
            i = 0
            while i < Zf.size:
                if np.isnan(Zf[i]):
                    grid2 = np.delete(grid2, i, 0)
                    Zf = np.delete(Zf, i, 0)
                else:
                    i += 1

    else:
        raise ValueError('Unkown interpolation method for grid2grid')

    # SPLINEs
    if spline:
        inter = sp_interpolate.SmoothBivariateSpline(grid2[:, 0],
                                                     grid2[:, 1],
                                                     Zf,
                                                     kx = 5,
                                                     ky = 5,
                                                     s = kwargs['smooth'])      # big number to work well, otherwise sh*t
        Zf = inter.ev(grid2[:, 0], grid2[:, 1])

    return grid2, Zf
