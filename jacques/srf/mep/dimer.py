"""
=======================================================================
  Dimer algorithm
=======================================================================

Functions
---------

    dimer

"""

import sys

import numpy as np
from scipy import optimize as sp_opt     #TODO: avoid scipy dependency

from jacques.umbrellaintegrator import umbrellaint
from ..vector                   import vproj, vnorm
from ..interpolation            import point_in, grid2point


def dimer(coord1, coord2, grid_d, grid, z, step=0.01, dimer_d=0.1,
          imethod='lowess', max_iter=1e4, fortran=True, **kwargs):
    '''
        Classic dimer algorithm to calculate a
        minimum energy path (MEP) between two points
        over an irregular grid

        Henkelman et al. J. Chem. Phys. 1999, 111(15), 7010


        Parameters
        ----------
        coord1 : ndarray(2)
            array of coordinates of the initial point
        coord2 : ndarray(2)
            array of coordinates of the final point
        grid_d : float
            distance between grid points
        grid : ndarray(n,2)
            array of grid coordinates
        z : ndarray(n)
            array of grid values
        step : float, optional
            distance to advance on every iteration (def: 0.01)
        dimer_d : float, optional
            length of each branch of the dimer (def: 0.1)
        imethod : {nearest, lowess, scipy}, optional
            interpolation method for the points (def: lowess)
            additional options will be pass to 'grid2point'
        max_iter : int, optional
            maximum number of iterations (def: 1e4)
        fortran : bool, optional
            use function written in Fortran (def: True)

        Returns
        -------
        mep_path : ndarray(n,2)
            array of coordinates of structures along the obtained MEP
    '''

    # angles between which the dimer is allowed to rotate
    BOUNDS = np.radians([-45., 45.])


    ##  internal functions  -------------------------------------------

    def _normal_vector(theta):
        '''Normalized vector pointing towards the final point and rotated by theta'''
        # rot_matrix = np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])    # clockwise
        rot_matrix = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])    # counterclockwise
        N = vnorm(coord2 - rc)
        return rot_matrix.dot(N)

    def _dimer_build(theta):
        '''Two components dimer vector rotated by theta'''
        N = _normal_vector(theta)
        # dimer as a matrix [[x1,y1],[x-1,y-1]]
        return np.array([rc+dimer_d*N,rc-dimer_d*N]), N

    def _dimer_energy(theta):
        '''Energy evaluation of the dimer at its location and rotated by theta'''
        dimer, N = _dimer_build(theta)
        if np.linalg.norm(dimer[0]-coord2) > np.linalg.norm(rc-coord2):
            return np.inf
        else:
            grad1 = np.array([grid2point(dimer[0],grid,dZ[:,0],imethod,fortran,**kwargs), grid2point(dimer[0],grid,dZ[:,1],imethod,fortran,**kwargs)])
            grad2 = np.array([grid2point(dimer[1],grid,dZ[:,0],imethod,fortran,**kwargs), grid2point(dimer[1],grid,dZ[:,1],imethod,fortran,**kwargs)])
            grad_p1 = grad1 - vproj(grad1, N)
            grad_p2 = grad2 - vproj(grad1, N)
            return np.linalg.norm(grad_p1) + np.linalg.norm(grad_p2)

    def _opt_theta(bounds, solver='scipy'):
        '''Calculate the theta angle that minimizes the dimer energy'''
        if solver == 'scipy':
            mini_result = sp_opt.minimize_scalar(fun    = _dimer_energy,
                                                 bounds = bounds,
                                                 method = 'Bounded')
            return mini_result.x
        elif solver == 'manual':
            angles = np.linspace(bounds[0],bounds[1],num=1000,endpoint=True)
            values = map(_dimer_energy, angles)
            return angles[np.argmin(values)]


    ##  dimer algorithm  ----------------------------------------------

    # initalizate list of structures
    coord1, coord2 = np.array(coord1), np.array(coord2)
    mep_path = [coord1]

    # calculate gradient
    grid_topol, grid_topol_d = umbrellaint.igrid_topol(grid_d, grid)
    dZ = umbrellaint.igrid_grad(z, grid_topol, grid_topol_d)

    for iter in range(int(max_iter)):
        rc = np.array(mep_path[-1])
        # find best rotated dimer
        theta = _opt_theta(BOUNDS, solver='scipy')
        dimer, N = _dimer_build(theta)
        # invert dimer if needed to point to coord2
        if np.linalg.norm(dimer[0]-coord2) > np.linalg.norm(dimer[1]-coord2):
            theta = -((np.pi/2)-theta)
        # next point
        rc_new = rc + step*_normal_vector(theta)
        if not point_in(rc_new,grid,grid_d) or np.linalg.norm(rc_new-coord2) < step:
            break
        else:
            mep_path.append(rc_new)
    else:
        sys.stdout.write("WARNING: Dimer search not converged\n")

    mep_path = np.asarray(mep_path)

    # return results
    return mep_path
