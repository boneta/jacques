"""
=======================================================================
  Nudged Elastic Band algorithm
=======================================================================

Functions
---------

    neb

"""

import sys

import numpy as np
import umbrellaint

from ..vector                   import vproj, vnorm
from ..interpolation            import point_in, mv_inside, grid2point

try:
    from . import neb_fortran
    _fortranization_local = True
except:
    _fortranization_local = False


def neb(coord1, coord2, grid_d, grid, z, mep_guess=None,
        nknots=150, step=0.01, spring=1., rescale=(0.,1.),
        imethod='lowess', max_iter=1e4, fortran=True, **kwargs):
    '''
        Nudged elastic band algorithm to calculate a
        minimum energy path (MEP) between two points
        over an irregular grid

        Jónsson et al. Nudged elastic band method for finding minimum energy paths of transitions. 1998, 59(8), 385–404
        Trygubenko et al. J. Chem. Phys. 2004, 120 (5), 2082–2094


        Parameters
        ----------
        coord1 : array(2)
            array of coordinates of the initial point
        coord2 : array(2)
            array of coordinates of the final point
        grid_d : float
            distance between grid points
        grid : ndarray(n,2)
            array of grid coordinates
        z : ndarray(n)
            array of grid values
        mep_guess : ndarray(n,2), optional
            array of coordinates of structures to use as initial guess
            otherwise a straight line is used (def: None)
        nknots : int, optional
            number of knots in the NEB (def: 150)
            ignored if a mep_guess is provided
        step : float, optional
            distance to advance on every iteration (def: 0.01)
        spring : float, optional
            force proporcionality of the spring in the NEB (def: 1.)
        rescale : tuple, optional
            rescale forces with respect maximum-minimum difference (def: (0.,1.))
        imethod : {nearest, lowess, scipy}, optional
            interpolation method for the points (def: lowess w/ span 0.01)
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

    # TODO: better step size
    # TODO: better optimization. newton-rhapson?

    # default kwargs options
    kwargs_def = {'span': 0.01}
    kwargs = {**kwargs_def, **kwargs}

    if not _fortranization_local:
        sys.stdout.write("WARNING: NEB fortran subroutines could not be imported\n")

    # calculate gradient
    grid_topol, grid_topol_d = umbrellaint.igrid_topol(grid_d, grid)
    dZ = umbrellaint.igrid_grad(z, grid_topol, grid_topol_d)

    # initial guess of path
    coord1, coord2 = np.array(coord1), np.array(coord2)
    if mep_guess is not None:
        nknots = mep_guess.shape[0]
    else:
        # straight line
        mep_guess = np.column_stack((np.linspace(coord1[0],coord2[0],num=nknots,endpoint=True),
                                     np.linspace(coord1[1],coord2[1],num=nknots,endpoint=True)))

    # move all guess points inside grid
    mep_guess = np.array([mv_inside(i, grid, grid_d) for i in mep_guess])


    # NEB algorithm
    if fortran and _fortranization_local:
        mep_path = neb_fortran.neb(mep_guess, grid, z, dZ, step, max_iter, grid_d, spring, kwargs['span'], rescale, 1)
    else:
        # initialize variables
        mep_path = mep_guess
        mep_last = mep_path
        T_path = np.zeros((nknots,2))
        P_path = np.zeros((nknots,2))
        S_path = np.zeros((nknots,2))
        P_mod = np.zeros(nknots)
        S_mod = np.zeros(nknots)
        max_step = step

        for iter in range(int(max_iter)):
            # energy of the MEP points
            E_path = np.array([grid2point(i, grid, z, imethod, fortran=fortran, **kwargs) for i in mep_path])

            for i in range(1,nknots-1):
                # tangent vector by direct interpolation
                # N = vnorm(mep_path[i+1,:] - mep_path[i-1,:])

                # tangent vector by bisection
                T_path[i,:] = vnorm( vnorm(mep_path[i,:] - mep_path[i-1,:]) + vnorm(mep_path[i+1,:] - mep_path[i,:]) )

                # gradient
                grad = - np.array([grid2point(mep_path[i, :], grid, dZ[:, 0], imethod, fortran=fortran, **kwargs),
                                   grid2point(mep_path[i, :], grid, dZ[:, 1], imethod, fortran=fortran, **kwargs)])

                # perpendicular force
                P_path[i,:] = grad - vproj(grad, T_path[i,:])

                # spring force
                S_path[i, :] = spring * ((mep_path[i+1, :] - mep_path[i, :]) +
                                         (mep_path[i-1, :] - mep_path[i, :]))

                # modulus
                P_mod[i] = np.linalg.norm(P_path[i,:])
                S_mod[i] = np.linalg.norm(S_path[i,:])

            P_mod_zero = np.isclose(P_mod, 0)
            S_mod_zero = np.isclose(S_mod, 0)

            # normalize forces to the maximum energy
            E_min = np.min(E_path[1:-1])
            E_max = np.max(E_path[1:-1])
            E_diff = E_max - E_min
            P_path = np.array([P_path[i]/P_mod[i] * ((E_path[i] - E_min)/E_diff * (rescale[1]-rescale[0]) + rescale[0])
                               if not P_mod_zero[i] else [0,0] for i in range(nknots)])
            S_path = np.array([S_path[i]/S_mod[i] if not S_mod_zero[i] else [0, 0] for i in range(nknots)])

            # adapt step size
            maxpoint = np.argmax(E_path[1:-1])
            step = min(max_step, abs(P_mod[maxpoint]/E_max))

            # final force vector
            F_path = P_path/2. + S_path/2.

            # move MEP points
            mep_path = mep_path + F_path * step

            # check if any point is out of the grid
            mep_path = np.array([mep_path[i] if point_in(mep_path[i], grid, grid_d/np.sqrt(2.)) else mep_last[i] for i in range(nknots)])

            # calculate convergence
            displacement = np.sum(np.linalg.norm(mep_path - mep_last, axis=1))/nknots

            # save mep path to compare
            mep_last = np.copy(mep_path)

            # print(displacement, iter, P_mod[maxpoint], P_mod[np.argmin(E_path[1:-1])], step)

            # check convergence
            if displacement < 1.e-4 or P_mod[maxpoint] < 1.e-3: break

    # return results
    return mep_path

