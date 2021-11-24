"""
=======================================================================
  Minimization
=======================================================================

Functions
---------

    minimize

"""

import sys

import numpy as np
import umbrellaint

from .interpolation import point_in, mv_inside, grid2point


def minimize(coord, grid_d, grid, z, opt='hc', step=0.01, imethod='lowess', max_iter=1e4, **kwargs):
    '''
        Analytical localization of minimum in a surface
        from some starting coordinates

        Parameters
        ----------
        coord : ndarray(2)
            array of coordinates of the point
        grid_d : float
            distance between grid points
        grid : ndarray(n,2)
            array of grid coordinates
        z : ndarray(n)
            array of grid values
        opt : str, optional
            algorithm to minimize (def: hc)
                sd: steepest descent
                hc: simple hill climbing
        step : float, optional
            distance to advance on every iteration (def: 0.01)
        imethod : {nearest, lowess, scipy}, optional
            interpolation method for the points (def: lowess)
            additional options will be pass to 'grid2point'
        max_iter : int, optional
            maximum number of iterations

        Returns
        -------
        min_path : ndarray(m,2)
            array of coordinates of structures during the minimization
    '''

    max_iter = int(max_iter)
    grid_d_red = grid_d/np.sqrt(2.)
    # grid_d = round(min({abs(grid[i0,0] - grid[i1,0]) for i0 in range(grid.shape[0]) for i1 in range(grid.shape[0])} - {0.}), 6)

    # move point inside grid
    coord = mv_inside(coord, grid, grid_d_red)

    # initalizate list of structures
    min_path = [coord]

    # steepest descent algorithm
    if opt.lower() == 'sd':
        # calculate gradient
        grid_topol, grid_topol_d = umbrellaint.igrid_topol(grid_d, grid)
        dZ = umbrellaint.igrid_grad(z, grid_topol, grid_topol_d)
        # minimize
        for iteration in range(1, max_iter+1):
            rc = min_path[-1]
            grad = [grid2point(rc, grid, dZ[:, i], imethod, **kwargs) for i in range(2)]
            grad_norm = np.linalg.norm(grad)
            rc_new = [rc[i] - grad[i]/grad_norm * step for i in range(2)]
            # break if out of grid or small grad
            if not point_in(rc_new, grid, grid_d_red) or grad_norm < 1e-3:
                break
            else:
                min_path.append(rc_new)
        if iteration >= max_iter:
            sys.stderr.write("WARNING: Steepest Descent minimization not converged\n")

    # simple hill climbing algorithm
    elif opt.lower() == 'hc':
        diag = step * np.sqrt(2.)
        for iteration in range(1, max_iter+1):
            rc = min_path[-1]
            # build sourrouding points to check
            around = [[rc[0]+step, rc[1]],
                      [rc[0]-step, rc[1]],
                      [rc[0],      rc[1]+step],
                      [rc[0],      rc[1]-step],
                      [rc[0]+diag, rc[1]+diag],
                      [rc[0]-diag, rc[1]+diag],
                      [rc[0]+diag, rc[1]-diag],
                      [rc[0]-diag, rc[1]-diag]]
            around_energy = np.array([grid2point(p, grid, z, imethod, **kwargs)
                                      if point_in(p, grid, grid_d_red) else np.inf for p in around])
            rc_new = around[np.argmin(around_energy)]
            min_path.append(rc_new)
            # check similarity of last structures
            if (iteration > 11 and np.std(np.array(min_path)[-10:, 0]) < step and np.std(np.array(min_path)[-10:, 1]) < step):
                break
        if iteration >= max_iter:
            sys.stderr.write("WARNING: Hill Climbing minimization not converged\n")

    else:
        raise NameError('Unkown minimization method')

    # return results
    return np.asarray(min_path)
