"""
=======================================================================
  Vector Operations
=======================================================================

Functions
---------

    vproj
    vnorm

"""

import numpy as np


def vproj(vec1, vec2):
    '''
        Vector projection of vector 1 over 2

        Parameters
        ----------
        vec1 : ndarray(n)
            array to be projected
        vec2 : ndarray(n)
            array to project on

        Returns
        -------
        vec3 : ndarray(n)
            array projected
    '''
    
    # vector projection: u*(v·u)/|u|**2
    return vec2 * np.dot(vec1, vec2) / np.dot(vec2, vec2)

def vnorm(vec):
    '''
        Vector normalization

        Parameters
        ----------
        vec : ndarray(n)
            array to be normalized

        Returns
        -------
        vec_n : ndarray(n)
            array normalized
    '''
    
    # return vec / np.linalg.norm(vec)      # equivalent but slower
    return vec / np.sqrt(np.dot(vec, vec))
