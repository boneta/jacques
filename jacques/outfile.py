"""
=======================================================================
  Distances/energies file manipulation (.out)
=======================================================================

  Classes
  -------

    OutFile

"""

import itertools

import numpy as np


class OutFile:
    """
        Class to manipulate distance/energies files (.out)

        ##   X  #  dist  Etot  Eqm  ndx  dist_ref

        Parameters
        ----------
        filename : str / list
            path of .out file(s) to read

        Attributes
        ----------
        dim : int
            dimensions, number of distances/constraints
        n : int
            number of points
        data : ndarray(n, dim+2)
            matrix of data
        dist : ndarray(n, dim)
            distances
        Etot : ndarray(n)
            total energies
        Eqm : ndarray(n)
            QM energies
        ndx : ndarray(n, dim)
            indexes of each distance/constraint
        dist_ref : ndarray(n, dim)
            distances of reference

        Properties
        ----------
        header : str
            header of .out file
    """

    def __init__(self, filename=None) -> None:
        self.dim = 0
        for i in ('_data', 'dist', 'Etot', 'Eqm', 'ndx', 'dist_ref'):
            setattr(self, i, np.array([]))
        if filename is not None:
            self.read(filename)

    def __bool__(self) -> bool:
        return bool(self.data.shape[0] > 0)

    @property
    def header(self) -> str:
        return f"##   {self.dim}  #  DIST  Etot  Eqm  INDX  DIST_REF"

    @property
    def data(self) -> 'np.ndarray':
        return self._data

    @data.setter
    def data(self, data:'np.ndarray') -> None:
        self._data = data
        # guess dimension of .out file based on number of columns
        dim = (self._data.shape[1] - 2.) / 3.
        if dim.is_integer():
            self.dim = int(dim)
        else:
            raise ValueError("Number of columns does not comply with convention")
        dim = self.dim
        # total number of points
        self.n = data.shape[0]
        # build references to fields
        self.dist = data[:, :dim+1]
        self.Etot = data[:, dim]
        self.Eqm = data[:, dim+1]
        self.ndx = data[:, dim+2:2*(dim+1)]
        self.dist_ref = data[:, 2*(dim+1):]

    def read(self, filename) -> None:
        """
            Read .out file

            Parameters
            ----------
            filename : str / list
                path of .out file(s) to read
        """
        if isinstance(filename, str):
            filename = [filename]
        out_data_all = []
        for out in filename:
            with open(out, 'r') as f:
                out_data_all.append(np.loadtxt(f, dtype=float, comments='#'))
        self.data = np.concatenate(out_data_all)

    def write(self, filename) -> None:
        """
            Write .out file

            Parameters
            ----------
            filename : str
                path of .out file to write
        """
        dim = self.dim
        fmt = r'%12.4f  '*dim + r'%20.10f  %20.10f  ' + r'%5d  '*dim + r'%12.4f  '*dim
        np.savetxt(filename, self.data, fmt=fmt, comments='#', header=self.header[1:])

    def zero(self, zero=0.) -> None:
        """
            Set minimum energy to zero

            Parameters
            ----------
            zero : float
                definition of zero to set the minimum
        """
        self.Etot -= self.Etot.min() - zero
        self.Eqm -= self.Eqm.min() - zero

    def unique(self, by='ndx') -> None:
        """
            Remove duplicates

            Parameters
            ----------
            by : str
                field to use to evaluate uniqueness (def: ndx)
        """
        by_data = getattr(self, by)
        self.data = self.data[np.unique(by_data, axis=0, return_index=True)[1]]

    def sort(self, by='ndx', reverse=False) -> None:
        """
            Sort by a field values in ascending order

            Parameters
            ----------
            by : str
                field to sort by (def: ndx)
            reverse : bool
                sort in descending order
        """
        by_data = getattr(self, by)
        if by_data.ndim == 1:
            self.data = self.data[by_data.argsort()]
        else:
            sort_keys = [by_data[:, i] for i in reversed(range(by_data.shape[1]))]
            self.data = self.data[np.lexsort(sort_keys, axis=0)]
        if reverse:
            self.data = self.data[::-1]

    def missing_ndx(self, max_ndx=None, min_ndx=0) -> 'np.ndarray':
        """
            Find missing indexes

            Parameters
            ----------
            max_ndx : int / list(dim)
                maximum value of index to look for
                if int, all dimensions are set to this value
                if list, each dimension is set to the corresponding value
                if None, the maximum found is used for each dimension
            min_ndx : int / list(dim)
                minimum value of index to look for (def: 0)
                if int, all dimensions are set to this value
                if list, each dimension is set to the corresponding value

            Returns
            -------
            ndarray of ndarray
                list of missing indexes
        """
        # check max and min extremes to look for
        if max_ndx is None:
            max_ndx = np.max(self.ndx, axis=0)
        elif isinstance(max_ndx, int):
            max_ndx = np.full(self.dim, max_ndx)
        if isinstance(min_ndx, int):
            min_ndx = np.full(self.dim, min_ndx)
        if len(max_ndx) != self.dim or len(min_ndx) != self.dim:
            raise ValueError("'max_ndx' and 'min_ndx' must have the same length as 'dim'")
        # sets of tuples of indexes
        ndx_all = set(itertools.product(*[np.arange(min_ndx[i], max_ndx[i]+1) for i in range(self.dim)]))
        ndx_data = set(tuple(i) for i in self.ndx.astype(int))
        ndx_missing = np.array(list(ndx_all - ndx_data), dtype=int)
        # sort
        if ndx_missing.size > 0:
            sort_keys = [ndx_missing[:, i] for i in reversed(range(ndx_missing.shape[1]))]
            ndx_missing = ndx_missing[np.lexsort(sort_keys, axis=0)]
        return ndx_missing
