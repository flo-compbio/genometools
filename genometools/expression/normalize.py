# Copyright (c) 2015, 2016 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Functions for normalizing expression data.

Note: Currently, only quantile normalization is implemented.
"""

import numpy as np

def _reorder_matrix(X, A):
    """Reorder matrix (column-by-column) according to an index matrix.

    For each column in ``X``, the column elements are re-ordered according to
    the indices contained in the corresponding column in ``A``.

    Parameters
    ----------
    X: 2-dimensional `numpy.ndarray`
        The matrix to be reordered.
    A: 2-dimensional `numpy.ndarray`
        The index matrix.

    """
    assert isinstance(X, np.ndarray)
    assert isinstance(A, np.ndarray)
    assert X.shape[1] == A.shape[1]

    n = X.shape[1]
    for j in range(n):
        X[:,j] = X[A[:,j],j]
    return X

def quantile_normalize(X):
    """Quantile-normalize a matrix.

    Performs quantile normalization as described in Bolstad et al.
    (PubMed ID: 12538238; DOI: 10.1093/bioinformatics/19.2.185).
    """

    assert isinstance(X, np.ndarray)

    # make a copy of X
    X = X.copy()

    # perform column-wise argsort on X
    A = np.argsort(X, axis = 0)

    # calculate inverse operation
    A_inv = np.argsort(A, axis = 0)
    
    # sort X column-wise
    X = _reorder(X, A)

    # calculate row-wise means
    mean = np.mean(X, axis = 1)
    n = X.shape[1]

    # set all values in each row to the row mean
    X = np.tile(mean, (n, 1)).T

    # apply inverse operation
    X = _reorder(X, A_inv)

    return X
