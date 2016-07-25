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

import logging

import numpy as np

from . import ExpMatrix

logger = logging.getLogger(__name__)


def quantile_normalize(matrix, inplace=False, target=None):
    """Quantile normalization, allowing for missing values (NaN).

    In case of nan values, this implementation will calculate evenly
    distributed quantiles and fill in the missing data with those values.
    Quantile normalization is then performed on the filled-in matrix,
    and the nan values are restored afterwards.

    Parameters
    ----------
    matrix: `ExpMatrix`
        The expression matrix (rows = genes, columns = samples).
    inplace: bool
        Whether or not to perform the operation in-place. [False]
    target: `numpy.ndarray`
        Target distribution to use. needs to be a vector whose first
        dimension matches that of the expression matrix. If ``None``,
        the target distribution is calculated based on the matrix
        itself. [None]

    Returns
    -------
    numpy.ndarray (ndim = 2)
        The normalized matrix.
    """
    assert isinstance(matrix, ExpMatrix)
    assert isinstance(inplace, bool)
    if target is not None:
        assert isinstance(target, np.ndarray) and \
               np.issubdtype(target.dtype, np.float)

    if not inplace:
        # make a copy of the original data
        matrix = matrix.copy()

    X = matrix.X
    _, n = X.shape
    nan = []
     # fill in missing values with evenly spaced quantiles
    for j in range(n):
        nan.append(np.nonzero(np.isnan(X[:, j]))[0])
        if nan[j].size > 0:
            q = np.arange(1, nan[j].size + 1, dtype=np.float64) / \
                (nan[j].size + 1.0)
            fill = np.nanpercentile(X[:, j], 100 * q)
            X[nan[j], j] = fill

    # generate sorting indices
    A = np.argsort(X, axis=0, kind='mergesort')  # mergesort is stable

    # reorder matrix
    for j in range(n):
        matrix.iloc[:, j] = matrix.X[A[:, j], j]

    # determine target distribution
    if target is None:
        # No target distribution is specified, calculate one based on the
        # expression matrix.
        target = np.mean(matrix.X, axis=1)
    else:
        # Use specified target distribution (after sorting).
        target = np.sort(target)

    # generate indices to reverse sorting
    A = np.argsort(A, axis=0, kind='mergesort')  # mergesort is stable

    # quantile-normalize
    for j in range(n):
        matrix.iloc[:, j] = target[A[:, j]]

    # set missing values to NaN again
    for j in range(n):
        if nan[j].size > 0:
            matrix.iloc[nan[j], j] = np.nan

    return matrix
