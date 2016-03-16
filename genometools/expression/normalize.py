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

logger = logging.getLogger(__name__)

def quantile_normalize(X, copy_matrix = True):
    """Quantile normalization, allowing for missing values (NaN).

    In case of nan values, this implementation will calculate evenly
    distributed quantiles and fill in the missing data with those values.
    Quantile normalization is then performed on the filled-in matrix,
    and the nan values are restored afterwards.

    Parameters
    ----------
    X: numpy.ndarray (ndim = 2)
        The expression matrix (rows = genes, columns = samples).
    copy_matrix: bool
        Whether or not to make a copy of the expression matrix. If set to
        False, the expression data will be modified in-place.
        
    Returns
    -------
    numpy.ndarray (ndim = 2)
        The normalized matrix.
    """

    assert isinstance(X, np.ndarray)
    assert isinstance(copy_matrix, bool)

    if copy_matrix:
        # make a copy of the original data
        X = X.copy()

    p, n = X.shape

    nan = []

     # fill in missing values with evenly spaced quantiles
    for j in range(n):
        nan.append(np.nonzero(np.isnan(X[:,j]))[0])
        if nan[j].size > 0:
            q = np.arange(1, nan[j].size + 1, dtype = np.float64) / (nan[j].size + 1.0)
            fill = np.nanpercentile(X[:,j], 100 * q)
            X[nan[j],j] = fill

    # generate sorting indices
    #A = np.argsort(X, axis = 0)
    A = np.argsort(X, axis = 0, kind = 'mergesort') # mergesort is stable

    # reorder matrix
    for j in range(n):
        X[:,j] = X[A[:,j],j]

    # calculate target distribution
    target = np.mean(X, axis = 1)
    
    # generate indices to reverse sorting
    #A = np.argsort(A, axis = 0)
    A = np.argsort(A, axis = 0, kind = 'mergesort') # mergesort is stable

    # quantile-normalize
    for j in range(n):
        X[:,j] = target[A[:,j]]

    # set missing values to NaN again
    for j in range(n):
        if nan[j].size > 0:
            X[nan[j],j] = np.nan

    return X

def quantile_normalize_exact(X, copy_matrix = True):
    """Quantile normalization, allowing for missing values.

    This implementation estimates exact quantiles for each column with missing
    values. This is VERY slow. DO NOT USE for matrices with more than a handful
    of rows.

    Parameters
    ----------
    X: numpy.ndarray (ndim = 2)
        The expression matrix (rows = genes, columns = samples).
    copy_matrix: bool
        Whether or not to make a copy of the expression matrix. If set to
        False, the expression data will be modified in-place.
        
    Returns
    -------
    numpy.ndarray (ndim = 2)
        The normalized matrix.
    """

    assert isinstance(X, np.ndarray)
    assert isinstance(copy_matrix, bool)

    logger.warning('The function quantile_normalize_exact() is painfully '
                   'slow! Use quantile_normalize() instead!')

    p,n = X.shape
    M = np.isnan(X) # indicator matrix for missing values
    num_missing = np.sum(M, axis = 0, dtype = np.int64)
    
    if copy_matrix:
        # make a copy of the original data
        X = X.copy()

    # generate sorting indices
    A = np.argsort(X, axis = 0)
    
    # calculate which percentiles (quantiles) to use
    # (we only need these for columns with missing values)
    perc = 100.0 * (np.arange(p, dtype = np.float64) / (p - 1.0))

    # calculate the quantiles for each column
    Q = np.empty((p, n), dtype = np.float64)
    for j in range(n):
        if num_missing[j] == 0:
            # no missing values, the percentiles correspond exactly to the sorted column vector as the quantiles
            Q[:,j] = X[A[:,j],j]
        else:
            # this column has missing values, we have to interpolate the quantiles
            Q[:,j] = np.nanpercentile(X[:,j], perc)

    target = np.mean(Q, axis = 1) # this is the target distribution
    
    for j in range(n):
        if num_missing[j] == 0:
            X[A[:,j],j] = target
        else:
            # - this column has missing values, we have to interpolate the quantiles
            # - we have to calculate the specific percentiles (quantiles) to use, given the number of missing values
            # - we then use interpolation to calculate those quantiles for the target distribution
            num_pres = p - num_missing[j]
            perc = 100.0 * (np.arange(num_pres, dtype = np.float64) / (num_pres - 1.0))
            ind = np.arange(num_pres, dtype = np.int64)
            X[A[:num_pres,j],j] = np.percentile(target, perc)
            
    return X
