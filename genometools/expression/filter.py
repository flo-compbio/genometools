# Copyright (c) 2016 Florian Wagner
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

"""Methods for filtering expression data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging

import numpy as np

from . import ExpMatrix

logger = logging.getLogger(__name__)


def filter_variance(matrix, top):
    """Filter genes in an expression matrix by variance.

    Parameters
    ----------
    matrix: ExpMatrix
        The expression matrix.
    top: int
        The number of genes to retain.

    Returns
    -------
    ExpMatrix
        The filtered expression matrix.
    """
    assert isinstance(matrix, ExpMatrix)
    assert isinstance(top, (int, np.integer))

    if top >= matrix.p:
        logger.warning('Variance filter has no effect '
                       '("top" parameter is >= number of genes).')
        return matrix.copy()

    var = np.var(matrix.X, axis=1, ddof=1)
    total_var = np.sum(var)  # total sum of variance
    a = np.argsort(var)
    a = a[::-1]
    sel = np.zeros(matrix.p, dtype=np.bool_)
    sel[a[:top]] = True

    lost_p = matrix.p - top
    lost_var = total_var - np.sum(var[sel])
    logger.info('Selected the %d most variable genes '
                '(excluded %.1f%% of genes, representing %.1f%% '
                'of total variance).',
                top, 100 * (lost_p / float(matrix.p)),
                100 * (lost_var / total_var))

    matrix = matrix.loc[sel]
    return matrix


def filter_mean(matrix, top):
    """Filter genes in an expression matrix by mean expression.

    Parameters
    ----------
    matrix: ExpMatrix
        The expression matrix.
    top: int
        The number of genes to retain.

    Returns
    -------
    ExpMatrix
        The filtered expression matrix.
    """
    assert isinstance(matrix, ExpMatrix)
    assert isinstance(top, int)

    if top >= matrix.p:
        logger.warning('Gene expression filter with `top` parameter that is '
                       '>= the number of genes!')
        top = matrix.p

    a = np.argsort(np.mean(matrix.X, axis=1))
    a = a[::-1]

    sel = np.zeros(matrix.p, dtype=np.bool_)
    sel[a[:top]] = True

    matrix = matrix.loc[sel]
    return matrix


def filter_percentile(matrix, top, percentile=50):
    """Filter genes in an expression matrix by percentile expression.

    Parameters
    ----------
    matrix: ExpMatrix
        The expression matrix.
    top: int
        The number of genes to retain.
    percentile: int or float, optinonal
        The percentile to use  Defaults to the median (50th percentile).

    Returns
    -------
    ExpMatrix
        The filtered expression matrix.
    """
    assert isinstance(matrix, ExpMatrix)
    assert isinstance(top, int)
    assert isinstance(percentile, (int, float))

    if top >= matrix.p:
        logger.warning('Gene expression filter with `top` parameter that is '
                       ' >= the number of genes!')
        top = matrix.p

    a = np.argsort(np.percentile(matrix.X, percentile, axis=1))
    a = a[::-1]

    sel = np.zeros(matrix.p, dtype=np.bool_)
    sel[a[:top]] = True

    matrix = matrix.loc[sel]
    return matrix