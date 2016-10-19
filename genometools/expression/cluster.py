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

"""Methods for clustering expression data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import sys
import logging

import pandas as pd
import numpy as np

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from . import ExpMatrix

logger = logging.getLogger(__name__)

_linkage_warn = ('Function scipy.cluster.hierarchy.linkage with "method = '
                 '\'average\'" sometimes crashes. If that happens, call '
                 'cluster_genes() with a different method (e.g., using '
                 '"method = \'median\'").')

def bicluster(
        matrix,
        gene_cluster_metric='correlation',
        sample_cluster_metric='euclidean',
        cluster_method='average',
        reverse_genes=False,
        reverse_samples=False):

    matrix = cluster_genes(
        matrix,
        metric=gene_cluster_metric,
        method=cluster_method,
        reverse=reverse_genes
    )

    matrix = cluster_samples(
        matrix,
        metric=sample_cluster_metric,
        method=cluster_method,
        reverse=reverse_samples
    )

    return matrix


def _cluster_ndarray(X, metric, method, reverse):
    assert isinstance(X, np.ndarray)
    assert isinstance(metric, (str, _oldstr))
    assert isinstance(method, (str, _oldstr))
    assert isinstance(reverse, bool)

    distxy = squareform(pdist(X, metric=metric))
    R = dendrogram(linkage(distxy, method=method), no_plot=True)
    order_rows = np.int64([int(l) for l in R['ivl']])
    if reverse:
        order_rows = order_rows[::-1]
    return order_rows


def cluster_genes(matrix, metric='correlation', method='average',
                  reverse=False):
    assert isinstance(matrix, ExpMatrix)

    # note: method = 'average' sometimes causes kernel to crash
    logger.debug(_linkage_warn)

    order_rows = _cluster_ndarray(matrix.X, metric=metric, method=method,
                                  reverse=reverse)
    matrix = matrix.iloc[order_rows]
    return matrix


def cluster_samples(matrix, metric='euclidean', method='average',
                    reverse=False):
    assert isinstance(matrix, ExpMatrix)

    # note: method = 'average' sometimes causes kernel to crash
    logger.debug(_linkage_warn)

    # workaround for scipy bug when supplied with a row of NaNs
    # see: https://github.com/scipy/scipy/issues/5142
    valid_rows = matrix.notnull().all(axis=1)
    filtered = matrix.loc[valid_rows]

    order_cols = _cluster_ndarray(filtered.X.T, metric=metric, method=method,
                                  reverse=reverse)
    matrix = matrix.iloc[:, order_cols]
    return matrix
