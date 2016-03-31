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
from builtins import *

import sys
import logging

import pandas as pd
import numpy as np

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from . import ExpMatrix

logger = logging.getLogger(__name__)

def bicluster(E, sample_cluster_metric = 'euclidean'):
    E, a_row = cluster_genes(E)
    E, a_col  = cluster_samples(E, metric = sample_cluster_metric)
    return E, a_row, a_col

def _cluster_ndarray(X, metric, method, reverse):
    assert isinstance(X, np.ndarray)
    assert isinstance(metric, str)
    assert isinstance(method, str)
    assert isinstance(reverse, bool)

    distxy = squareform(pdist(X, metric = metric))
    R = dendrogram(linkage(distxy, method = method), no_plot = True)
    order_rows = np.int64([int(l) for l in R['ivl']])
    if reverse:
        order_rows = order_rows[::-1]
    return order_rows

def cluster_genes(E, metric = 'correlation', method = 'average', reverse = False):
    assert isinstance(E, ExpMatrix)

    # note: method = 'average' sometimes causes kernel to crash
    logger.warning('Function scipy.cluster.hierarchy.linkage with "method = \'average\'" '
                   'sometimes crashes. If that happens, call cluster_genes() with a '
                   'different method (e.g., using "method = \'median\'").')

    order_rows = _cluster_ndarray(E.X, metric = metric, method = method, reverse = reverse)
    E = E.iloc[order_rows]
    return E, order_rows

def cluster_samples(E, metric = 'euclidean', method = 'average', reverse = False):
    assert isinstance(E, ExpMatrix)    

    # note: method = 'average' sometimes causes kernel to crash
    logger.warning('Function scipy.cluster.hierarchy.linkage with "method = \'average\'" '
                   'sometimes crashes. If that happens, call cluster_samples() with a '
                   'different method (e.g., using "method = \'median\'").')

    order_cols = _cluster_ndarray(E.X.T, metric = metric, method = method, reverse = reverse)
    E = E.iloc[:,order_cols]
    return E, order_cols
