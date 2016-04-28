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

"""Module containing the `GSEResult` class."""

# TO-DO:
# * GSEResult should inherit from "mHGResult" (new class in XL-mHG package)?
# * E-score calculation should be part of XL-mHG package

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging

import numpy as np
from scipy.stats import hypergeom

from ..basic import GeneSet

logger = logging.getLogger(__name__)


class GSEResult(object):
    """Result of an XL-mHG-based test for gene set enrichment in a ranked list.

    Parameters
    ----------
    n_star: int
        The cutoff at which the XL-mHG test statistic was attained.
    stat: float
        The XL-mHG test statistic.
    pval: float
        The XL-mHG p-value.
    N: int
        The total number of genes in the analysis.
    X: int
        The XL-mHG X parameter.
    L: int
        The XL-mHG L parameter.
    gene_set `genometools.basic.GeneSet`
        The gene set.
    indices: np.ndarray of integers
        The indices of the gene set genes in the ranked list.
    genes: list or tuple of str
        The gene names in order of their appearance in the ranked list.
    """
    def __init__(self, n_star, stat, pval, N, X, L,
                 gene_set, indices, genes):

        assert isinstance(n_star, int)
        assert isinstance(stat, float)
        assert isinstance(pval, float)
        assert isinstance(N, int)
        assert isinstance(X, int)
        assert isinstance(L, int)
        assert isinstance(gene_set, GeneSet)
        assert isinstance(indices, np.ndarray)
        assert isinstance(genes, (tuple, list))
        for g in genes:
            assert isinstance(g, str)

        self.n_star = n_star
        self.stat = stat
        self.pval = pval
        self.N = N
        self.X = X  # XL-mHG "X" parameter
        self.L = L  # XL-mHG "L" parameter

        self.indices = np.int32(indices).copy()
        self.indices.flags.writeable = False  # makes it hashable

        self.gene_set = gene_set
        self.genes = tuple(genes)

        # strength of enrichment
        self.escore_pval_thresh = None
        self.escore = None

    def __repr__(self):
        return '<%s object (gene_set_id=%s; pval=%.1e; hash=%d)' \
                % (self.__class__.__name__, self.gene_set.id, self.pval,
                   hash(self))

    def __str__(self):
        return '<%s object (gene_set=%s; pval=%.1e)>' \
                % (self.__class__.__name__, str(self.gene_set), self.pval)

    def __hash__(self):
        data = [
            self.n_star,
            self.stat,
            self.pval,
            self.N,
            self.X,
            self.L,
            self.indices.data,
            self.gene_set,
            self.genes
        ]
        return hash(tuple(data))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __setstate__(self, d):
        self.__dict__ = d
        self.indices.flags.writeable = False

    @property
    def K(self):
        return self.indices.size

    @property
    def k(self):
        return int(np.sum(self.indices < self.n_star))

    def calculate_escore(self, pval_thresh):
        """ Calculate XL-mHG E-score.  """
        N = self.N
        K = self.K
        X = self.X
        L = self.L
        indices = self.indices
        
        if K == 0 or L == N or K < X:
            return 0

        # k_max = 0
        fe_max = 0.0
        k = 1
        # pval = 1.0

        while k <= K and indices[k-1] < L:
            if k >= X:
                n = indices[k-1] + 1
                if pval_thresh == 1.0 or \
                        hypergeom.sf(k-1, N, K, n) <= pval_thresh:
                    fe = k / (K * (n / float(N)))
                    if fe >= fe_max:
                        fe_max = fe
                        # k_max = k
            k += 1

        self.escore_pval_thresh = pval_thresh
        self.escore = fe_max

    def get_pretty_format(self, omit_param=True, max_name_length=0):
        # TO-DO: clean up, commenting
        gs_name = self.gene_set.name
        if max_name_length > 0 and len(gs_name) > max_name_length:
            assert max_name_length >= 3
            gs_name = gs_name[:(max_name_length-3)] + '...'
        gs_str = gs_name + ' (%d / %d @ %d)' % \
                (self.k, len(self.genes), self.n_star)
        param_str = ''
        if not omit_param:
            param_str = ' [X=%d,L=%d,N=%d]' % (self.X, self.L, self.N)
        escore_str = ''
        if self.escore is not None:
            escore_str = ', e=%.1fx' % self.escore
        details = ', p=%.1e%s%s' % (self.pval, escore_str, param_str)
        return '%s%s' % (gs_str, details)
        
    def get_pretty_GO_format(self, GO, omit_acc=False, omit_param=True,
                             max_name_length=0):
        # accepts a GOParser object ("GO")
        # TO-DO: clean up, commenting
        term = GO.terms[self.gene_set.id]
        term_name = term.get_pretty_format(omit_acc=omit_acc,
                                           max_name_length=max_name_length)
        term_str = term_name + ' (%d)' % (len(self.genes))
        param_str = ''
        if not omit_param:
            param_str = ' [X=%d,L=%d,N=%d]' % (self.X, self.L, self.N)
        escore_str = ''
        if self.escore is not None:
            escore_str = ', e=%.1fx' % self.escore
        details = ', p=%.1e%s%s' % (self.pval, escore_str, param_str)
        return '%s%s' % (term_str, details)
