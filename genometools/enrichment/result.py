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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
import hashlib
from collections import Iterable

import numpy as np
from scipy.stats import hypergeom

from xlmhg import mHGResult
from ..basic import GeneSet

logger = logging.getLogger(__name__)


class GSEResult(mHGResult):
    """Result of an XL-mHG-based test for gene set enrichment.

    Parameters
    ----------
    gene_set `genometools.basic.GeneSet`
        The gene set.
    N: int
        The total number of genes in the ranked list.
    indices: np.ndarray of integers
        The indices of the gene set genes in the ranked list.
    ind_genes: list of str
        The names of the genes corresponding to the indices.
    X: int
        The XL-mHG X parameter.
    L: int
        The XL-mHG L parameter.
    stat: float
        The XL-mHG test statistic.
    cutoff: int
        The cutoff at which the XL-mHG test statistic was attained.
    pval: float
        The XL-mHG p-value.
    pval_thresh: float, optional
        The p-value threshold used in the analysis. [None]
    escore_pval_thresh: float, optional
        The hypergeometric p-value threshold used for calculating the E-score.
        If not specified, the XL-mHG p-value will be used, resulting in a
        conservative E-score. [None]
    escore_tol: float, optional
        The tolerance used for calculating the E-score. [None]
    """
    def __init__(self, gene_set, N, indices, ind_genes, X, L,
                 stat, cutoff, pval,
                 pval_thresh=None, escore_pval_thresh=None, escore_tol=None):

        # call parent constructor
        mHGResult.__init__(self, N, indices, X, L, stat, cutoff, pval,
                           pval_thresh=pval_thresh,
                           escore_pval_thresh=escore_pval_thresh,
                           escore_tol=escore_tol)

        # type checks
        assert isinstance(gene_set, GeneSet)
        assert isinstance(ind_genes, Iterable)

        if len(ind_genes) != indices.size:
            raise ValueError('The number of genes must match the number of '
                             'indices.')

        self.gene_set = gene_set
        self.ind_genes = list(ind_genes)

    def __repr__(self):
        return '<%s object (N=%d, gene_set_id="%s", hash="%s">' \
                % (self.__class__.__name__,
                   self.N, self.gene_set._id, self.hash)

    def __str__(self):
        return '<%s object (gene_set=%s, pval=%.1e)>' \
                % (self.__class__.__name__, str(self.gene_set), self.pval)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return self.hash == other.hash
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def hash(self):
        data_str = ';'.join(
            [super().hash] +
            [str(repr(var)) for var in [self.gene_set, self.ind_genes]])
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    @property
    def genes_above_cutoff(self):
        return self.ind_genes[:self.k]

    def get_pretty_format(self, omit_param=True, max_name_length=0):
        # TO-DO: clean up, commenting
        gs_name = self.gene_set._name
        if max_name_length > 0 and len(gs_name) > max_name_length:
            assert max_name_length >= 3
            gs_name = gs_name[:(max_name_length-3)] + '...'
        gs_str = gs_name + ' (%d / %d @ %d)' % \
                (self.k, len(self.ind_genes), self.cutoff)
        param_str = ''
        if not omit_param:
            param_str = ' [X=%d,L=%d,N=%d]' % (self.X, self.L, self.N)
        escore_str = ''
        if self.escore is not None:
            escore_str = ', e=%.1fx' % self.escore
        details = ', p=%.1e%s%s' % (self.pval, escore_str, param_str)
        return '%s%s' % (gs_str, details)

    def get_pretty_GO_format(self, GO, omit_acc=False, omit_param=True,
                             max_name_length=0): # pragma: no cover
        # accepts a GOParser object ("GO")
        # TO-DO: clean up, commenting
        term = GO.terms[self.gene_set._id]
        term_name = term.get_pretty_format(omit_acc=omit_acc,
                                           max_name_length=max_name_length)
        term_str = term_name + ' (%d)' % (len(self.ind_genes))
        param_str = ''
        if not omit_param:
            param_str = ' [X=%d,L=%d,N=%d]' % (self.X, self.L, self.N)
        escore_str = ''
        if self.escore is not None:
            escore_str = ', e=%.1fx' % self.escore
        details = ', p=%.1e%s%s' % (self.pval, escore_str, param_str)
        return '%s%s' % (term_str, details)
