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

"""Module containing the `GSEAnalysis` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
from math import ceil

import numpy as np

from xlmhg import xlmhg_test

# from ..basic import GeneSet, GeneSetDB
from ..basic import GeneSetDB
from ..expression import ExpGenome
from . import GSEResult

logger = logging.getLogger(__name__)


class GSEAnalysis(object):
    """Test ranked gene lists for gene set enrichment using the XL-mHG test.

    Parameters
    ----------
    genome: `ExpGenome` object
        See :attr:`genome` attribute.
    gene_set_db: `GeneSetDB` object
        See :attr:`geneset_db` attribute.

    Attributes
    ----------
    genome: `ExpGenome` object
        The universe of genes.
    gene_set_db: `GeneSetDB` object
        The gene set database.

    Notes
    -----
    The class is initialized with a set of valid gene names (an `ExpGenome`
    object), as well as a set of gene sets (a `GeneSetDB` object). During
    initialization, a binary "gene-by-gene set" matrix is constructed,
    which stores information about which gene is contained in each gene set.
    This matrix is quite sparse, and requires a significant amount of memory.
    As an example, for a set of p = 10,000 genes and n = 10,000 gene sets, this
    matrix is of size 100 MB in the memory (i.e., p x n bytes).

    Once the class has been initialized, the function `get_enriched_gene_sets`
    can be called with a ranked list of genes, a significance threshold, and a
    set of parameters for the XL-mHG test. This function returns a list of
    `GSEResult` objects, one for each gene set that was found to be
    significantly enriched.
    """
    def __init__(self, genome, gene_set_db):

        assert isinstance(genome, ExpGenome)
        assert isinstance(gene_set_db, GeneSetDB)

        self.genome = genome
        self.gene_set_db = gene_set_db

        # generate annotation matrix by going over all gene sets
        logger.info('Generating gene-by-gene set membership matrix...')

        self._A = np.zeros((len(genome), gene_set_db.n), dtype=np.uint8)
        for j, gs in enumerate(self.gene_set_db.gene_sets):
            for g in gs.genes:
                try:
                    idx = self.genome.index(g)
                except ValueError:
                    pass
                else:
                    self._A[idx, j] = 1

    def __repr__(self):
        return '<%s object (genome=%s; gene_set_db=%s)>' \
               % (self.__class__.__name__,
                  repr(self.genome), repr(self.gene_set_db))

    def __str__(self):
        return '<%s object (%d genes in genome; %d gene sets)>' \
               % (self.__class__.__name__,
                  len(self.genome), len(self.gene_set_db))

    @property
    def genes(self):
        return self.genome.genes

    def get_enriched_gene_sets(
            self, ranked_genes, pval_thresh, X_frac, X_min, L,
            escore_pval_thresh=None, gene_set_ids=None, mat=None):
        """Tests gene set enrichment given a ranked list of genes.

        This function also calculates XL-mHG E-scores for the enriched gene
        sets, using ``escore_pval_thresh`` as the p-value threshold "psi".
        """
        if isinstance(X_frac, (int, np.integer)):
            X_frac = float(X_frac)

        # checks
        assert isinstance(ranked_genes, (list, tuple))
        for g in ranked_genes:
            assert isinstance(g, str)
        assert isinstance(pval_thresh, float)
        assert isinstance(X_frac, float)
        assert isinstance(X_min, int)
        assert isinstance(L, int)

        if escore_pval_thresh is not None:
            assert isinstance(escore_pval_thresh, float)
        if gene_set_ids is not None:
            assert isinstance(gene_set_ids, (list, tuple))
            for id_ in gene_set_ids:
                assert isinstance(id_, str)
        if mat is not None:
            assert isinstance(mat, np.ndarray)
        
        gene_set_db = self.gene_set_db
        A = self._A

        if escore_pval_thresh is None:
            # if no separate E-score p-value threshold is specified, use the
            # p-value threshold (this results in very conservative E-scores)
            logger.warning('Setting the E-score p-value threshold to the '
                           'global significance threshold results in '
                           'conservative E-scores.')
            escore_pval_thresh = pval_thresh

        # test only some terms?
        if gene_set_ids is not None:
            gs_indices = np.int64([self.gene_set_db.index(id_)
                                  for id_ in gene_set_ids])
            gene_sets = [gene_set_db[id_] for id_ in gene_set_ids]
            gene_set_db = GeneSetDB(gene_sets)
            A = A[:, gs_indices]  # not a view!

        # reorder rows in annotation matrix to match the given gene ranking
        # also exclude genes not in the ranking
        gene_indices = np.int64([self.genome.index(g) for g in ranked_genes])

        A = A[gene_indices, :]  # not a view either!

        # determine largest K
        K_lim = np.sum(A[:L, :], axis=0, dtype=np.int64)
        K_rem = np.sum(A[L:, :], axis=0, dtype=np.int64)
        K = K_lim + K_rem
        K_max = np.amax(K)

        # prepare matrix for XL-mHG p-value calculation
        p, m = A.shape
        if mat is None:
            mat = np.empty((K_max+1, p+1), dtype=np.longdouble)

        # find enriched GO terms
        logger.info('Testing %d gene sets for enrichment...', m)
        logger.debug('(N=%d, X_frac=%.2f, X_min=%d, L=%d; K_max=%d)',
                     len(ranked_genes), X_frac, X_min, L, K_max)

        enriched = []
        tested = 0  # number of tests conducted
        N, m = A.shape
        for j in range(m):
            # determine gene set-specific value for X (based on K[j])
            X = max(X_min, int(ceil(X_frac*float(K[j]))))

            # determine significance of gene set enrichment using XL-mHG test
            # (only if there are at least X gene set genes in the list)
            if K[j] >= X:
                tested += 1

                # we only need to perform the XL-mHG test if there are enough
                # gene set genes at or above L'th rank (otherwise, pval = 1.0)
                if K_lim[j] >= X:

                    v = np.ascontiguousarray(A[:, j])  # copy
                    stat, n_star, pval = xlmhg_test(
                        v, X, L, K=int(K[j]), table=mat
                    )

                    # check if gene set is significantly enriched
                    if pval <= pval_thresh:
                        # generate GSEResult
                        sel = np.nonzero(A[:, j])[0]  # indices of all the 1's
                        # k_n = np.sum(sel < n)
                        sel_genes = [ranked_genes[i] for i in sel]
                        result = GSEResult(stat, n_star, pval, N, X, L,
                                           gene_set_db[j], sel, sel_genes,
                                           escore_pval_thresh)
                        enriched.append(result)

        # report results
        q = len(enriched)
        ignored = m-tested
        if ignored > 0:
            logger.debug('%d / %d gene sets (%.1f%%) had less than X genes '
                         'annotated with them and were ignored.',
                         ignored, m, 100*(ignored/float(m)))

        logger.info('%d / %d gene sets were found to be significantly '
                    'enriched (p-value <= %.1e).', q, m, pval_thresh)

        return enriched
