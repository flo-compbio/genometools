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

"""Module containing the `GeneSetEnrichmentAnalysis` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
from math import ceil
import copy
from collections import Iterable

import numpy as np
from scipy.stats import hypergeom

import xlmhg

# from ..basic import GeneSet, GeneSetDB
from ..basic import GeneSetDB
from ..expression import ExpGenome
from . import StaticGSEResult, RankBasedGSEResult

logger = logging.getLogger(__name__)


class GeneSetEnrichmentAnalysis(object):
    """Test a set of genes or a ranked list of genes for gene set enrichment.

    Parameters
    ----------
    genome: `ExpGenome` object
        See :attr:`_genome` attribute.
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

    Once the class has been initialized, the function `get_static_enrichment`
    can be used to test a set of genes for gene set enrichment, and the
    function `get_rank_based_enrichment` can be used to test a ranked list of
    genes for gene set enrichment.
    """

    def __init__(self, genome, gene_set_db):

        assert isinstance(genome, ExpGenome)
        assert isinstance(gene_set_db, GeneSetDB)

        self._genome = genome
        self._gene_set_db = gene_set_db

        # generate annotation matrix by going over all gene sets
        logger.info('Generating gene-by-gene set membership matrix...')
        gene_memberships = np.zeros((len(genome), gene_set_db.n),
                                    dtype=np.uint8)
        for j, gs in enumerate(self._gene_set_db.gene_sets):
            for g in gs.genes:
                try:
                    idx = self._genome.index(g)
                except ValueError:
                    pass
                else:
                    gene_memberships[idx, j] = 1
        self._gene_memberships = gene_memberships

    def __repr__(self):
        return '<%s object (_genome=%s; gene_set_db=%s)>' \
               % (self.__class__.__name__,
                  repr(self._genome), repr(self._gene_set_db))

    def __str__(self):
        return '<%s object (%d genes in _genome; %d gene sets)>' \
               % (self.__class__.__name__,
                  len(self._genome), len(self._gene_set_db))

    @property
    def genes(self):
        return self._genome.genes

    @property
    def genome(self):
        return copy.deepcopy(self._genome)

    def get_static_enrichment(
            self, genes, pval_thresh, adjust_pval_thresh=True, X=3,
            gene_set_ids=None):
        """Find enriched gene sets in a set of genes."""
        assert isinstance(genes, set)
        assert isinstance(pval_thresh, (float, np.float))
        assert isinstance(X, (int, np.integer))
        if gene_set_ids is not None:
            assert isinstance(gene_set_ids, Iterable)

        gene_set_db = self._gene_set_db
        gene_sets = gene_set_db.gene_sets
        gene_memberships = self._gene_memberships
        sorted_genes = sorted(genes)

        # test only some terms?
        if gene_set_ids is not None:
            gs_indices = np.int64([self._gene_set_db.index(id_)
                                   for id_ in gene_set_ids])
            gene_sets = [gene_set_db[id_] for id_ in gene_set_ids]
            # gene_set_db = GeneSetDB(gene_sets)
            gene_memberships = gene_memberships[:, gs_indices]  # not a view!

        # determine K's
        K_vec = np.sum(gene_memberships, axis=0, dtype=np.int64)

        # exclude terms with too few genes
        sel = np.nonzero(K_vec >= X)[0]
        K_vec = K_vec[sel]
        gene_sets = [gene_sets[j] for j in sel]
        gene_memberships = gene_memberships[:, sel]

        # determine k's, ignoring unknown genes
        unknown = 0
        sel = []
        filtered_genes = []
        logger.debug('Looking up indices for %d genes...', len(sorted_genes))
        for i, g in enumerate(sorted_genes):
            assert isinstance(g, str)
            try:
                idx = self._genome.index(g)
            except ValueError:
                unknown += 1
            else:
                sel.append(idx)
                filtered_genes.append(g)

        sel = np.int64(sel)
        gene_indices = np.int64(sel)
        # gene_memberships = gene_memberships[sel, :]
        k_vec = np.sum(gene_memberships[sel, :], axis=0, dtype=np.int64)
        if unknown > 0:
            logger.warn('%d / %d unknown genes (%.1f %%), will be ignored.',
                        unknown, len(genes),
                        100 * (unknown / float(len(genes))))

        # determine n and N
        n = len(filtered_genes)
        N, m = gene_memberships.shape
        logger.info('Conducting %d tests.', m)

        # correct p-value threshold, if specified
        final_pval_thresh = pval_thresh
        if adjust_pval_thresh:
            final_pval_thresh /= float(m)
            logger.info('Using Bonferroni-corrected p-value threshold: %.1e',
                        final_pval_thresh)

        # calculate p-values and get significantly enriched gene sets
        enriched = []

        logger.debug('N=%d, n=%d', N, n)
        sys.stdout.flush()
        for j in range(m):
            pval = hypergeom.sf(k_vec[j] - 1, N, K_vec[j], n)
            if pval <= final_pval_thresh:
                # found significant enrichment
                # sel_genes = [filtered_genes[i] for i in np.nonzero(gene_memberships[:, j])[0]]
                sel_genes = [self._genome[i] for i in
                             np.nonzero(gene_memberships[gene_indices, j])[0]]
                enriched.append(
                    StaticGSEResult(N, gene_sets[j], n, set(sel_genes), pval))

        return enriched

    def get_rank_based_enrichment(
            self, ranked_genes, pval_thresh, X_frac, X_min, L,
            adjust_pval_thresh=True, escore_pval_thresh=None,
            gene_set_ids=None, table=None):
        """Find enriched gene sets in a ranked list of genes.

        This function uses the XL-mHG test to identify enriched gene sets.

        This function also calculates XL-mHG E-scores for the enriched gene
        sets, using ``escore_pval_thresh`` as the p-value threshold "psi".
        """
        if isinstance(X_frac, (int, np.integer)):
            X_frac = float(X_frac)

        # type checks
        assert isinstance(ranked_genes, Iterable)
        assert isinstance(pval_thresh, (float, np.float))
        assert isinstance(X_frac, (float, np.float))
        assert isinstance(X_min, (int, np.integer))
        assert isinstance(L, (int, np.integer))
        assert isinstance(adjust_pval_thresh, bool)

        if escore_pval_thresh is not None:
            assert isinstance(escore_pval_thresh, (float, np.float))
        if gene_set_ids is not None:
            assert isinstance(gene_set_ids, Iterable)
        if table is not None:
            assert isinstance(table, np.ndarray) and \
                   np.issubdtype(table.dtype, np.longdouble)

        gene_set_db = self._gene_set_db
        gene_memberships = self._gene_memberships

        # postpone this
        if escore_pval_thresh is None:
            # if no separate E-score p-value threshold is specified, use the
            # p-value threshold (this results in very conservative E-scores)
            logger.warning('No E-score p-value threshold supplied Setting the '
                           'E-score '
                           'p-value threshold to the '
                           'global significance threshold results in '
                           'conservative E-scores.')

        # test only some terms?
        if gene_set_ids is not None:
            gs_indices = np.int64([self._gene_set_db.index(id_)
                                   for id_ in gene_set_ids])
            gene_sets = [gene_set_db[id_] for id_ in gene_set_ids]
            gene_set_db = GeneSetDB(gene_sets)
            gene_memberships = gene_memberships[:, gs_indices]  # not a view!

        # reorder rows in annotation matrix to match the given gene ranking
        # also exclude genes not in the ranking
        unknown = 0
        L_adj = L
        sel = []
        filtered_genes = []
        logger.debug('Looking up indices for %d genes...' % len(ranked_genes))
        for i, g in enumerate(ranked_genes):
            assert isinstance(g, str)
            try:
                idx = self._genome.index(g)
            except ValueError:
                unknown += 1
                # adjust L if the gene was above the original L cutoff
                if i < L:
                    L_adj -= 1
            else:
                sel.append(idx)
                filtered_genes.append(g)
        sel = np.int64(sel)
        #gene_indices = np.int64(sel)
        logger.debug('Adjusted L: %d', L_adj)

        # the following also copies the data (not a view)
        gene_memberships = gene_memberships[sel, :]
        N, m = gene_memberships.shape
        if unknown > 0:
            # Some genes in the ranked list were unknown (i.e., not present in
            # the specified genome).
            logger.warn('%d / %d unknown genes (%.1f %%), will be ignored.',
                        unknown, len(genes),
                        100 * (unknown / float(len(genes))))


        # Determine the number of gene set genes above the L'th cutoff,
        # for all gene sets. This quantity is useful, because we don't need
        # to perform any more work for gene sets that have less than X genes
        # above the cutoff.
        k_above_L = np.sum(gene_memberships[:L_adj, :], axis=0, dtype=np.int64)

        # Determine the number of genes below the L'th cutoff, for all gene
        # sets.
        k_below_L = np.sum(gene_memberships[L_adj:, :], axis=0, dtype=np.int64)

        # Determine the total number K of genes in each gene set that are
        # present in the ranked list (this is equal to k_above_L + k_below_L)
        K_vec = k_above_L + k_below_L

        # Determine the largest K across all gene sets.
        K_max = np.amax(K_vec)

        # Determine X for all gene sets.
        X = np.amax(
            np.c_[np.tile(X_min, m), np.int64(np.ceil(X_frac * K_vec))],
            axis=1)

        # Determine the number of tests (we do not conduct a test if the
        # total number of gene set genes in the ranked list is below X).
        num_tests = np.sum(K_vec-X >= 0)
        logger.info('Conducting %d tests.', num_tests)

        # determine Bonferroni-corrected p-value, if desired
        final_pval_thresh = pval_thresh
        if adjust_pval_thresh and num_tests > 0:
            final_pval_thresh /= float(num_tests)
            logger.info('Using Bonferroni-corrected p-value threshold: %.1e',
                        final_pval_thresh)


        # Prepare the matrix that holds the dynamic programming table for
        # the calculation of the XL-mHG p-value.
        if table is None:
            table = np.empty((K_max+1, N+1), dtype=np.longdouble)
        else:
            if table.shape[0] < K_max+1 or table.shape[1] < N+1:
                raise ValueError(
                    'The supplied array is too small (%d x %d) to hold the '
                    'entire dynamic programming table. The required size is'
                    '%d x %d (rows x columns).'
                    % (table.shape[0], table.shape[1], K_max+1, N+1))

        # find enriched GO terms
        logger.info('Testing %d gene sets for enrichment...', m)
        logger.debug('(N=%d, X_frac=%.2f, X_min=%d, L=%d; K_max=%d)',
                     len(ranked_genes), X_frac, X_min, L, K_max)

        enriched = []
        num_tests = 0  # number of tests conducted
        for j in range(m):
            # determine gene set-specific value for X
            X = max(X_min, int(ceil(X_frac * float(K_vec[j]))))

            # Determine significance of gene set enrichment using the XL-mHG
            # test (only if there are at least X gene set genes in the list).
            if K_vec[j] >= X:
                num_tests += 1

                # We only need to perform the XL-mHG test if there are enough
                # gene set genes above the L'th cutoff (otherwise, pval = 1.0).
                if k_above_L[j] >= X:
                    # perform test

                    # Determine the ranks of the gene set genes in the
                    # ranked list.
                    indices = np.uint16(np.nonzero(gene_memberships[:, j])[0])
                    res = xlmhg.get_xlmhg_test_result(
                        N, indices, X, L, table=table)

                    # check if gene set is significantly enriched
                    if res.pval <= pval_thresh:
                        # generate RankedGSEResult
                        ind_genes = [ranked_genes[i] for i in indices]
                        gse_result = RankBasedGSEResult(
                            gene_set_db[j], N, indices, ind_genes,
                            X, L, res.stat, res.cutoff, res.pval,
                            escore_pval_thresh=escore_pval_thresh
                        )
                        enriched.append(gse_result)

        # report results
        q = len(enriched)
        ignored = m - num_tests
        if ignored > 0:
            logger.debug('%d / %d gene sets (%.1f%%) had less than X genes '
                         'annotated with them and were ignored.',
                         ignored, m, 100 * (ignored / float(m)))

        logger.info('%d / %d gene sets were found to be significantly '
                    'enriched (p-value <= %.1e).', q, m, pval_thresh)

        return enriched